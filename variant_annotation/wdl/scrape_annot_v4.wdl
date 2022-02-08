
task join_annot {

	Array[File] files
	File? external_annot
	Int local_disk=500
  	String docker

	runtime {
		docker: "${docker}"
		memory:"8G"
		disks: "local-disk ${local_disk} SSD"
		cpu:"1"
		zones: "europe-west1-b europe-west1-c europe-west1-d"
		preemptible: 1
		noAddress: true
	}

	command <<<
		set -eux
		#save headercheck to file
		cat > headercheck.awk <<'CODE'
		{
		for(i = 1; i <= NF; i++) {
			h[$i] = i
			if($i==value){print i}
		}
		exists=value in h
		if(!exists)
		{
			system("echo ERROR header did not contain column "value"  >&2")
			exit 1
		}
		exit 0
		}
		CODE

		#if there are no external annots, then just simply concat files
		if [[ -z "${external_annot}" ]]; then
			cat <(head -n 1 ${files[0]}) <(awk 'FNR>1' ${sep=" " files}) | bgzip > annotated_variants.gz
		else
			#else sort VEP and annotation files in similar order, then join
			
			#headercheck: check that the column is in header, else abort with error 1
			vep_varcol=$(zcat ${external_annot} |head -n1|awk -v value="variant" -F "\t" -f headercheck.awk)
			vep_conscol=$(zcat ${external_annot}|head -n1|awk -v value="most_severe" -F "\t" -f headercheck.awk )
			vep_genecol=$(zcat ${external_annot}|head -n1|awk -v value="gene_most_severe" -F "\t" -f headercheck.awk)
			cat <(zcat ${external_annot}|head -n1|cut -f $vep_varcol,$vep_conscol,$vep_genecol) <(zcat ${external_annot}|cut -f $vep_varcol,$vep_conscol,$vep_genecol|sort -V -k 1,1)|bgzip > sorted_vep.gz
			#sort annotation
			cat <(head -n 1 ${files[0]}) <(awk 'FNR>1' ${sep=" " files}) | bgzip > annotated_variants_1.gz
			cat <(zcat annotated_variants_1.gz|head -n1)  <(zcat annotated_variants_1.gz|tail -n+2|sort -V -k 1,1 -T ./) |bgzip  > sorted_annotated.gz
			#join
			join -t $'\t' -1 1 -2 1 -e "NA" --header -a 1 <(zcat sorted_annotated.gz) <(zcat sorted_vep.gz) |bgzip  > annotated_variants.gz
		fi

		tabix -b 3 -e 3 -s 2 annotated_variants.gz

	>>>



	output {
		File out="annotated_variants.gz"
		File index="annotated_variants.gz.tbi"
	}


}


task small{
	File annotation
	Int local_disk=100
	String docker
	runtime {
		docker: "${docker}"
		memory:"8G"
		disks: "local-disk ${local_disk} HDD"
		cpu:"1"
		zones: "europe-west1-b europe-west1-c europe-west1-d"
		preemptible: 1
		noAddress: true
	}

	command <<<
		set -eux
		cat > cols.awk <<'AWK'
		BEGIN {
		FS=OFS="\t"
		len_head=split("#variant chr pos ref alt INFO AF AC_Het AC_Hom most_severe gene_most_severe rsid", wanted, " ")
		}
		NR==1 {
		for(i=1;i<=NF;i++) h[$i]=i
		for(i in wanted) if (!(wanted[i] in h)) {print wanted[i]" expected in header">>"/dev/stderr"; exit 1}
		for(i=1;i<=len_head;i++) { printf $h[wanted[i]]"\t" } printf "index\n"
		}
		NR >1 {
		for(i=1;i<=len_head;i++) { printf $h[wanted[i]]"\t" } printf NR-2"\n"
		}
		AWK

		zcat ${annotation}|awk -f cols.awk |bgzip  > annotated_variants.small.gz
		tabix -b 3 -e 3 -s 2 annotated_variants.small.gz
	>>>

	output {
		File small_out="annotated_variants.small.gz"
		File small_index="annotated_variants.small.gz.tbi"
	}
}

task extract {

	File vcf
	Int local_disk=200
	String docker
	String outfile=basename(vcf) + ".annot.gz"

	runtime {
		docker: "${docker}"
		memory:"2G"
		disks: "local-disk ${local_disk} SSD"
		cpu:"1"
		zones: "europe-west1-b europe-west1-c europe-west1-d"
		preemptible: 1
		noAddress: true
	}

	command <<<

        python3 <<EOF

        import gzip
        import re

        def read_info_fields(file):
            fields = []
            for line in file:
                if line.startswith("#CHROM"):
                    break
                if line.startswith("##INFO=<"):
                    line = line.rstrip("\n").replace(', ', '-') # next line splits on commas also within Description so avoid breaking
                    dat = { elem[0]:elem[1] for elem in map(lambda x: x.split("="), re.sub("^##INFO=<|>$","",line).split(",")) }
                    fields.append(dat)

            return fields

        with open('${outfile}', 'w') as out:
            with gzip.open('${vcf}',mode='rt') as infile:
                fields = read_info_fields(infile)
                fields.sort( key=lambda x: x["ID"]  )

                items = map(lambda x: x["ID"], fields )

                out.write( "#variant\tchr\tpos\tref\talt\t" + "\t".join(items) + "\n" )
                chromsanitize = lambda c: c.replace("chr","").replace("X","23").replace("Y","24").replace("MT","25").replace("M","25")
                for line in infile:
                    dat = line.split("\t")[0:8]
                    var = chromsanitize(dat[0]) + ":" + dat[1] + ":" + dat[3] + ":" + dat[4]
                    info = { d[0]:(d[1] if len(d)>1 else "") for d in list(map(lambda x: x.split("="), dat[7].split(";")) ) }

                    def getdat( elem, info):
                        if elem["Type"]=="Flag":
                            return "1" if elem["ID"] in info else "0"
                        else:
                            return info[elem["ID"]] if elem["ID"] in info else "NA"

                    out.write( var + "\t" + chromsanitize(dat[0]) + "\t" + dat[1] + "\t" + dat[3] + "\t" + dat[4] + "\t" + "\t".join([  getdat(elem, info) for elem in fields  ] ) + "\n")

        EOF

	>>>

	output {
		File out="${outfile}"
	}

}

# Add rsids
task add_rsids {

	File annotation_file

	File ref_file
	String docker


    command <<<

        set -euxo pipefail

        echo "`date` Adding rsids"

        python3 <<EOF | bgzip > annotation.gz

        import gzip

        fp_ref = gzip.open('${ref_file}', 'rt')
        ref_has_lines = True
        ref_chr = 1
        ref_pos = 0
        ref_line = fp_ref.readline()
        while ref_line.startswith("##"):
            ref_line = fp_ref.readline()
        if ref_line.startswith('#'):
            assert ref_line.rstrip('\r\n').split('\t') == '#CHROM POS ID REF ALT QUAL FILTER INFO'.split(), repr(ref_line)
        ref_h_idx = {h:i for i,h in enumerate(ref_line.rstrip('\r\n').split('\t'))}

        with gzip.open('${annotation_file}', 'rt') as f:
            header = f.readline().strip()
            h_idx = {h:i for i,h in enumerate(header.split('\t'))}
            print(header + '\trsid')
            for line in f:
                line = line.strip()
                s = line.split('\t')
                chr = int(s[h_idx['#CHR']])
                pos = int(s[h_idx['POS']])
                ref = s[h_idx['REF']]
                alt = s[h_idx['ALT']]
                ref_vars = []
                while ref_has_lines and int(ref_chr) < chr or (int(ref_chr) == chr and ref_pos < pos):
                    ref_line = fp_ref.readline().rstrip('\r\n').split('\t')
                    try:
                        ref_chr = ref_line[ref_h_idx['#CHROM']]
                        ref_pos = int(ref_line[ref_h_idx['POS']])
                    except ValueError:
                        ref_has_lines = False
                while ref_has_lines and int(ref_chr) == chr and ref_pos == pos:
                    ref_vars.append(ref_line)
                    ref_line = fp_ref.readline().strip().split('\t')
                    try:
                        ref_chr = ref_line[ref_h_idx['#CHROM']]
                        ref_pos = int(ref_line[ref_h_idx['POS']])
                    except ValueError:
                        ref_has_lines = False

                rsid = 'NA'
                for r in ref_vars:
                    if r[ref_h_idx['REF']] == ref and alt in r[ref_h_idx['ALT']].split(','):
                        rsid = r[ref_h_idx['ID']]
                        break

                print(line + '\t' + rsid)

        EOF

        echo "`date` tabixing"
        tabix -s 1 -b 2 -e 2 annotation.gz
        echo "`date` done"

    >>>

    output {
        File rsid_out = "annotation.gz"
        File rsid_tbi = "annotation.gz.tbi"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk 200 SSD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }

}

workflow scrape_annots {

	File vcfs

	Array[String] files = read_lines(vcfs)

	String docker


	## write readme and specify annot: output will be in the same order as input files.... external annot will be matched by first col. id chr:pos:ref:alt
	## join external annotation already in the extract for speed!!

	File? external_annot

	scatter(file in files) {
		call extract {
			input: vcf=file,
      		docker=docker
		}
	}
	call join_annot {
		input: files=extract.out,
		external_annot=external_annot,
		docker=docker
	}

	call add_rsids {
		input: docker=docker,
		annotation_file = join_annot.out
	}

	call small {
		input: annotation = add_rsids.rsid_out,
		docker=docker
	}

	output {
		File annotation = add_rsids.rsid_out
		File annotatation_index = add_rsids.rsid_tbi
		File annotation_small = small.small_out
		File annotation_small_index = small.small_index
	}

}
