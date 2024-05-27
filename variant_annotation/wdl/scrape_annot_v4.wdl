
task join_annot {

	Array[File] files
	File? external_annot
	Int local_disk=500
  	String docker

	runtime {
		docker: "${docker}"
		memory:"8G"
		disks: "local-disk ${local_disk} HDD"
		cpu:"1"
		zones: "europe-west1-b europe-west1-c europe-west1-d"
		preemptible: 0
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

        #helper function
        body ()
        {
            IFS= read -r header;
            printf '%s\n' "$header";
            "$@"
        }

        awk 'BEGIN{FS=OFS="\t"} NR==1 { for(i=1;i<=NF;i++) sub(/_chr1$/, "", $i); print } FNR>1' ${sep=" " files} | bgzip > annotated_variants.gz.tmp

		#if there are no external annots, then just simply concat files
		if [[ ! -z "${external_annot}" ]]; then
			#headercheck: check that the column is in header, else abort with error 1
			vep_varcol=$(zcat ${external_annot} | head -n1 | awk -v value="variant" -F "\t" -f headercheck.awk)
			vep_conscol=$(zcat ${external_annot}| head -n1 | awk -v value="most_severe" -F "\t" -f headercheck.awk )
			vep_genecol=$(zcat ${external_annot}| head -n1 | awk -v value="gene_most_severe" -F "\t" -f headercheck.awk)

			zcat ${external_annot} | cut -f $vep_varcol,$vep_conscol,$vep_genecol | body sort -V -k 1,1 -T ./ | bgzip > sorted_vep.gz
			#sort annotation
			zcat annotated_variants.gz.tmp | body sort -V -k 1,1 -T ./ | bgzip  > sorted_annotated.gz
			#join
			join --nocheck-order -t $'\t' -1 1 -2 1 -e "NA" --header -a 1 <(zcat sorted_annotated.gz) <(zcat sorted_vep.gz) | bgzip > annotated_variants.gz
		fi

		tabix -b 3 -e 3 -s 2 annotated_variants.gz

	>>>

	output {
		File out="annotated_variants.gz"
		File index="annotated_variants.gz.tbi"
	}

}


task small {
	File annotation
	String wanted_columns
	Int local_disk = ceil(size(annotation, "G")) + 10
	String docker
	runtime {
		docker: "${docker}"
		memory:"8G"
		disks: "local-disk ${local_disk} HDD"
		cpu:"1"
		zones: "europe-west1-b europe-west1-c europe-west1-d"
		preemptible: 0
		noAddress: true
	}

	command <<<
		set -eux
		cat > cols.awk <<'AWK'
		BEGIN {
		FS=OFS="\t"
		len_head=split("${wanted_columns}", wanted, " ")
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
	Int local_disk = ceil(size(vcf, "G")) + 100
	String docker
	String outfile=basename(vcf) + ".annot.gz"

	runtime {
		docker: "${docker}"
		memory:"2G"
		disks: "local-disk ${local_disk} HDD"
		cpu:"1"
		zones: "europe-west1-b europe-west1-c europe-west1-d"
		preemptible: 0
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

    Int local_disk = (ceil(size(annotation_file, "G")) + ceil(size(ref_file, "G"))) * 3

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
            #variants that are on same c:p
            ref_vars = []
            for line in f:
                #get variant values
                line = line.strip()
                s = line.split('\t')
                chrom = int(s[h_idx['chr']])
                pos = int(s[h_idx['pos']])
                ref = s[h_idx['ref']]
                alt = s[h_idx['alt']]
                #check if ref vars c:p has moved forward
                if ref_vars:
                    if int(ref_vars[0][ref_h_idx['#CHROM']]) != chrom or int(ref_vars[0][ref_h_idx['POS']]) != pos:
                        ref_vars = []
                while ref_has_lines and ref_chr < chrom or (ref_chr == chrom and ref_pos < pos):
                    ref_line = fp_ref.readline().rstrip('\r\n').split('\t')
                    try:
                        ref_chr = int(ref_line[ref_h_idx['#CHROM']])
                        ref_pos = int(ref_line[ref_h_idx['POS']])
                    except ValueError:
                        ref_has_lines = False
                while ref_has_lines and int(ref_chr) == chrom and ref_pos == pos:
                    ref_vars.append(ref_line)
                    ref_line = fp_ref.readline().strip().split('\t')
                    try:
                        ref_chr = int(ref_line[ref_h_idx['#CHROM']])
                        ref_pos = int(ref_line[ref_h_idx['POS']])
                    except ValueError:
                        ref_has_lines = False

                rsid = 'NA'
                #if chrom and pos match and there are items in the ref vars, we try to match.
                # If chrom and pos don't match, we reset ref vars
                if ref_vars:
                    if (chrom == int(ref_vars[0][ref_h_idx["#CHROM"]]) ) and (pos == int(ref_vars[0][ref_h_idx["POS"]]) ):
                        for r in ref_vars:
                            if r[ref_h_idx['REF']] == ref and alt in r[ref_h_idx['ALT']].split(','):
                                rsid = r[ref_h_idx['ID']]
                                break
                    else:
                        #chrom and pos have advanced further than ref vars, clear ref vars
                        ref_vars = []
                
                print(line + '\t' + rsid)

        EOF

        echo "`date` done"

    >>>

    output {
        File rsid_out = "annotation.gz"
    }

    runtime {
        docker: "${docker}"
        cpu: "2"
        memory: "8 GB"
        disks: "local-disk ${local_disk} HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 0
        noAddress: true
    }

}

task add_gnomad {
	File annotation_file

	File gnomad_genomes
	File gnomad_exomes
    Array[String] exome_cols
    Array[String] genome_cols

    Int local_disk = (ceil(size(annotation_file, "G")) + ceil(size(gnomad_genomes, "G")) + ceil(size(gnomad_exomes, "G"))) * 3

	String docker

	command <<<
        set -euxo pipefail

        echo "`date` Adding gnomad annotations"
        
        python3 <<EOF | bgzip > annotation_gnomad.gz
        
        import gzip
        class Reference:
            def __init__(self, fp):
                self.fp = fp
                self.chrom = 1
                self.pos = 0
                self.has_lines=True
                self.line = self.fp.readline()
                self.header = self.line.rstrip('\r\n').split('\t')
                self.h_idx = {h:i for i,h in enumerate(self.line.rstrip('\r\n').split('\t'))}
                self.has_lines = True
                self.has_b37_coords = 'CHROM_37' in self.header

        def format_chrom(chr: str):
            return int(chr.replace("chr", "").replace("X", "23").replace("Y", "24").replace("M", "25").replace("MT", "25"))

        fp_genome = gzip.open('${gnomad_genomes}', 'rt')
        fp_exome = gzip.open('${gnomad_exomes}', 'rt')

        gen_ref = Reference(fp_genome)
        exo_ref = Reference(fp_exome)

        exome_cols = ["${sep='\",\"' exome_cols}"]
        genome_cols = ["${sep='\",\"' genome_cols}"]

        assert all([col in gen_ref.header for col in genome_cols]), "All given exome columns not found from exome header"
        assert all([col in exo_ref.header for col in exome_cols]), "All given exome columns not found from exome header"

        exome_cols_headernames = [f"EXOME_{val}" for val in exome_cols]
        genome_cols_headernames = [f"GENOME_{val}" for val in genome_cols]

        exome_rename_dict = {a:b for a,b in zip(exome_cols,exome_cols_headernames)}
        genome_rename_dict = {a:b for a,b in zip(genome_cols,genome_cols_headernames)}

        with gzip.open('${annotation_file}', 'rt') as f:
            header = f.readline().strip()
            h_idx = {h:i for i,h in enumerate(header.split('\t'))}
            ex_head = '\t'.join(exome_cols_headernames)
            gen_head = '\t'.join(genome_cols_headernames)
            out_header = f"{header}\tb37_coord\t{ex_head}\t{gen_head}" if gen_ref.has_b37_coords or exo_ref.has_b37_coords else f"{header}\t{ex_head}\t{gen_head}"
            print(out_header)

            #create vars for each of the gen/exomes
            gen_vars = []
            exo_vars=[]
            for line in f:
                oline = line.rstrip("\r\n")
                s = oline.split('\t')
                chrom = format_chrom(s[h_idx['chr']])
                pos = int(s[h_idx['pos']])
                ref = s[h_idx['ref']]
                alt = s[h_idx['alt']]
                #check if gen and exo vars have different c:p than variant. If yes, empty them.
                if gen_vars:
                    if format_chrom(gen_vars[0][gen_ref.h_idx['#CHROM']]) != chrom or int(gen_vars[0][gen_ref.h_idx['POS']]) != pos:
                        gen_vars = []

                if exo_vars:
                    if format_chrom(exo_vars[0][exo_ref.h_idx['#CHROM']]) != chrom or int(exo_vars[0][exo_ref.h_idx['POS']]) != pos:
                        exo_vars = []

                #genome
                while gen_ref.has_lines and gen_ref.chrom < chrom or ( gen_ref.chrom == chrom and gen_ref.pos < pos):
                    gen_ref.line = gen_ref.fp.readline().rstrip('\r\n').split('\t')
                    try:
                        gen_ref.pos = int(gen_ref.line[gen_ref.h_idx["POS"]])
                        gen_ref.chrom = format_chrom(gen_ref.line[gen_ref.h_idx["#CHROM"]])
                    except ValueError:
                        gen_ref.has_lines = False
                    
                #exome
                while exo_ref.has_lines and exo_ref.chrom < chrom or (exo_ref.chrom == chrom and exo_ref.pos < pos):
                    exo_ref.line = exo_ref.fp.readline().rstrip('\r\n').split('\t')
                    try:
                        exo_ref.pos = int(exo_ref.line[exo_ref.h_idx["POS"]])
                        exo_ref.chrom = format_chrom(exo_ref.line[exo_ref.h_idx["#CHROM"]])
                    except ValueError:
                        exo_ref.has_lines = False
                
                #drain all lines with same chrom & pos to gen_vars
                while gen_ref.has_lines and ( gen_ref.chrom == chrom and gen_ref.pos == pos ):
                    gen_vars.append(gen_ref.line)
                    gen_ref.line = gen_ref.fp.readline().rstrip("\r\n").split('\t')
                    try:
                        gen_ref.chrom = format_chrom(gen_ref.line[gen_ref.h_idx['#CHROM']])
                        gen_ref.pos = int(gen_ref.line[gen_ref.h_idx['POS']])
                    except ValueError:
                        gen_ref.has_lines = False
                
                while exo_ref.has_lines and ( exo_ref.chrom == chrom and exo_ref.pos == pos ):
                    exo_vars.append(exo_ref.line)
                    exo_ref.line = exo_ref.fp.readline().rstrip("\r\n").split('\t')
                    try:
                        exo_ref.chrom = format_chrom(exo_ref.line[exo_ref.h_idx['#CHROM']])
                        exo_ref.pos = int(exo_ref.line[exo_ref.h_idx['POS']])
                    except ValueError:
                        exo_ref.has_lines = False

                #if the chrom and pos are the same, we try to match the variants. Else, we reset the respective {gen,exo}_vars
                gen_values = ["NA"]*len(genome_cols)
                if gen_vars:
                    if chrom == format_chrom(gen_vars[0][gen_ref.h_idx['#CHROM']]) and pos == int(gen_vars[0][gen_ref.h_idx['POS']]):
                        for genvar in gen_vars:
                            match = (ref == genvar[gen_ref.h_idx['REF']]) and (alt == genvar[gen_ref.h_idx['ALT']])
                            if match:
                                gen_values = [genvar[gen_ref.h_idx[col]] for col in genome_cols]
                                break
                    else: #reset
                        gen_vars = []

                exo_values = ["NA"]*len(exome_cols)
                if exo_vars:
                    if chrom == format_chrom(exo_vars[0][exo_ref.h_idx['#CHROM']]) and pos == int(exo_vars[0][exo_ref.h_idx['POS']]):
                        for exvar in exo_vars:
                            match = (ref == exvar[exo_ref.h_idx['REF']]) and (alt == exvar[exo_ref.h_idx['ALT']])
                            if match:
                                exo_values = [exvar[exo_ref.h_idx[col]] for col in exome_cols]
                                break
                    else: #reset
                        exo_vars=[]
                
                exvals = '\t'.join(exo_values)
                gevals = '\t'.join(gen_values)
                #b37 coord
                b37_coord = "NA"
                if exo_ref.has_b37_coords:
                    for exo in exo_vars:
                        match = (chrom == format_chrom(exo[exo_ref.h_idx['#CHROM']]) ) and (pos == int(exo[exo_ref.h_idx['POS']]) ) and (ref == exo[exo_ref.h_idx['REF']]) and (alt == exo[exo_ref.h_idx['ALT']])
                        if match:
                            #CHROM_37	POS_37	REF_37	ALT_37
                            b37_coord = f"{exo[exo_ref.h_idx['CHROM_37']]}:{exo[exo_ref.h_idx['POS_37']]}:{exo[exo_ref.h_idx['REF_37']]}:{exo[exo_ref.h_idx['ALT_37']]}"
                            break
                if gen_ref.has_b37_coords:
                    for geno in gen_vars:
                        if b37_coord != "NA":
                            break
                        match = (chrom == format_chrom(geno[gen_ref.h_idx['#CHROM']]) ) and (pos == int(geno[gen_ref.h_idx['POS']]) ) and (ref == geno[gen_ref.h_idx['REF']]) and (alt == geno[gen_ref.h_idx['ALT']])
                        if match:
                            b37_coord = f"{geno[gen_ref.h_idx['CHROM_37']]}:{geno[gen_ref.h_idx['POS_37']]}:{geno[gen_ref.h_idx['REF_37']]}:{geno[gen_ref.h_idx['ALT_37']]}"
                            break
                outline = f"{oline}\t{b37_coord}\t{exvals}\t{gevals}" if gen_ref.has_b37_coords or exo_ref.has_b37_coords else f"{oline}\t{exvals}\t{gevals}"
                print(outline)
        gen_ref.fp.close()
        exo_ref.fp.close()
        EOF

        echo "`date` tabixing"
        tabix -s 2 -b 3 -e 3 annotation_gnomad.gz
	>>>
	output {
		File gnomad_joined_out = "annotation_gnomad.gz"
		File gnomad_joined_out_tbi = "annotation_gnomad.gz.tbi"
	}
	runtime {
        docker: "${docker}"
        cpu: "2"
        memory: "8 GB"
        disks: "local-disk ${local_disk} HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 0
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

    call add_gnomad {
        input: docker=docker,
        annotation_file = add_rsids.rsid_out,
    }

	call small {
		input: annotation = add_gnomad.gnomad_joined_out,
		docker=docker
	}

	output {
		File annotation = add_gnomad.gnomad_joined_out
		File annotatation_index = add_gnomad.gnomad_joined_out_tbi
		File annotation_small = small.small_out
		File annotation_small_index = small.small_index
	}

}
