
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
		#if there are no external annots, then just simply concat files
		if [[ -z "${external_annot}" ]]; then
			cat <(head -n 1 ${files[0]}) <(awk 'FNR>1' ${sep=" " files}) | bgzip > annotated_variants.gz
			tabix -S 2 -b 3 -e 3 -s 1 annotated_variants.gz
		else
			#else sort VEP and annotation files in similar order, then join
			#sort VEP
			vep_varcol=$(zcat ${external_annot}|head -n1|tr "\t" "\n"|grep -n "variant"|cut -f 1)
			vep_conscol=$(zcat ${external_annot}|head -n1|tr "\t" "\n"|grep -n "most_severe"|cut -f 1)
			vep_genecol=$(zcat ${external_annot}|head -n1|tr "\t" "\n"|grep -n "gene_most_severe"|cut -f 1)
			cat <(zcat ${external_annot}|head -n1|cut -f $vep_varcol,$vep_conscol,$vep_genecol) <(zcat ${external_annot}|cut -f $vep_varcol,$vep_conscol,$vep_genecol|sort -V -k 1,1)|bgzip > sorted_vep.gz
			#sort annotation
			cat <(head -n 1 ${files[0]}) <(awk 'FNR>1' ${sep=" " files}) | bgzip > annotated_variants_1.gz
			cat <(zcat annotated_variants_1.gz|head -n1)  <(zcat annotated_variants_1.gz|tail -n+2|sort -V -k 1,1 -T ./) |bgzip  > sorted_annotated.gz
			#join
			join -t $'\t' -1 1 -2 1 -e "NA" --header -a 1 <(zcat sorted_annotated.gz) <(zcat sorted_vep.gz) |bgzip  > annotated_variants.gz
			tabix -b 3 -e 3 -s 2 annotated_variants.gz
		fi

	>>>



	output {
		File out="annotated_variants.gz"
		File index="annotated_variants.gz.tbi"
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

}
