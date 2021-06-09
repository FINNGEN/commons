
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

		if [[ -z "${external_annot}" ]]; then
			cat <(head -n 1 ${files[0]}) <(awk 'FNR>1' ${sep=" " files}) | bgzip > annotated_variants.gz
			tabix -S 2 -b 3 -e 3 -s 1 annotated_variants.gz
		else
            n_join1=$(( `head -n 1 ${files[0]}| awk 'BEGIN{ FS=OFS="\t"} { print NF};'` + `gunzip -c "${external_annot}" | head -n 1 | awk 'BEGIN{ FS=OFS="\t"} { print NF-1};'`))
            join -a 1 -t $'\t' -1 2 -2 1 \
            <(cat <(head -n 1 ${files[0]} | awk 'BEGIN{ FS=OFS="\t"} { print "#"$0}' ) \
            <(awk 'BEGIN{ FS=OFS="\t"} FNR>1' ${sep=" " files} | awk 'BEGIN{ FS=OFS="\t"} { gsub("_", ":", $1); sub("chr", "", $1); sub("X", 23, $1); sub("Y", 24, $1); sub("MT", 25, $1); sub("M", 25, $1); gsub("chr", "", $2); sub("X", 23, $2); sub("Y", 24, $2); sub("MT", 25, $2); sub("M", 25, $2); print $0 }') \
            | nl -nln | sort -b -k 2,2 ) \
            <(gunzip -c "${external_annot}" | awk 'BEGIN{ FS=OFS="\t"} NR==1{ print $0  } \
            NR>1{ gsub("_", ":", $1); sub("chr", "", $1); sub("X", 23, $1); sub("Y", 24, $1); sub("MT", 25, $1); sub("M", 25, $1); gsub("chr", "", $2); sub("X", 23, $2); sub("Y", 24, $2); sub("MT", 25, $2); sub("M", 25, $2); print $0}' \
            | sort -b -k 1,1  ) | sort -g -k 2,2 | cut -f 1,3- \
            | awk -v nf=$n_join1 'BEGIN{ FS=OFS="\t"} { printf $0; if(NF<nf) { for(i=NF+1;i<=nf;i++) { printf "\tNA"  };  } printf "\n" }' | bgzip > annotated_variants.gz
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

                out.write( "variant\tchr\tpos\t" + "\t".join(items) + "\n" )

                for line in infile:
                    dat = line.split("\t")[0:8]
                    var = dat[0] + "_" + dat[1] + "_" + dat[3] + "_" + dat[4]
                    info = { d[0]:(d[1] if len(d)>1 else "") for d in list(map(lambda x: x.split("="), dat[7].split(";")) ) }

                    def getdat( elem, info):
                        if elem["Type"]=="Flag":
                            return "1" if elem["ID"] in info else "0"
                        else:
                            return info[elem["ID"]] if elem["ID"] in info else "NA"

                    out.write( var + "\t" + dat[0] + "\t" + dat[1] + "\t" + "\t".join([  getdat(elem, info) for elem in fields  ] ) + "\n")

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
