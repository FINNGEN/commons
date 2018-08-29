
task join_annot {

	Array[File] files
	File? external_annot
	Int local_disk=500
  	String docker
  	String? outputfile

	runtime {
		docker: "${docker}"
		memory:"8G"
		disks: "local-disk ${local_disk} SSD"
		cpu:"1"
	}

	command <<<

		if [[ -z "${external_annot}" ]]; then
			cat <(head -n 1 ${files[0]}) <(awk 'FNR>1' ${sep=" " files}) | bgzip > annotated_variants.gz
			tabix -S 2 -b 3 -e 3 -s 1 annotated_variants.gz
		else
			n_join1=$(( `head -n 1 ${files[0]}| awk 'BEGIN{ FS="\t"} { print NF};'` + `head -n 1 "${external_annot}" | awk 'BEGIN{ FS="\t"} { print NF-1};'`))
			join -a 1 -t $'\t' -1 2 -2 1 <(cat <(head -n 1 ${files[0]} | awk '{ print "#"$0}' ) <(awk 'FNR>1' ${sep=" " files}) | nl -nln | sort -b -k 2,2 ) <( awk 'BEGIN{ FS="\t"} NR==1{ print "#"$0  } NR>1{ print $0}' "${external_annot}" | sort -b -k 1,1  ) | sort -g -k 2,2 | cut -f 1,3- | awk -v nf=$n_join1 'BEGIN{ FS="\t";OFS="\t"} { printf $0; if(NF<nf) { for(i=NF+1;i<=nf;i++) { printf "\tNA"  };  } printf "\n" }' | bgzip > annotated_variants.gz
			tabix -b 3 -e 3 -s 2 annotated_variants.gz
		fi

	    if [[ -n "${outputfile}" ]]; then
	      gsutil cp annotated_variants.gz ${outputfile}
	      gsutil cp annotated_variants.gz.tbi ${outputfile}".tbi"
	    fi

	>>>



	output {
		File out="annotated_variants.gz"
		File index="annotated_variants.gz.tbi"
	}


}

task extract {
	File vcf
	Int local_disk=100
  	String docker

    String outfile=basename(vcf) + ".annot.gz"


	runtime {
		docker: "${docker}"
		memory:"2G"
		disks: "local-disk ${local_disk} SSD"
		cpu:"1"
	}

	command <<<
        scrape_vcf_info.py ${vcf} ${outfile}
	>>>

	output {
		File out="${outfile}"
	}

}

workflow scrape_annots {

	File vcfs

	Array[String] files = read_lines(vcfs)

	String docker
	String memory= "3GB"
	Int cpu = 1

	## write readme and specify annot: output will be in the same order as input files.... external annot will be matched by first col. id chr:pos:ref:alt
	## join external annotation already in the extract for speed!!

	File? external_annot
	String? outputfile

	scatter(file in files) {
		call extract {
			input: vcf=file,
      		docker=docker
		}
	}

	call join_annot {
		input: files=extract.out,
		external_annot=external_annot,
    	docker=docker,
    	outputfile=outputfile
	}

}
