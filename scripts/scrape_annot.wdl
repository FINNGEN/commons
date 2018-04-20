
task join_annot {

	Array[File] files
	File header
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
			cat <(head -n 1 ${header}) <(zcat ${sep=" " files}) | bgzip > annotated_variants.gz
			tabix -S 1 -b 2 -e 2 -s 1 annotated_variants.gz
		else
			n_join1=$(( `zcat ${files[0]} | head -n 1 | awk 'BEGIN{ FS="\t"} { print NF};'` + `head -n 1 "${external_annot}" | awk 'BEGIN{ FS="\t"} { print NF-1};'`))
			join -a 1 -t $'\t' -1 2 -2 1 <(cat <(head -n 1 ${header} | awk '{ print "#"$0}' ) <(zcat ${sep=" " files}) | nl -nln | sort -b -k 2,2 ) <( awk 'BEGIN{ FS="\t"} NR==1{ print "#"$0  } NR>1{ print $0}' "${external_annot}" | sort -b -k 1,1  ) | sort -g -k 2,2 | cut -f 1,3- | awk -v nf=$n_join1 'BEGIN{ FS="\t";OFS="\t"} { printf $0; if(NF<nf) { for(i=NF+1;i<=nf;i++) { printf "\tNA"  };  } printf "\n" }' | bgzip > annotated_variants.gz			
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

	runtime {
		docker: "${docker}"
		memory:"8G"
		disks: "local-disk ${local_disk} SSD"
		cpu:"1"
	}
	
	command <<<
		#set -o pipefail
    	zcat ${vcf} | cut -f 1-8 | gawk 'BEGIN{ FS="\t"; OFS="\t"; header=0} $1!~/^#/{ delete out; delete names; elems=split($8,a,";"); for(i=1;i<=elems;i++) { split(a[i],b,"="); out[b[1]]=b[2]; names[i]=b[1] }; asorti(out,outdest); s=out[outdest[1]];h=outdest[1]; for(ind=2;ind<=elems; ind++) { h=h"\t"outdest[ind];  s=s"\t"out[outdest[ind]] }; print $1":"$2":"$4":"$5,$1,$2,$4,$5,s   }' | bgzip > ${basename(vcf)}.annot.gz
    	zcat ${vcf} | grep -v "^#"  | head -n 1 > first
		gawk ' BEGIN{ FS="\t"; OFS="\t"; header=0} { delete out; delete names; elems=split($8,a,";"); for(i=1;i<=elems;i++) { split(a[i],b,"="); out[b[1]]=b[2]; names[i]=b[1] }; asorti(out,outdest); s=out[outdest[1]];h=outdest[1]; for(ind=2;ind<=elems; ind++) { h=h"\t"outdest[ind];  s=s"\t"out[outdest[ind]] }; if(NR==1) print "variant","chr","pos","ref","alt",h;} ' first > ${basename(vcf)}.annot.header
    	echo "Done converting"
	>>>

	output {
		File out="${basename(vcf)}.annot.gz"
		File header="${basename(vcf)}.annot.header"
	}

}

workflow scrape_annots {

	File vcfs 

	Array[String] files = read_lines(vcfs)

	String docker
	String memory= "8GB"
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
		header=extract.header[0],
		external_annot=external_annot,
    	docker=docker,
    	outputfile=outputfile
	}

}
