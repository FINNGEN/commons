workflow fgfactory_pass_samples {
	String docker
	Array[File] VCF_ARRAY
	File CHIP_INFO
	String TO_EXCLUDE 
	String RELEASE
	
	scatter(VCF in VCF_ARRAY) {
		call get_values {input: docker=docker, VCF=VCF, CHIP_INFO=CHIP_INFO, TO_EXCLUDE=TO_EXCLUDE}
	}
	
	call concatenate_values {input: docker=docker, values_array=get_values.out_values}
	
	scatter(VCF in VCF_ARRAY) {
		call calculate_vcf_sizes {input: docker=docker, VCF=VCF}
	}
	
	call sum_all_file_sizes {input: docker=docker, size_array = calculate_vcf_sizes.out_size_array}
	
	call samplewise_summary {input: docker=docker, SUMMARY_DATA = concatenate_values.out_all_values, VCF_ARRAY=VCF_ARRAY, TOTAL_VCF_FILE_SIZE = sum_all_file_sizes.out_total_vcf_file_size, RELEASE=RELEASE}
	
	output {
		
		Array[File] out_values = get_values.out_values
		File out_all_values = concatenate_values.out_all_values
		File out_fgfactory_pass_samples = samplewise_summary.out_fgfactory_pass_samples
		Array[File] outLists = samplewise_summary.outLists
		
	}
}

task get_values {
	String docker
	String VCF
	File CHIP_INFO
	String TO_EXCLUDE
	String DATA_NAME = sub(basename(VCF), "${TO_EXCLUDE}", "")
	
	command <<<
	
	file_name=`basename ${VCF}`
	chip=`cat ${CHIP_INFO} | grep ${DATA_NAME} | cut -f1`
	variants=`cat ${CHIP_INFO} | grep ${DATA_NAME} | cut -f3`
	release=`cat ${CHIP_INFO} | grep ${DATA_NAME} | cut -f4`
	cohort=`cat ${CHIP_INFO} | grep ${DATA_NAME} | cut -f5`
	echo $file_name$'\t'$chip$'\t'$variants$'\t'$release$'\t'$cohort$'\t'${DATA_NAME} > ${DATA_NAME}_values.txt
	
	>>>
	
	output {
		
		File out_values = "${DATA_NAME}_values.txt"
		
	}
	
	runtime {
		docker: docker
		memory: "2 GB"
		cpu: 1
		disks: "local-disk 200 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
	}
}

task concatenate_values {
	
	String docker
	Array[File] values_array
	
	command <<<
	
	cat ${sep=' ' values_array} >> summary_data.txt
	
	>>>
	
	output {
		File out_all_values = "summary_data.txt"
	}
	
	runtime {
		docker: docker
		memory: "2 GB"
		cpu: 1
		disks: "local-disk 200 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
	}
}

task calculate_vcf_sizes {
	File VCF
	String docker
	
	command { }
	
	output {
		
		Float out_size_array = size("${VCF}", "GB")
		
	}
	
	runtime {
		docker: docker
		memory: "2 GB"
		cpu: 1          
		disks: "local-disk 200 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
	}
}

task sum_all_file_sizes {
	String docker
	Array[Float] size_array 
	
	command <<<
	
	python -c "print(${sep="+" size_array})"
	
	>>>
	
	output {
		
		Float out_total_vcf_file_size = read_float(stdout())
		
	}
	
	runtime {
		docker: docker
		memory: "2 GB"
		cpu: 1
		disks: "local-disk 200 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
	}              
}

task samplewise_summary {
	String docker
	Array[File] VCF_ARRAY
	Float TOTAL_VCF_FILE_SIZE
	File SUMMARY_DATA
	String RELEASE
	Int allocated_disk_size = ceil(TOTAL_VCF_FILE_SIZE*2)
	String dollar = "$"
	
	command <<<
	
		set -e
		
		for file in ${sep=" " VCF_ARRAY}; do mv $file .; done
		
		# Create sample-wise summary data file 
		
		# Run as:bash Samples_Chips_Variants_array_v2.0.sh input_data/R2_dataset_summary_table.txt ouput_data
		#Table of summary data with columns like
		#/fs/ffworkspace/fgfact_prod_R2_0718/output/AxiomGT1_b01_26072018/AxiomGT1_b01.calls_tags_edited_drem_drem_chr22.vcf.gz  Axiom_FinnGen1.r1       507250  R1       COHORT  AxiomGT1_b01
		#/fs/ffworkspace/fgfact_prod_R2_0718/output/AxiomGT1_b02_26072018/AxiomGT1_b02.calls_tags_edited_drem_drem_chr22.vcf.gz  Axiom_FinnGen1.r1       510673  R2       COHORT2 AxiomGT1_b02

		#extract bash arrays from summary datatable
		input_chr=($(awk 'BEGIN{FS=OFS="\t"} {print $1}' ${SUMMARY_DATA}))
		chip_type=($(awk 'BEGIN{FS=OFS="\t"} {print $2}' ${SUMMARY_DATA}))
		variants_passQC=($(awk 'BEGIN{FS=OFS="\t"} {print $3}' ${SUMMARY_DATA}))
		Cohort_dataset=($(awk 'BEGIN{FS=OFS="\t"} {print $6}' ${SUMMARY_DATA}))
		release_dataset=($(awk 'BEGIN{FS=OFS="\t"} {print $4}' ${SUMMARY_DATA}))
		original_cohort=($(awk 'BEGIN{FS=OFS="\t"} {print $5}' ${SUMMARY_DATA}))


		#Run a loop to join Finngen IDs and agregate data from the summary table
		for ((i = 0; i < "${dollar}{#input_chr[@]}" && i < "${dollar}{#chip_type[@]}" && i < "${dollar}{#variants_passQC[@]}" && i < "${dollar}{#Cohort_dataset[@]}"; i++))
		do
		zcat "${dollar}{input_chr[i]}" | grep -m 1 '#CHROM' |  tr '\t' '\n' \
		| eval "sed s/^/${dollar}{chip_type[i]}':'/" \
		| eval "sed s/^/${dollar}{variants_passQC[i]}':'/" \
		| eval "sed s/^/${dollar}{release_dataset[i]}':'/" \
		| eval "sed s/^/${dollar}{original_cohort[i]}':'/" \
		| eval "sed s/^/${dollar}{Cohort_dataset[i]}':'/" | tail -n +10 >> samplewise_summary_data.txt
		done
		
		sort -t":" -k1,1 samplewise_summary_data.txt > fgfactory_pass_samples_${RELEASE}_no_header.txt
		
		# Add a header line 
		
		echo "BATCH:COHORT:RELEASE:VARIANTS:CHIP:FINNGENID" > header.txt
		cat header.txt fgfactory_pass_samples_${RELEASE}_no_header.txt > fgfactory_pass_samples_${RELEASE}.txt 
		
		
		# Also prepare batch-wise inclusion list for use in chipd pipeline and next release production 
		# The samplelists should end with string .samplelist 
		
		batches=`cut -f1 -d":" fgfactory_pass_samples_${RELEASE}_no_header.txt | sort | uniq`
		for i in $batches; do grep $i fgfactory_pass_samples_${RELEASE}.txt | cut -f6 -d":" > pass_samples_${RELEASE}_$i.samplelist; done
		
	>>>
		
	output {
		File out_fgfactory_pass_samples = "fgfactory_pass_samples_${RELEASE}.txt"
		Array[File] outLists = glob("pass_samples_*")
	}
	
	runtime {
		docker: docker
		memory: "2 GB"
		cpu: 1
		disks: "local-disk ${allocated_disk_size} HDD" 
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
	}
}
