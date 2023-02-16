version 1.0

workflow fgfactory_pass_samples {

  input {
    File CHIP_INFO
    Array[File] VCF_ARRAY
    String docker
    String TO_EXCLUDE
    String RELEASE
  }

  scatter (VCF in VCF_ARRAY) {
    call get_values {
      input:
        docker = docker,
        VCF = VCF,
        CHIP_INFO = CHIP_INFO,
        TO_EXCLUDE = TO_EXCLUDE
    }
  }

  call concatenate_values {
    input:
      docker = docker,
      values_array = get_values.out_values
  }

  call samplewise_summary {
    input:
      docker = docker,
      VCF_ARRAY = VCF_ARRAY,
      SUMMARY_DATA = concatenate_values.out_all_values,
      RELEASE = RELEASE
  }

  output {
    Array[File] out_values = get_values.out_values
    File out_fgfactory_pass_samples = samplewise_summary.out_fgfactory_pass_samples
    Array[File] outLists = samplewise_summary.outLists
    File out_all_values = concatenate_values.out_all_values
  }

}

task get_values {

  input {
    String docker
    String VCF
    File CHIP_INFO
    String TO_EXCLUDE
    String DATA_NAME = sub(basename(VCF), "~{TO_EXCLUDE}", "")
  }

  command <<<
    
    file_name=`basename ~{VCF}`
    dataset_name=$(echo ~{DATA_NAME} | sed 's/_relift$//;s/_nosymmetric$//') # Remove extension to guarantee matching
    chip=`cat ~{CHIP_INFO} | grep $dataset_name | cut -f1`
    variants=`cat ~{CHIP_INFO} | grep $dataset_name | cut -f3`
    release=`cat ~{CHIP_INFO} | grep $dataset_name | cut -f4`
    cohort=`cat ~{CHIP_INFO} | grep $dataset_name | cut -f5`
    echo $file_name$'\t'$chip$'\t'$variants$'\t'$release$'\t'$cohort$'\t'~{DATA_NAME} > ~{DATA_NAME}_values.txt
    
  >>>

  output {
    File out_values = "${DATA_NAME}_values.txt"
  }

  runtime {
    preemptible: 2
    disks: "local-disk 200 HDD"
    docker: docker
    cpu: 1
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    memory: "2 GB"
  }

}

task concatenate_values {

  input {
    String docker
    Array[File] values_array
  }

  command <<<

    cat ~{sep=" "  values_array} >> summary_data.txt
    
  >>>

  output {
    File out_all_values = "summary_data.txt"
  }

  runtime {
    preemptible: 2
    disks: "local-disk 200 HDD"
    docker: docker
    cpu: 1
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    memory: "2 GB"
  }

}

task samplewise_summary {

  input {
    String docker
    Array[File] VCF_ARRAY
    File SUMMARY_DATA
    String RELEASE
    Int allocated_disk_size = ceil(size(VCF_ARRAY, "GB") * 2)
  }

  command <<<

    set -e
    
    for file in ~{sep=" " VCF_ARRAY}; do mv $file .; done
    
    # Create sample-wise summary data file
    
    # Run as:bash Samples_Chips_Variants_array_v2.0.sh input_data/R2_dataset_summary_table.txt ouput_data
    #Table of summary data with columns like
    #/fs/ffworkspace/fgfact_prod_R2_0718/output/AxiomGT1_b01_26072018/AxiomGT1_b01.calls_tags_edited_drem_drem_chr22.vcf.gz  Axiom_FinnGen1.r1       507250  R1       COHORT  AxiomGT1_b01
    #/fs/ffworkspace/fgfact_prod_R2_0718/output/AxiomGT1_b02_26072018/AxiomGT1_b02.calls_tags_edited_drem_drem_chr22.vcf.gz  Axiom_FinnGen1.r1       510673  R2       COHORT2 AxiomGT1_b02

    #extract bash arrays from summary datatable
    input_chr=($(awk 'BEGIN{FS=OFS="\t"} {print $1}' ~{SUMMARY_DATA}))
    chip_type=($(awk 'BEGIN{FS=OFS="\t"} {print $2}' ~{SUMMARY_DATA}))
    variants_passQC=($(awk 'BEGIN{FS=OFS="\t"} {print $3}' ~{SUMMARY_DATA}))
    Cohort_dataset=($(awk 'BEGIN{FS=OFS="\t"} {print $6}' ~{SUMMARY_DATA}))
    release_dataset=($(awk 'BEGIN{FS=OFS="\t"} {print $4}' ~{SUMMARY_DATA}))
    original_cohort=($(awk 'BEGIN{FS=OFS="\t"} {print $5}' ~{SUMMARY_DATA}))


    #Run a loop to join Finngen IDs and agregate data from the summary table
    for ((i = 0; i < "${#input_chr[@]}" && i < "${#chip_type[@]}" && i < "${#variants_passQC[@]}" && i < "${#Cohort_dataset[@]}"; i++))
    do
    zcat "${input_chr[i]}" | grep -m 1 '#CHROM' |  tr '\t' '\n' \
      | eval "sed s/^/${chip_type[i]}':'/" \
      | eval "sed s/^/${variants_passQC[i]}':'/" \
      | eval "sed s/^/${release_dataset[i]}':'/" \
      | eval "sed s/^/${original_cohort[i]}':'/" \
      | eval "sed s/^/${Cohort_dataset[i]}':'/" | tail -n +10 >> samplewise_summary_data.txt
    done
    
    sort -t":" -k1,1 samplewise_summary_data.txt > fgfactory_pass_samples_~{RELEASE}_no_header.txt
    
    # Add a header line 
    
    echo "BATCH:COHORT:RELEASE:VARIANTS:CHIP:FINNGENID" > header.txt
    cat header.txt fgfactory_pass_samples_~{RELEASE}_no_header.txt > fgfactory_pass_samples_~{RELEASE}.txt 
    
    
    # Also prepare batch-wise inclusion list for use in chipd pipeline and next release production 
    # The samplelists should end with string .samplelist 
    
    batches=`cut -f1 -d":" fgfactory_pass_samples_~{RELEASE}_no_header.txt | sort | uniq`
    for i in $batches; do grep $i fgfactory_pass_samples_~{RELEASE}.txt | cut -f6 -d":" > pass_samples_~{RELEASE}_$i.samplelist; done
    
  >>>

  output {
    File out_fgfactory_pass_samples = "fgfactory_pass_samples_~{RELEASE}.txt"
    Array[File] outLists = glob("pass_samples_*")
  }

  runtime {
    preemptible: 2
    disks: "local-disk ~{allocated_disk_size} HDD"
    docker: docker
    cpu: 1
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    memory: "2 GB"
  }

}
