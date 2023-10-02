task concat {

  Int chr
  File header_file
  Array[File] files

  command <<<
    zcat ${header_file} ${sep=" " files} | \
    bgzip > daly_finnish_callset_gnomad_v4_chr${chr}.vcf.gz
    tabix -p vcf daly_finnish_callset_gnomad_v4_chr${chr}.vcf.gz
  >>>

  runtime {
    docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8"
    cpu: 2
    disks: "local-disk 600 HDD"
    memory: "2 GB"
    zones: "us-east1-b us-east1-c us-east1-d"
    preemptible: 1
  }

  output {
    File vcf = "daly_finnish_callset_gnomad_v4_chr${chr}.vcf.gz"
    File tbi = "daly_finnish_callset_gnomad_v4_chr${chr}.vcf.gz.tbi"
  }
}

workflow concat_vcf {

  String header_file
  File chr_part_files_loc
  Array[Array[String]] chr_part_files = read_tsv(chr_part_files_loc)

  scatter (i in range(length(chr_part_files))) {
    call concat {
      input: chr=i+1, header_file=header_file, files=chr_part_files[i]
    }
  }
}
