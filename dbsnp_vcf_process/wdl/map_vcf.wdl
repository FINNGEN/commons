workflow map_vcf{
    call vcf_task{
    }
}

task vcf_task{
    String vcf_dl_link
    File chrom_map
    String docker

    String vcf_index_dl_link = vcf_dl_link+".tbi"
    Array[Array[String]] chrom_arr = read_tsv(chrom_map)
    Array[Array[String]] transposed = transpose(chrom_arr)
    Array[String] old_chr_names = transposed[0]
    

    command <<<
        set -euxo pipefail
        #download
        curl -o vcf_file.vcf.gz ${vcf_dl_link}
        curl -o vcf_file.vcf.gz.tbi ${vcf_index_dl_link}
        #rename & resctrict chromosomes
        bcftools view -r ${sep=',' old_chr_names} vcf_file.vcf.gz|bcftools annotate --rename-chr ${chrom_map}|bgzip -@4 > output.vcf.gz 
        tabix -p vcf output.vcf.gz
    >>>

    runtime {
        docker: docker
        cpu: "4"
        memory: "8 GB"
        disks:"local-disk 100 HDD"
        preemptible: 0
    }

    output {
        File out_vcf = "output.vcf.gz"
        File out_tbi = "output.vcf.gz.tbi"    
    }
}