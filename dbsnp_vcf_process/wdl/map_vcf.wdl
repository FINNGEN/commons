workflow map_vcf{
    call vcf_task{
    }
}

task vcf_task{
    String vcf_dl_link
    File chrom_map
    File fasta
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
        cat ${chrom_map}|sed "s/23$/X/g;s/24$/Y/g;s/25$/M/g;s/\t/\tchr/g" > fasta_compliant_chrom_map
        paste <(cut -f 2 fasta_compliant_chrom_map) <(cut -f 2 ${chrom_map}) > chrtonum
        #rename & resctrict chromosomes
        bcftools view   vcf_file.vcf.gz -m 2 -M 2 -r ${sep=',' old_chr_names} -Ou| \
        bcftools annotate --rename-chr fasta_compliant_chrom_map -Ou | \
        bcftools norm -f ${fasta} -c wx  -Ou| \
        bcftools annotate --rename-chr chrtonum -Ov|bgzip -@4 > output.vcf.gz 
        tabix -p vcf output.vcf.gz
    >>>

    runtime {
        docker: docker
        cpu: "4"
        memory: "8 GB"
        disks:"local-disk 100 HDD"
        preemptible: 0
        zones: "europe-west1-b europe-west1-c europe-west1-d"
    }

    output {
        File out_vcf = "output.vcf.gz"
        File out_tbi = "output.vcf.gz.tbi"    
    }
}