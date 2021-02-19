task filter_af {

    File file
    String af_col
    Float min_af
    String outfile = basename(file)
    String docker

    command <<<
        gunzip -c ${file} | awk '
        NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print $0}
        NR>1 && $a["${af_col}"]>${min_af}&&(1-$a["${af_col}"])>${min_af}' | bgzip > ${outfile}
    >>>

    output {
        File out = outfile
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }
}

workflow filter_gwas {

    File sumstats_loc
    Array[String] sumstat_files = read_lines(sumstats_loc)

    scatter (file in sumstat_files) {
        call filter_af {
            input: file=file
        }
    }
}
