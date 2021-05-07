task snpstats {

    File bgenfile

    command {
        qctool -g ${bgenfile} -snp-stats -osnp ${basename(bgenfile)}.snp_stats.txt
    }

    output {
        File out = basename(bgenfile) + ".snp_stats.txt"
    }

    runtime {
        docker: "gcr.io/finngen-refinery-dev/bioinformatics:0.5"
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 50 HDD"
        zones: "europe-west1-b"
        preemptible: 1
    }
}

task combine {

    Array[File] snpstatfiles
    String outfile

    command <<<

        for file in ${sep=" " snpstatfiles}; do
            if [[ $file == *.gz ]]
            then
                gunzip -c $file | grep -v '^#' | sed 's/ /\t/g' > `basename $file`"DATAUNZIP"
            else
                grep -v '^#' $file | sed 's/ /\t/g' > `basename $file`"DATAUNZIP"
            fi
        done

        cat <(head -n 1 `basename ${snpstatfiles[0]}"DATAUNZIP"` | sed 's/alternate_ids/#alternate_ids/') \
        <(awk 'FNR>1' \
        `find *DATAUNZIP | sort -V | tr '\n' ' '`) | bgzip > ${outfile} && \
        tabix -S 1 -s 3 -b 4 -e 4 ${outfile}

    >>>

    output {
        File out = outfile
    }

    runtime {
        docker: "gcr.io/finngen-refinery-dev/bioinformatics:0.5"
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 20 HDD"
        zones: "europe-west1-b"
        preemptible: 1
    }
}

workflow qctool_snpstats {

    File bgenlistfile
    Array[String] bgenfiles = read_lines(bgenlistfile)

    scatter (bgenfile in bgenfiles) {
        call snpstats {
            input: bgenfile=bgenfile
        }
    }

    call combine {
        input: snpstatfiles=snpstats.out
    }
}
