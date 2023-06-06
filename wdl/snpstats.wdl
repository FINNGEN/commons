version 1.0

task snpstats {

    input {
        File genofile
    }

    command {
        qctool -g ${genofile} -snp-stats -osnp ${basename(genofile)}.snp_stats.txt
    }

    output {
        File out = basename(genofile) + ".snp_stats.txt"
    }

    runtime {
        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8"
        cpu: 1
        memory: "3 GB"
        disks: "local-disk " + 2 * (ceil(size(genofile, "G"))) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 0
    }
}

task combine {

    input {
        Array[File] snpstatfiles
        String outfile
    }

    command <<<

        for file in ~{sep=" " snpstatfiles}; do
            if [[ $file == *.gz ]]
            then
                gunzip -c $file | grep -v '^#' | sed 's/ /\t/g' > `basename $file`"DATAUNZIP"
            else
                grep -v '^#' $file | sed 's/ /\t/g' > `basename $file`"DATAUNZIP"
            fi
        done

        cat <(head -n 1 `basename ~{snpstatfiles[0]}"DATAUNZIP"` | sed 's/alternate_ids/#alternate_ids/') \
        <(awk 'FNR>1' \
        `find *DATAUNZIP | sort -V | tr '\n' ' '`) | bgzip > ~{outfile} && \
        tabix -S 1 -s 3 -b 4 -e 4 ~{outfile}

    >>>

    output {
        File out = outfile
        File tbi = outfile + ".tbi"
    }

    runtime {
        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8"
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 20 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 1
    }
}

workflow qctool_snpstats {

    input {
        File genolistfile
        Array[String] genofiles = read_lines(genolistfile)
    }

    scatter (genofile in genofiles) {
        call snpstats {
            input: genofile=genofile
        }
    }

    call combine {
        input: snpstatfiles=snpstats.out
    }
}
