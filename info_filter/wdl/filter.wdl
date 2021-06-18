task filter {

    File variant_file
    File sumstat
    String base = basename(sumstat,".gz")

    command <<<

        python3 - <<EOF > ${base}

        import sys
        import gzip
        variant_file="${variant_file}"
        sumstat="${sumstat}"
        variants = {}
        with gzip.open(variant_file, 'rt') as f:
            for line in f:
                variants[line.strip().replace('chr', '').replace('X', '23').replace('Y', '24').replace('MT', '25').replace('M', '25')] = True
        print("${variant_file}",file=sys.stderr)
        print("${sumstat}",file=sys.stderr)
        print("${base}",file=sys.stderr)
        with gzip.open(sumstat, 'rt') as f:
            l = f.readline().strip()
            print(l)
            print(l,file=sys.stderr)
            for line in f:
                line = line.strip()
                s = line.split('\t')
                chr = s[0].replace('chr', '').replace('X', '23').replace('Y', '24').replace('MT', '25').replace('M', '25')
                id = chr + ':' + s[1] + ':' + s[2] + ':' + s[3]
                if id in variants:
                    print(line)
        EOF
        bgzip ${base}
        tabix -s1 -b2 -e2 ${base}.gz

    >>>

    output {
        File out = base+".gz"
        File out_tbi = base + ".gz.tbi"
    }

    runtime {

        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.7"
        cpu: 1
        # 40M variants in variant_file to look up takes about 4G
        memory: "6 GB"
        disks: "local-disk 100 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 0
        noAddress: true
    }
}

workflow filter_sumstats_by_variants {

    File sumstats_loc
    Array[String] sumstats = read_lines(sumstats_loc)

    scatter (sumstat in sumstats) {
        call filter {
            input: sumstat=sumstat
        }
    }
    output{
        Array[File] sumstat = filter.out
        Array[File] sumstat_tbi = filter.out_tbi
    }
}
