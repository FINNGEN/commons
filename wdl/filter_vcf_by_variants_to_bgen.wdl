task filter_to_bgen {

    File variant_file
    File vcf
    File tbi = vcf + ".tbi"
    String base_vcf = basename(vcf)
    String genofield="GP"
    String ofiletype="bgen_v1.2"
    Int precision=8
    Float input_rounding_error=0.005

    command <<<

        ## depending on the number of variants and samples it may be faster to first tabix but i get a couple of duplicate variants for some reason
        mv ${vcf} ${base_vcf}
        #gunzip -c ${variant_file} | awk 'BEGIN{OFS="\t"}NR>1{split($1, a, ":"); print "chr"a[1],a[2],a[2]}' | sed 's/chr23/chrX/' | uniq > variants.tabix
        #tabix -hR variants.tabix ${vcf} | uniq | bgzip > ${base_vcf}

        python3 - ${variant_file} ${base_vcf} <<EOF
        import sys
        import gzip

        variants = {}
        with gzip.open(sys.argv[1], 'rt') as f:
            for line in f:
                variants[line.strip().replace('X', '23').replace('Y', '24').replace('MT', '25').replace('M', '25')] = True

        out_vcf = gzip.open(sys.argv[2] + '.filtered.vcf.gz', 'wt')

        with gzip.open(sys.argv[2], 'rt') as f:
            line = f.readline().strip()
            out_vcf.write(line + '\n')
            while line.startswith('##'):
                line = f.readline().strip()
                out_vcf.write(line + '\n')
            for line in f:
                line = line.strip()
                # TODO no need to split the whole row
                s = line.split('\t')
                chr = s[0].replace('chr', '').replace('X', '23').replace('Y', '24').replace('MT', '25').replace('M', '25')
                id = chr + ':' + s[1] + ':' + s[3] + ':' + s[4]
                if id in variants:
                    out_vcf.write(line + '\n')

        out_vcf.close()
        EOF

        qctool -g ${base_vcf}.filtered.vcf.gz -vcf-genotype-field ${genofield} -filetype vcf -og ${base_vcf}.filtered.bgen -bgen-compression zlib -ofiletype ${ofiletype} -bgen-bits ${precision} -bgen-permitted-input-rounding-error ${input_rounding_error}
    >>>

    output {
        File out_bgen = base_vcf + '.filtered.bgen'
    }

    runtime {

        docker: "gcr.io/finngen-refinery-dev/bioinformatics:0.5"
        cpu: 1
        # 40M variants in variant_file to look up takes about 4G
        memory: "2 GB"
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b"
        preemptible: 0
        noAddress: true
    }
}

task combine {

    Array[File] bgenfiles
    String outfile

    command <<<
        cat-bgen -g ${sep=" " bgenfiles} -og ${outfile} && \
        bgenix -g ${outfile} -index
    >>>

    runtime {
        docker: "gcr.io/finngen-refinery-dev/bioinformatics:0.5"
        cpu: 1
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b"
        preemptible: 0
    }

    output {
        File out = outfile
        File out_index = outfile + ".bgi"
    }
}

workflow filter_vcf_by_variants_to_bgen {

    File vcf_files
    Array[String] vcfs = read_lines(vcf_files)

    scatter (vcf in vcfs) {
        call filter_to_bgen {
            input: vcf=vcf
        }
    }

    call combine {
        input: bgenfiles=filter_to_bgen.out_bgen
    }
}
