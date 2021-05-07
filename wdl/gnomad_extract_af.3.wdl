task extract_af {

    File vcf
    String base_vcf = basename(vcf)

    command <<<

        mv ${vcf} ${base_vcf}

        python3 - ${base_vcf} <<EOF

        import sys
        import gzip

        #POPS = ['all']
        POPS = ['all', 'afr', 'ami', 'amr', 'asj', 'eas', 'fin', 'nfe', 'sas', 'oth']

        out_tsv = [gzip.open('ref_${base_vcf}_' + pop + '.gz', 'wt') for pop in POPS]
        for i,out in enumerate(out_tsv):
            if POPS[i] == 'all':
                out.write('#chr\tpos\tref\talt\taf_alt\tfilter\tan\trsid\n')
            else:
                out.write('#chr\tpos\tref\talt\taf_alt\tfilter\tan\n')

        af_fields = [('AF_' + pop).replace('_all', '') for pop in POPS]
        with gzip.open(sys.argv[1], 'rt') as f:
            line = f.readline().strip()
            while line.startswith('##'):
                line = f.readline().strip()
            for line in f:
                line = line.strip()
                s = line.split('\t')
                chrom = s[0].replace('chr', '').replace('X', '23')
                cpra = chrom + '\t' + s[1] + '\t' + s[3] + '\t' + s[4]
                info = {dat.split('=')[0]: dat.split('=')[1] for dat in s[7].split(';') if '=' in dat}
                for i,out in enumerate(out_tsv):
                    af = float(info[af_fields[i]]) if af_fields[i] in info else 'NA'
                    if af_fields[i] == 'AF': # for overall AF also write rsid
                        out.write(cpra + '\t' + str(af) + '\t' + s[6] + '\t' + info['AN'] + '\t' + s[2].replace('.', 'NA') + '\n')
                    else:
                        out.write(cpra + '\t' + str(af) + '\t' + s[6] + '\t' + info['AN'] + '\n')

        for out in out_tsv:
            out.close()

        EOF

    >>>

    output {
        Array[File] out = glob('ref_*')
    }

    runtime {
        docker: "gcr.io/finngen-refinery-dev/bioinformatics:0.6"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 200 SSD"
        zones: "europe-west1-b"
        preemptible: 0
        noAddress: true
    }
}

task combine {

    Array[Array[File]] files2D
    Array[File] files = flatten(files2D)
    #Array[String] pops = ["all"]
    Array[String] pops = ["all", "afr", "ami", "amr", "asj", "eas", "fin", "nfe", "sas", "oth"]
    String dollar = "$"

    command <<<

        for file in ${sep=" " files}; do
            mv $file `basename $file`
        done
        while read pop; do
            files=(`ls *_$pop* | sort -V`)
            cat <(gunzip -c ${dollar}{files[0]} | head -1) \
                <(for file in ${dollar}{files[@]}; do gunzip -c $file | tail -n+2; done) \
            | bgzip > gnomad_v3_b38_ref_$pop.gz
            tabix -s 1 -b 2 -e 2 gnomad_v3_b38_ref_$pop.gz
        done < ${write_lines(pops)}

    >>>

    output {
        Array[File] out = glob("*.gz")
        Array[File] tbi = glob("*.tbi")
    }

    runtime {
        docker: "gcr.io/finngen-refinery-dev/bioinformatics:0.6"
        cpu: 1
        memory: "3 GB"
        disks: "local-disk 200 SSD"
        zones: "europe-west1-b"
        preemptible: 0
        noAddress: true
    }
}

workflow gnomad_extract_af {

    File vcf_files
    Array[String] vcfs = read_lines(vcf_files)

    scatter (vcf in vcfs) {
        call extract_af {
            input: vcf=vcf
        }
    }

    call combine {
        input: files2D=extract_af.out
    }
}
