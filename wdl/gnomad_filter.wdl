task filter {

    File variant_file
    File vcf
    String base_vcf = basename(vcf)
    String out_suffix = ".filtered.tsv.gz"
    Int build

    command <<<

        mv ${vcf} ${base_vcf}

        python3 - ${variant_file} ${base_vcf} <<EOF
        import sys
        import gzip

        FIELDS = []
        if ${build} == 37:
            #gnomad 2.1.1 (build 37)
            FIELDS = sorted(['AF', 'AF_nfe_seu', 'AF_afr', 'AF_nfe_onf', 'AF_amr', 'AF_eas', 'AF_nfe_nwe', 'AF_nfe_est', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth'])
        elif ${build} == 38:
            #gnomad 3 (build 38) - doesn't have European subpops yet in gnomad 3.0
            FIELDS = sorted(['AF', 'AF_ami', 'AF_oth', 'AF_afr', 'AF_sas', 'AF_raw', 'AF_asj', 'AF_fin', 'AF_amr', 'AF_nfe', 'AF_eas', 'AF_popmax'])

        variants = {}
        with gzip.open(sys.argv[1], 'rt') as f:
            for line in f:
                variants[line.strip().replace('X', '23').replace('Y', '24').replace('MT', '25').replace('M', '25')] = True

        out_tsv = gzip.open(sys.argv[2] + '${out_suffix}', 'wt')

        out_tsv.write('#CHR\tPOS\tREF\tALT\tRSID\t' + '\t'.join(FIELDS) + '\tconsequence\tgene\tensg\n')

        with gzip.open(sys.argv[2], 'rt') as f:
            line = f.readline().strip()
            while line.startswith('##'):
                line = f.readline().strip()
            for line in f:
                line = line.strip()
                s = line.split('\t')
                chr = s[0].replace('chr', '').replace('X', '23').replace('Y', '24').replace('MT', '25').replace('M', '25')
                id = chr + ':' + s[1] + ':' + s[3] + ':' + s[4]
                if id in variants:
                    info = {dat.split('=')[0]: dat.split('=')[1] for dat in s[7].split(';') if '=' in dat}
                    vep = info['vep'].split('|')
                    consequence = vep[1].split('&')[0]
                    if consequence == '' or consequence == '.':
                        consequence = 'NA'
                    gene = vep[3] if vep[3] != '' else 'NA'
                    ensg = vep[4] if vep[3] != '' else 'NA'
                    if s[2] == '.':
                        s[2] = 'NA'
                    out_tsv.write('\t'.join([chr, s[1], s[3], s[4], s[2]]) + '\t' + '\t'.join(['{:.2e}'.format(float(info[key])) if key in info else 'NA' for key in FIELDS]) + '\t' + '\t'.join([consequence, gene, ensg]) + '\n')

        out_tsv.close()
        EOF

    >>>

    output {
        File out_tsv = base_vcf + out_suffix
    }

    runtime {

        docker: "gcr.io/finngen-refinery-dev/bioinformatics:0.5"
        cpu: 1
        # 40M variants in variant_file to look up takes about 4G
        memory: "2 GB"
        disks: "local-disk 200 SSD"
        zones: "europe-west1-b"
        preemptible: 1
        noAddress: true
    }
}

task combine {

    Array[File] files
    String out_name

    command <<<

        cat <(gunzip -c ${files[0]} | head -1) \
        <(for file in ${sep=" " files}; do gunzip -c $file | tail -n+2; done) \
        | bgzip > ${out_name}

        tabix -s 1 -b 2 -e 2 ${out_name}
    >>>

    output {
        File anno = out_name
        File anno_tbi = out_name + ".tbi"
    }

    runtime {

        docker: "gcr.io/finngen-refinery-dev/bioinformatics:0.5"
        cpu: 1
        memory: "3 GB"
        disks: "local-disk 100 SSD"
        zones: "europe-west1-b"
        preemptible: 1
        noAddress: true
    }
}

workflow gnomad_filter {

    File vcf_files
    Array[String] vcfs = read_lines(vcf_files)

    scatter (vcf in vcfs) {
        call filter {
            input: vcf=vcf
        }
    }

    call combine {
        input: files=filter.out_tsv
    }
}
