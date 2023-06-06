task filter_to_bgen {

    File variant_file
    File vcf
    String base_vcf = basename(vcf)
    String genofield
    String ofiletype="bgen_v1.2"
    Int precision=8
    Float input_rounding_error=0.005

    command <<<

        awk 'NR==FNR{a[$1]} NR>FNR && ($0~"^#" || $3 in a)' ${variant_file} <(zcat -f ${vcf}) | \
        sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' | \
	qctool -g - -vcf-genotype-field ${genofield} -filetype vcf -og ${base_vcf}.filtered.bgen -os ${base_vcf}.filtered.bgen.sample \
	-bgen-compression zlib -ofiletype ${ofiletype} -bgen-bits ${precision} -bgen-permitted-input-rounding-error ${input_rounding_error}
    >>>

    output {
        File out_bgen = base_vcf + ".filtered.bgen"
        File out_sample = base_vcf + ".filtered.bgen.sample"
    }

    runtime {
        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8"
        cpu: 3
        memory: "2 GB"
        disks: "local-disk 300 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
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
        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8"
        cpu: 1
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 0
	noaddress: true
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
