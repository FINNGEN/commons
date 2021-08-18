task filter {

    File vcf
    Float af
    Float af_ = 1-af
    String genofield
    String ofiletype
    Int precision
    Float input_rounding_error

    command <<<

        bcftools filter -i'AF<${af} || AF>${af_}' ${vcf} | \
        qctool -g - -filetype vcf -vcf-genotype-field ${genofield} -og ${basename(vcf)}.bgen -bgen-compression zlib -ofiletype ${ofiletype} -bgen-bits ${precision} -bgen-permitted-input-rounding-error ${input_rounding_error}

    >>>

    runtime {
        docker: "gcr.io/finngen-refinery-dev/bioinformatics:0.6"
        cpu: 2
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
    }

    output {
        File outbgen = "${basename(vcf)}.bgen"
    }
}

task combine {

    Array[File] bgenfiles
    String outfile

    command <<<
        cat-bgen -g ${sep=" " bgenfiles} -og ${outfile}
        bgenix -g ${outfile} -index
    >>>

    runtime {
        docker: "gcr.io/finngen-refinery-dev/bioinformatics:0.6"
        cpu: 1
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
    }

    output {
        File out = outfile
        File out_index = outfile + ".bgi"
    }
}

workflow filter_vcf_to_bgen {

    File files_to_conv
    Float af
    String genofield="GT"
    String ofiletype="bgen_v1.2"
    Int precision=8
    Float input_rounding_error=0.005

    Array[String] files = read_lines(files_to_conv)

    scatter( file in files) {

        call filter {
            input: vcf=file,
            af=af,
            genofield=genofield,
            ofiletype=ofiletype,
            precision=precision,
            input_rounding_error=input_rounding_error
        }
    }

    call combine {
        input: bgenfiles=filter.outbgen
    }
}
