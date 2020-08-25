task split {

    File bgen_file
    String base = basename(bgen_file, ".bgen")
    File bgi = bgen_file + ".bgi"
    Int n_variants_per_chunk_approx
    String docker

    command <<<
        bgenix -g ${bgen_file} -list | tail -n+3 | cut -f2 > ${base}.variants
        n_var=`wc -l ${base}.variants | cut -d' ' -f1`
        n_chunks=$((n_var / ${n_variants_per_chunk_approx}))
        n_variants_per_chunk_actual=$((n_var / n_chunks + 1))
        n_cpu=`nproc --all`
        echo ${bgen_file}: $n_var variants, $n_chunks chunks, $n_variants_per_chunk_actual variants per chunk, $n_cpu cpus
        split -d -a 3 -l $n_variants_per_chunk_actual ${base}.variants ${base}_
        echo `date` start split
        ls ${base}_* | xargs -n 1 -P $n_cpu bash -c 'bgenix -g ${bgen_file} -incl-rsids $0 > $0.bgen'
        echo `date` end split
        ls -1 *.bgen
    >>>

    output {
        Array[File] chunks = glob("*.bgen")
    }

    runtime {
        docker: "${docker}"
        cpu: 4 # with 2T HDD, more than 4 CPUs doesn't help. with 2T SSD, can increase to 8 CPUs but HDD is fast enough for e.g. full UKB
        disks: "local-disk 2000 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }
}

workflow split_bgen {

    File bgen_files_loc
    Array[String] bgen_files = read_lines(bgen_files_loc)

    scatter (bgen_file in bgen_files) {
        call split {
            input: bgen_file=bgen_file
        }
    }
}
