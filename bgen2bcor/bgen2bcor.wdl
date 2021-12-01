version 1.0

workflow bgen2bcor {

    input {
        File bgen
    }

    call split2chr {
        input:
            bgen = bgen
    }

    scatter (chr in range(23)) {
        call convert2bcor {
            input:
                bgen = split2chr.chr_bgens[chr],
                bgi = split2chr.chr_bgis[chr]
        }
    }

    output {
        Array[File] chr_bgens = split2chr.chr_bgens
        Array[File] chr_bgis = split2chr.chr_bgis
        Array[File] chr_bcors = convert2bcor.bcor
    }
}


task split2chr {

    input {
        File bgen

        String docker
        Int cpu
        Int mem

        File bgi = bgen + ".bgi"
        String base = basename(bgen, ".bgen")
    }

    command <<<

        set -euxo pipefail

        for chr in {1..23}
        do
            if [ ${chr} == 23 ]; then chr=X; fi
            bgenix -g ~{bgen} -incl-range chr${chr}:0- > ~{base}_chr${chr}.bgen
            bgenix -index -g ~{base}_chr${chr}.bgen
        done

    >>>

    output {
        Array[File] chr_bgens = glob("*_chr*.bgen")
        Array[File] chr_bgis = glob("*_chr*.bgen.bgi")
    }

    runtime {
        docker: "~{docker}"
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disks: "local-disk " + 3*ceil(size(bgen, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 0
        noAddress: true
    }
}


task convert2bcor {

    input {
        File bgen
        File bgi

        String docker
        Int cpu
        Int mem
        Int variant_window_size
        Float ld_thold
        String accuracy

        String base = basename(bgen, ".bgen")
    }

    command <<<

        set -euxo pipefail

        n_threads=`grep -c ^processor /proc/cpuinfo`
        
        # ldstore v1.1
        ldstore --bgen ~{bgen} --variant-window-size ~{variant_window_size} --ld-thold ~{ld_thold} --accuracy ~{accuracy} --n-threads ${n_threads} --bcor ~{base}.bcor
        ldstore --bcor ~{base}.bcor --merge ${n_threads}

    >>>

    output {
        File bcor = base + ".bcor"
    }

    runtime {
        docker: "~{docker}"
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disks: "local-disk " + 50 * ceil(size(bgen, "G") + size(bgi, "G")) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 0
        noAddress: true
    }
}
