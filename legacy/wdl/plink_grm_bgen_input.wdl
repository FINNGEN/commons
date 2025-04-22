task filter_prune {

    File bgenfile
    String docker

    File include_variants
    File include_samples
    Float geno_missing
    Float maf
    String ld
    String genome_build = "hg38"

    File samplefile = bgenfile + ".sample"
    String out = basename(sub(bgenfile, ".bgen$", "")) + ".filtered"

    Int mem = ceil(size(bgenfile, "G") * 1.5)
    Int local_disk = ceil(size(bgenfile, "G") * 4)

    command <<<

        sed 's/^chr//;s/^/chr/;s/^chr23/chrX/' ${include_variants} > include_variants.txt

        plink2 \
        --memory 100000 \
        --allow-extra-chr \
        --snps-only \
        --bgen ${bgenfile} ref-first \
        --sample ${samplefile} \
        --split-par ${genome_build} \
        --extract include_variants.txt \
        --keep ${include_samples} \
        --geno ${geno_missing} \
        --maf ${maf} \
        --indep-pairwise ${ld} \
        --make-bed \
        --out ${out} && \
        plink2 \
        --memory 100000 \
        --allow-extra-chr \
        --bfile ${out} \
        --extract ${out}.prune.in \
        --keep ${include_samples} \
        --make-bed \
        --out ${out}.grm
    >>>

    output {
        File log1 = out + ".log"
        File bed = out + ".grm.bed"
        File bim = out + ".grm.bim"
        File fam = out + ".grm.fam"
        File log2 = out + ".grm.log"
        File pin = out + ".prune.in"
        File pout = out + ".prune.out"
    }

    # 1.5 times the bgen size not enough memory for small chromosomes so give at minimum 32GB to be sure
    runtime {
        docker: "${docker}"
        memory: if mem < 32 then "32 GB" else mem + " GB" 
        disks: "local-disk ${local_disk} HDD"
        preemptible: 0
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        noAddress: true
    }
}

task merge_plink {

    Array[File] bedfiles
    Array[File] bimfiles
    Array[File] famfiles
    String docker

    String out

    command <<<

        for file in ${sep=" " bedfiles}; do
            echo $file | sed 's/\.bed$//' >> mergelist.txt
        done

        plink \
        --allow-extra-chr \
        --merge-list mergelist.txt \
        --make-bed \
        --out ${out}

    >>>

    output {
        File log = out + ".log"
        File bed = out + ".bed"
        File bim = out + ".bim"
        File fam = out + ".fam"
    }

    runtime {
        docker: "${docker}"
        memory: "8 GB"
        disks: "local-disk 200 HDD"
        preemptible: 0
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        noAddress: true
    }
}

workflow plink_grm_bgen_input {

    File bgenfilelist
    String docker
    Array[String] bgenfiles = read_lines(bgenfilelist)

    scatter (bgenfile in bgenfiles) {
        call filter_prune {
            input:
                bgenfile = bgenfile,
                docker = docker
        }
    }

    call merge_plink {
        input:
            bedfiles = filter_prune.bed,
            bimfiles = filter_prune.bim,
            famfiles = filter_prune.fam,
            docker = docker
    }
}
