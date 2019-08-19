task filter_prune {

    File bedfile
    File bimfile = sub(bedfile, ".bed$", ".bim")
    File famfile = sub(bedfile, ".bed$", ".fam")
    String out = basename(sub(bedfile, ".bed$", "")) + ".filtered"
    File include_variants
    File include_samples
    Float geno_missing
    Float maf
    String ld

    command <<<
        plink2 \
        --memory 40000 \
        --allow-extra-chr \
        --snps-only \
        --bfile ${sub(bedfile, ".bed$", "")} \
        --extract ${include_variants} \
        --keep ${include_samples} \
        --geno ${geno_missing} \
        --maf ${maf} \
        --indep-pairwise ${ld} \
        --make-bed \
        --out ${out} && \
        plink2 \
        --memory 40000 \
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

    # R4 chr1 original plink file is 56G
    # R4 chr1 required 26G RAM in LD pruning
    runtime {
        docker: "gcr.io/finngen-refinery-dev/bioinformatics:0.3"
        memory: "42 GB"
        disks: "local-disk 200 HDD"
        preemptible: 0
    }
}

task merge_plink {

    Array[File] bedfiles
    Array[File] bimfiles
    Array[File] famfiles
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
        docker: "gcr.io/finngen-refinery-dev/bioinformatics:0.3"
        memory: "8 GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}

workflow plink_grm {

    File bedfilelist
    Array[String] bedfiles = read_lines(bedfilelist)

    scatter (bedfile in bedfiles) {
        call filter_prune {
            input: bedfile=bedfile
        }
    }

    call merge_plink {
        input: bedfiles=filter_prune.bed, bimfiles=filter_prune.bim, famfiles=filter_prune.fam
    }
}
