task filter_merge {

    File vcf1
    File vcf2
    String base1 = basename(vcf1, ".vcf.gz")
    String base2 = basename(vcf2, ".vcf.gz")
    File sample_pairs
    String outname = basename(vcf1, ".vcf.gz") + "_" + basename(vcf2, ".vcf.gz")

    command <<<

        set -euxo pipefail

        cut -f1 ${sample_pairs} | awk '{print 0,$0}' > keep_vcf1_samples
        cut -f2 ${sample_pairs} | awk '{print 0,$0}' > keep_vcf2_samples

        # convert vcfs to bed keeping wanted samples
        plink2 --allow-extra-chr --vcf ${vcf1} --keep keep_vcf1_samples --freq --make-bed --out ${base1}
        plink2 --allow-extra-chr --vcf ${vcf2} --keep keep_vcf2_samples --freq --make-bed --out ${base2}

        # filter beds so they will have the same variants
        plink2 --allow-extra-chr --extract <(cut -f2 ${base2}.bim) --bfile ${base1} --make-bed --out ${base1}_shared
        plink2 --allow-extra-chr --extract <(cut -f2 ${base1}.bim) --bfile ${base2} --make-bed --out ${base2}_shared

        # underscores cause trouble with vcf conversions
        sed -i 's/_/-/g' ${base1}_shared.fam
        sed -i 's/_/-/g' ${base2}_shared.fam

        # in case sample names are shared between datasets, add suffix to one dataset's samples
        awk 'BEGIN{OFS="\t"} {$2=$2"-DATA1"} 1' ${base1}_shared.fam > temp && mv temp ${base1}_shared.fam
        awk 'BEGIN{OFS="\t"} {$1=$1"-DATA1"} 1' ${sample_pairs} > sample_pairs_suffixed

        # merge datasets
        plink --keep-allele-order --allow-extra-chr --bfile ${base1}_shared --bmerge ${base2}_shared --make-bed --out merged

        # convert to vcf
        plink2 --bfile merged --recode vcf-iid bgz --output-chr chrM --out ${outname}
        #tabix -p vcf ${outname}.vcf.gz

    >>>

    runtime {
        cpu: 4
        memory: "3G"
        disks: "local-disk 200 HDD"
        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.7"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 1
    }

    output {
        File out_vcf = outname + ".vcf.gz"
        File out_sample_pairs = "sample_pairs_suffixed"
    }
}

task concordance {

    File vcf
    File sample_pairs

    command <<<

        set -euxo pipefail

        # underscores were converted to dashes in the filter_merge task, do the same here to match samples right
        sed -i 's/_/-/g' ${sample_pairs}

        python3 <<EOF | bgzip > concordance.gz

        """
        computes variant-wise concordance stats
        inputs: a vcf file and a text file with a list of corresponding samples one pair of samples on a row
        prints stats to stdout and messsages to stderr
        """

        import gzip
        import sys

        def eprint(*args, **kwargs):
            print(*args, file=sys.stderr, **kwargs)

        def read_sample_index(vcf):
            with gzip.open(vcf, 'rt') as f:
                line = f.readline().strip()
                while line.startswith('##'):
                    line = f.readline().strip()
                if not line.startswith('#CHROM'):
                    eprint('weird vcf, quit')
                    # sorry
                    quit(2)
                return {sample: i for i,sample in enumerate(line.split('\t'))}

        def count_alt_alleles(gt):
            if gt == './.':
                return -1
            s = gt.split('/')
            return int(s[0]) + int(s[1])

        def concordance(vcf, sample_pairs):
            with open(sample_pairs, 'rt') as f:
                pairs = [[line.strip().split('\t')[0], line.strip().split('\t')[1]] for line in f]
            eprint(f'{len(pairs)} sample pairs read')
            sample_index = read_sample_index(vcf)
            eprint(f'{len(sample_index)-9} samples in vcf')

            for pair in pairs:
                if pair[0] not in sample_index:
                    eprint(f'sample {pair[0]} not found in vcf, quit')
                    quit(3)
                if pair[1] not in sample_index:
                    eprint(f'sample {pair[1]} not found in vcf, quit')
                    quit(3)

            with gzip.open(vcf, 'rt') as f:
                line = f.readline().strip()
                while line.startswith('##'):
                    line = f.readline().strip()
                print('\t'.join(['variant', 'af0', 'miss0', 'miss1',
                                 'ref0_ref1', 'ref0_het1', 'ref0_alt1',
                                 'het0_ref1', 'het0_het1', 'het0_alt1',
                                 'alt0_ref1', 'alt0_het1', 'alt0_alt1',
                                 'het0_het1_prop', 'hom0_hom1_prop']))
                eprint('printing counts to stdout')
                for line in f:
                    s = line.strip().split('\t')
                    # counts are: [0 ref_ref, 1 ref_het, 2 ref_alt, 3 het_ref, 4 het_het, 5 het_alt, 6 alt_ref, 7 alt_het, 8 alt_alt]
                    counts = [0 for i in range(0,9)]
                    miss0=miss1=0
                    for pair in pairs:
                        all0 = count_alt_alleles(s[sample_index[pair[0]]])
                        all1 = count_alt_alleles(s[sample_index[pair[1]]])
                        if all0 == -1:
                            miss0 = miss0 + 1
                        if all1 == -1:
                            miss1 = miss1 + 1
                        if all0 == -1 or all1 == -1:
                            continue
                        count_index = 3 * all0 + all1
                        counts[count_index] = counts[count_index] + 1

                    het_total = counts[1]+counts[3]+counts[4]+counts[5]+counts[7]
                    het_prop = counts[4]/het_total if het_total > 0 else -1

                    # calculate hom_total and hom_prop for decide minor allele
                    n_alt0 = 2*sum(counts[6:9]) + sum(counts[3:6])
                    n_ref0 = 2*len(pairs) - n_alt0
                    if n_alt0 < n_ref0:
                        hom_total = counts[2]+counts[5]+counts[6]+counts[7]+counts[8]
                        hom_prop = counts[8]/hom_total if hom_total > 0 else -1
                    else:
                        hom_total = counts[0]+counts[1]+counts[2]+counts[3]+counts[6]
                        hom_prop = counts[0]/hom_total if hom_total > 0 else -1

                    af = n_alt0/2/len(pairs)
                    print('\t'.join([s[2], str(af), str(miss0), str(miss1)] + [str(val) for val in counts + [het_prop, hom_prop]]))

            eprint('done')

        concordance('${vcf}', '${sample_pairs}')
        EOF

    >>>

    runtime {
        cpu: 2
        memory: "3G"
        disks: "local-disk 200 HDD"
        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.7"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 1
    }

    output {
        File conc = "concordance.gz"
    }
}

task concat_concordance {

    Array[File] conc_files

    command <<<

        zcat ${sep=" " conc_files} | awk 'NR==1||$0~"^chr"' | bgzip > concordance_allchr.gz

    >>>

    runtime {
        cpu: 2
        memory: "3G"
        disks: "local-disk 200 HDD"
        docker: "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.7"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 1
    }

    output {
        File conc_allchr = "concordance_allchr.gz"
    }

}

workflow variant_concordance {

    File vcf_vcf_loc
    Array[Array[String]] vcf_pairs = read_tsv(vcf_vcf_loc)
    String sample_pairs

    scatter (vcf_vcf in vcf_pairs) {
        call filter_merge {
            input: vcf1=vcf_vcf[0], vcf2=vcf_vcf[1], sample_pairs=sample_pairs
        }
        call concordance {
            input: vcf=filter_merge.out_vcf, sample_pairs=filter_merge.out_sample_pairs
        }
    }

    call concat_concordance {
        input: conc_files = concordance.conc
    }
}

