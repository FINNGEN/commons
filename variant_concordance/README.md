## Compute concordance of genotype calls between two datasets

This workflow computes basic variant concordance statistics between two given datasets using given matched samples

### Inputs:

`vcf_vcf_loc` is gs location of a file that contains on each line two vcf file gs locations tab-separated. Each vcf file can contain variants and samples that are not present in the other dataset but the stats are computed only for shared variants and based on sample pairs given in `sample_pairs`. If the datasets are small (e.g. just a handful of variants), this input can contain just one line with a pair of vcf files. If the datasets are large, each line can contain e.g. one chromosome

`sample_pairs` is gs location of a file that contains on each line two sample ids tab-separated. The first sample id refers to the sample id in the first file in `vcf_vcf_loc` and the second sample id refers to the sample id in the second file. Each line (pair) of these sample ids should be the same individiuals. Stats are computed using these sample pairs only

### Outputs:

`variant_concordance.concat_concordance.conc_allchr` is a gzipped tab-separated text file that contains the concordance stats, one variant per line. On the header line 0 refers to the first given vcf file and 1 refers to the second one. E.g. `het0_alt1` means the number of individuals that are heterozygous in the first file but alt homozygous in the second one. `het0_het1_prop` is `het0_het1` divided by `ref0_het1+het0_ref1+het0_het1+het0_alt1+alt0_het1`. `hom0_hom1_prop` is `alt0_alt1` divided by `ref0_alt1+het0_alt1+alt0_ref1+alt0_het1+alt0_alt1` if alt0 is the minor allele and `ref0_ref1` divided by `ref0_ref1+ref0_het1+ref0_alt1+het0_ref1+alt0_ref1` otherwise
