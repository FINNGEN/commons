# WDLs unrelated to any specific project

[wdl/plink_grm.wdl](wdl/plink_grm.wdl)
creates an LD-pruned genotype dataset we use for GRM calculation with SAIGE  
input: per-chromosome Plink .bed datasets (FinnGen release Plink files), list of high-quality variants to use, list of samples to be included (so we can subset e.g. population outliers out), max genotype missingess, MAF threshold, LD pruning parameters, output filename  
output: An LD-pruned Plink .bed dataset with high quality variants only

[filter_gwas.wdl](wdl/filter_gwas.wdl)  
filters a GWAS summary stat file by allele frequency  
input: a list of .gz summary stats, AF column name, min AF  
output: .gz summary stats filtered to AF > min AF

[filter_vcf_by_variants_to_bgen.wdl](wdl/filter_vcf_by_variants_to_bgen.wdl)  
filters VCF files to given variants and creates one BGEN file with those variants  
input: a list of .vcf.gz files, .gz list of variants  
output: one BGEN file + .bgi index with the given variants

[gnomad_extract_af.3.wdl](wdl/gnomad_extract_af.3.wdl)  
scrapes AF, FILTER, AN and rsid from gnomAD 3.0 VCF files per population  
input: a list of gnomAD 3.0 .vcf.gz files  
output: .gz tsv files + .tbi index per population with all gnomAD variants. Columns: cpra, AF_population, FILTER, total gnomAD AN, and rsid (rsid only for the "all" population)

[snpstats.wdl](wdl/snpstats.wdl)  
calculates qctool snpstats  
input: a list of bgen files  
output: a .gz tsv file + .tbi index with snpstats for all variants in the bgen files

[split_bgen.wdl](wdl/split_bgen.wdl)  
splits given bgen files (e.g. per chromosome) into smaller bgen files with the given number of variants in each file  
input: a list of bgen files  
output: smaller bgen files

[tomahawk_convert.wdl](wdl/tomahawk_convert.wdl)  
converts VCF files to Tomahawk .twk format and creates a mapping file between actual variant and its "Tomahawk position" (original Tomahawk doesn't handle multiallelics so we force all variants to be in different positions for Tomahawk)  
input: a list of .vcf.gz files  
output: one .twk file per .vcf file, a .gz mapping file + .tbi index
