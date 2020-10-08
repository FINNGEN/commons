# R4 ANALYSIS NOTES

## Get variants with info >0.95 in every batch

```
cd /mnt/disks/r4/geno
gsutil cp gs://fg-datateam-analysisteam-share/fromdatateam/R4-data/stats/R4_info_values.txt .
awk ‘NR>1 { m=$3; for (i=4;i<=NF;i++) if ($i<m) m=$i; if (m > 0.95) print $1 }’ R4_info_values.txt > R4_variants_info_allbatches_0.95.txt
```

## Impute sex

```
mkdir -p plink_missingness

gsutil -mq cp gs://fg-cromwell/convert_plink_merge/e3a8ac87-3f39-43a4-a52f-e01d70d9bfaf/call-chrom_convert/shard-22/results/R4_missingness_0.95_X/R4_missingness_0.95_X.* plink_missingness/

sudo gcloud docker -- pull gcr.io/finngen-refinery-dev/bioinformatics:0.3
sudo docker run -v`pwd`:/data -it gcr.io/finngen-refinery-dev/bioinformatics:0.3

cd /data/plink_missingness
# this used 10G mem - plink 1.9 requires way more than plink 2
time plink2 --memory 63000 --allow-extra-chr --snps-only --bfile R4_missingness_0.95_X --extract ../R4_variants_info_allbatches_0.95.txt --geno 0.03 --maf 0.05 --indep-pairwise 1000000 1000 0.1 --make-bed --out X_0.95
# real    4m15.888s

plink --allow-extra-chr --bfile X_0.95 --extract X_0.95.prune.in --impute-sex 0.4 0.8 --make-bed --out R4_seximputed

PLINK v1.90b6.7 64-bit (2 Dec 2018)            www.cog-genomics.org/plink/1.9/
(C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to R4_seximputed.log.                                                
Options in effect:                                                           
  --allow-extra-chr                                                          
  --bfile X_0.95                                                             
  --extract X_0.95.prune.in                                                                                                            
  --impute-sex 0.4 0.8                                                        
  --make-bed                                                                  
  --out R4_seximputed                                                        
                                                                             
65431 MB RAM detected; reserving 32715 MB for main workspace.                
105132 variants loaded from .bim file.                                
183694 people (0 males, 0 females, 183694 ambiguous) loaded from .fam.       
Ambiguous sex IDs written to R4_seximputed.nosex .                           
--extract: 812 variants remaining.                                           
Using 1 thread (no multithreaded calculations invoked).                      
Before main variant filters, 183694 founders and 0 nonfounders present.      
Calculating allele frequencies... done.                                      
Total genotyping rate is 0.988725.                                           
812 variants and 183694 people pass filters and QC.                          
Note: No phenotypes present.                                                   
--impute-sex: 812 Xchr and 0 Ychr variant(s) scanned, all sexes imputed. Report
written to R4_seximputed.sexcheck .                                                                                
--make-bed to R4_seximputed.bed + R4_seximputed.bim + R4_seximputed.fam ...  
done.                                                                                                            

# all sexes should be imputed - if not, to check which ones are not:
awk '{if ($4 != 1 && $4 != 2) print $0}' R4_seximputed.sexcheck
```

## Create covariate/phenotype file

```
gsutil -mq cp gs://thl-incoming-data/from_THL_registerteam/finngen_R4_phenotype_v2/finngen_R4_phenotype_data/*.gz .
gunzip -c finngen_R4_endpoint.gz | awk -f scripts/clean_pheno.awk | gzip > finngen_R4_endpoint_pheno_only.gz
```

Gather required files:
```
COHORT_FILE <- "../FGIDlist_DF3-DF4_COHORT-27-02-2019.txt"
BATCH_FILE <- "../../geno/fgfactory_R4_passSamples_10-07-2019.txt"
PCA_FILE <- "../../geno/PCA/R4_final.eigenvec"
PHENO_FILE <- "finngen_R4_endpoint_pheno_only.gz"
MINIMUM_FILE <- "finngen_R4_minimum.gz"
SEX_FILE <- "../R4_seximputed.sexcheck"
OUT_FILE <- "R4_COV_PHENO_V1.txt"
```

and run `Rscript scripts/create_r4_cov_pheno_file.R`, creates `R4_COV_PHENO_V1.txt.gz`

## Plink dataset for GRM

```
gsutil cp R4_variants_info_allbatches_0.95.txt gs://r4_data/grm/
gsutil ls gs://fg-cromwell/convert_plink_merge/e3a8ac87-3f39-43a4-a52f-e01d70d9bfaf/call-chrom_convert/shard-*/**/R4_missingness_0.95_*/*.bed | sort -V > R4_missingness_set_bedfiles.txt
gsutil cp R4_missingness_set_bedfiles.txt gs://r4_data/grm/

# get included samples after creating cov/pheno file with the R script
gunzip -c R4_COV_PHENO_V1.txt.gz | tail -n+2 | awk '{print $1,$1,0,0,0,-9}' > R4_COV_PHENO_V1.fam
gsutil cp R4_COV_PHENO_V1.fam gs://r4_data/grm/
```

GRM Plink dataset creation workflow `fb58aa0f-d03c-4914-8c1f-40ff7a73793b`

Inputs:

```
{
    "plink_grm.bedfilelist": "gs://r4_data/grm/R4_missingness_set_bedfiles.txt",
    "plink_grm.filter_prune.include_variants": "gs://r4_data/grm/R4_variants_info_allbatches_0.95.txt",
    "plink_grm.filter_prune.include_samples": "gs://r4_data/grm/R4_COV_PHENO_V1.fam",
    "plink_grm.filter_prune.geno_missing": 0.03,
    "plink_grm.filter_prune.maf": 0.01,
    "plink_grm.filter_prune.ld": "1000000 1000 0.1",
    "plink_grm.merge_plink.out": "R4_GRM_V1"
}
```

The WDL for GRM Plink creation is in [https://github.com/FINNGEN/commons/blob/master/wdl/plink_grm.wdl](https://github.com/FINNGEN/commons/blob/master/wdl/plink_grm.wdl)

```
gsutil cat gs://fg-cromwell/plink_grm/a22515bf-f716-4943-9005-ca93b790147f/call-merge_plink/R4_GRM_V0.bim | wc -l
# 58892 variants

gsutil -m cp gs://fg-cromwell/plink_grm/a22515bf-f716-4943-9005-ca93b790147f/call-merge_plink/R4_GRM_V0.* gs://r4_data/grm/
```

## SAIGE workflows

### V0

Demo phenos - most nulls ok (preemptibles not restarted) tests not (1G memory is not enough, 2G seems to be fine)
3fa51ecd-48c9-4cae-8c4e-d8ff1fe680e4

Demo phenos - rest of the nulls and all tests (update to cromwell 44 in between, tests run with n1-standard-1, and fixed SAIGE docker post munging)
006d578b-c88d-4630-88b8-828ec77b301d

## PheWeb import

### V0

```
gsutil -m cp gs://fg-cromwell/saige/006d578b-c88d-4630-88b8-828ec77b301d/**/*.pheweb.gz* gs://r4_data/saige/summary_stats_for_pheweb/
gsutil ls gs://r4_data/saige/summary_stats_for_pheweb/*.gz | sort > R4_DEMO_ETC_ENDPOINTS_summary_files.txt
gsutil cp R4_DEMO_ETC_ENDPOINTS_summary_files.txt gs://r4_data/pheno/
```

fdcd0274-bbab-43aa-8085-2afaeaeb3c6d  
^ some pheno tasks keep running, custom 1 cpu/3g  
e668e9d2-9cfc-47a1-985d-caf60217560f  
^ changed to n1-standard-1, all ok

```
cd /mnt/r4
mkdir -p generated-by-pheweb/pheno_gz generated-by-pheweb/manhattan generated-by-pheweb/sites generated-by-pheweb/qq generated-by-pheweb/cache

gsutil -mq cp gs://fg-cromwell/pheweb_import/e668e9d2-9cfc-47a1-985d-caf60217560f/call-pheno/**/*.gz* generated-by-pheweb/pheno_gz/
gsutil -mq cp gs://fg-cromwell/pheweb_import/fdcd0274-bbab-43aa-8085-2afaeaeb3c6d/call-pheno/**/*.gz* generated-by-pheweb/pheno_gz/
gsutil -mq cp gs://fg-cromwell/pheweb_import/e668e9d2-9cfc-47a1-985d-caf60217560f/call-pheno/**/manhattan/* generated-by-pheweb/manhattan/
gsutil -mq cp gs://fg-cromwell/pheweb_import/fdcd0274-bbab-43aa-8085-2afaeaeb3c6d/call-pheno/**/manhattan/* generated-by-pheweb/manhattan/
gsutil -mq cp gs://fg-cromwell/pheweb_import/e668e9d2-9cfc-47a1-985d-caf60217560f/call-pheno/**/qq/* generated-by-pheweb/qq/
gsutil -mq cp gs://fg-cromwell/pheweb_import/fdcd0274-bbab-43aa-8085-2afaeaeb3c6d/call-pheno/**/qq/* generated-by-pheweb/qq/

gsutil -mq cp gs://fg-cromwell/pheweb_import/e668e9d2-9cfc-47a1-985d-caf60217560f/**/top* generated-by-pheweb/
gsutil -mq cp gs://fg-cromwell/pheweb_import/e668e9d2-9cfc-47a1-985d-caf60217560f/**/best* generated-by-pheweb/
gsutil -mq cp gs://fg-cromwell/pheweb_import/e668e9d2-9cfc-47a1-985d-caf60217560f/**/pheno-list.json .

gsutil -mq cp gs://fg-cromwell/pheweb_import/fdcd0274-bbab-43aa-8085-2afaeaeb3c6d/**/sites/* generated-by-pheweb/sites/
gsutil -mq cp gs://fg-cromwell/pheweb_import/fdcd0274-bbab-43aa-8085-2afaeaeb3c6d/**/gene*b38* generated-by-pheweb/cache/

touch generated-by-pheweb/*/*.tbi
touch generated-by-pheweb/*.tbi
```

### Add metadata to pheno-list.json

```
mv pheno-list.json pheno-list.json.orig

cat ~/pheweb/js/phenolist_conf.json
{
    "phenolist": "/mnt/r4/pheno-list.json.orig",
    "pheno": "/home/jkarjala/R4_COV_PHENO_V0.txt",
    "phenoname": "/home/jkarjala/FINNGEN_ENDPOINTS_DF4_2019-06-27.txt",
    "category": "/home/jkarjala/FG_class_definitions.txt",
    "nomesco": "",
    "html": ""
}

node --max-old-space-size=8192 ~/pheweb/js/augmentPhenolist.js pheweb/js/phenolist_conf.json > /mnt/r4/pheno-list.json

cat pheno-list.log
info: 36 phenotypes, 3952 phenotypes' annotations and 14 categories read
info: Unexpected 'category': ALCOHOL, using 'Other'
info: Unexpected 'category': C, using 'Other'
info: No-category phenocode: IPF, using 'Other'
info: Unexpected 'category': RHEUMA, using 'Other'
info: Unexpected 'category': RX, using 'Other'
info: Unexpected 'category': RX, using 'Other'
info: No-category phenocode: SSRI, using 'Other'
info: No-category phenocode: T1D, using 'Other'
info: No-category phenocode: T2D, using 'Other'

cat descriptions_error.log
cat descriptions_combined.log
# empty
```

### FinnGen and gnomAD annotations

```
mkdir -p annotations/finngen && cd annotations/finngen
gsutil -m cp gs://r4_data/annotations/* .
mv R4_annotated_variants_v0.gz annotated_variants.gz
mv R4_annotated_variants_v0.gz.tbi annotated_variants.gz.tbi
touch *.tbi

cp -r /mnt/r3_1/annotations/gnomad ..
```

## Variant annotations

fcc15974-bf2c-4027-9f49-4880da9db25d

After workflow:

```
gsutil -m cp gs://fg-cromwell/scrape_annots/fcc15974-bf2c-4027-9f49-4880da9db25d/**/annotated_variants.gz* .
# get rid of extra columns
gunzip -c annotated_variants.gz | cut -f1-322,325-326,328- | bgzip -@9 > R4_annotated_variants_v0.gz && tabix -s2 -b3 -e3 R4_annotated_variants_v0.gz
gsutil -m cp R4_annotated_variants_v0.gz* gs://r4_data/annotations/
```

## Create variant lists with max 100k variants per chunk for smaller BGEN creation

```
gsutil -m cp gs://fg-cromwell/convert_plink_merge/e3a8ac87-3f39-43a4-a52f-e01d70d9bfaf/call-chrom_convert/**/R4_missingness_0.95_*.bim .
for file in *.bim; do split -d -l 100000 <(cut -f2 $file) temp/"${file/R4_missingness_0\.95/chr}_"; done
cd temp
for file in *; do mv "$file" "R4_${file/\.bim/}"; done
gsutil -m cp * gs://r4_data/bgen_split/
```

