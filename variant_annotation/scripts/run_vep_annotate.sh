#!/usr/bin/env bash

clustername="fgveps"
build="GRCh38" #GRCh37
region="europe-west1"
zone="europe-west1-b"

project="finngen-refinery-dev"

network="projects/$project/regions/$region/subnetworks/default"

## can add more workers for full genomes. VEP is slow.
workers=2

## ## start hail cluster with vep options.
## requestes pays needed to read vep config file from hail buckets.
hailctl dataproc start $clustername --vep $build --region $region --requester-pays-allow-all --num-workers $workers --max-idle 30m --subnet $network --zone $zone

# note that double wildcards ** dont work here
vcf_in="gs://finngen-imputation-panel/sisu4/vcf/snpid/*.vcf.gz"
out_prefix="gs://finngen-imputation-panel/sisu4/hail/sisu4_annot"
#chip:
#vcf_in="gs://r9-data/chip/vcf/r9_axiom_chip_chr*_subset.vcf.gz"
#out_prefix="gs://r9-data/chip/anno/r9_chip_vep_annot_update"

## submit annotation job
#
# add --no_maf in the end if VCFs dont have MAF INFO field
#
# if you get IllegalArgumentException: requirement failed
# while running, uncomment hl._set_flags(no_whole_stage_codegen='1') in vep_annotate.py
#
hailctl dataproc submit $clustername --region $region --pyfiles $(dirname -- "$0")/hail_functions.py \
	$(dirname -- "$0")/vep_annotate.py $vcf_in $out_prefix --overwrite \
	--vep_conf gs://r9-data/vep95-GRCh38-loftee-gcloud.json
