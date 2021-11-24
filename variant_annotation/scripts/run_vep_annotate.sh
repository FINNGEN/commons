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

vcf_in="gs://finngen-imputation-panel/sisu4/vcf/snpid/*.vcf.gz"
out_prefix="gs://finngen-imputation-panel/sisu4/hail/sisu4_annot"

## submit annotation job
# add --no_maf in the end if VCFs dont have MAF INFO field
hailctl dataproc submit $clustername --pyfiles $(dirname -- "$0")/hail_functions.py \
	$(dirname -- "$0")/vep_annotate.py $vcf_in $out_prefix --overwrite
