### R1 annotation
cluster start fganno --zone europe-west1-b --init gs://hail-common/vep/vep/GRCh38/vep85-GRCh38-init.sh --vep --version 0.1
cluster submit fganno scripts/annotate_variants.py --args 'gs://r1_data/scripts/fullfiles gs://finngen-production-library-green/R1_variant_annotation.tsv'


## R2 annotation
cluster start fganno --zone europe-west1-b --init gs://hail-common/vep/vep/GRCh38/vep85-GRCh38-init.sh --vep --version 0.1
cluster submit fganno scripts/annotate_variants.py --args 'gs://r2_data/annotation/r2_vcfs gs://r2_data/annotation/R2_var_annotation -vep-conf /vep/vep-gcloud-grch38.properties'

