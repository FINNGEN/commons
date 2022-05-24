## Scrape info from release VCF files workflow

Pipeline for scraping all info fields from released VCF files and joining with external (e.g. VEP) annotation.
Prepare configuration file by editing wdl/scrape.annot.wdl.json:
```
  "scrape_annots.vcfs": "gs://r7_data/R7_vcf_v0.txt" # vcf files in chromosomal order
  "scrape_annots.external_annot": "gs://finngen_commons/annotations/R3_vep_annot.tsv.gz", ## optional external annotation. FIRST column MUST be called "variant" containing variant id in CHR:POS:REF:ALT format
  "scrape_annots.docker": "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.6" # which docker to use
```

**NOTE: use the scrape_annots_v4.wdl for creating variant annotations for releases imputed with sisuv4.**

**IMPORTANT: Give the scrape_annots.vcfs files in chromosomal order as the chunks are joined to a single file in given order.**

Run wdl/scrape_annot.wdl with the prepared configuration file


## Annotation variants with VEP/Hail

First install hail `pip install hail` to install `hailctl` script.

Check/modify arguments in scripts/run_vep_annotate.sh for examples how to run variant annotations.


## Query OpenTargets

This script (scripts/query_opentargets.py) can be used to annotate a tab-separated file with left-aligned and minimized chr, pos, ref, alt in build 38 with OpenTargets phewas
and tag proxies.

Needs Python >= 3.7.2

Install requirements using `pip install -r variant_annotation/opentargets_requirements.txt`

Run scripts/query_opentargets.py for command line instructions

Example command with genomewide signficiant pvalue and r2 0.2 thresholds.  Additionally ignores FINNGEN_R5 study from results and searches all credible set variants in Autoreporting format.

`scripts/query_opentargets.py  input.tsv --p_threshold 5e-8 --ignore_studies FINNGEN_R5  --r_threshold 0.2 --parse_credsets_col credible_set_variants > input.with.OpenTargets.tsv
