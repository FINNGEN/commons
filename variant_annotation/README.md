## Scrape info from release VCF files workflow

Pipeline for scraping all info fields from released VCF files and joining with external (e.g. VEP) annotation.
Prepare configuration file by editing wdl/scrape.annot.wdl.json:
```
  "scrape_annots.vcfs": "gs://r7_data/R7_vcf_v0.txt" # vcf files in chromosomal order
  "scrape_annots.cpu": "2", how many cpus to use per node.
  "scrape_annots.memory": "8G", # how much memory to use per node
  "scrape_annots.external_annot": "gs://finngen_commons/annotations/R3_vep_annot.tsv.gz", ## optional external annotation. FIRST column MUST be called "variant" containing variant id in CHR:POS:REF:ALT format
  "scrape_annots.docker": "eu.gcr.io/finngen-refinery-dev/bioinformatics:0.6" # which docker to use
```

**IMPORTANT: Give the scrape_annots.vcfs files in chromosomal order as the chunks are joined to a single file in given order.**

Run wdl/scrape_annot.wdl with the prepared configuration file


## Annotation variants with VEP/Hail

First install hail `pip install hail` to install `hailctl` script.

Check/modify arguments in scripts/run_vep_annotate.sh for examples how to run variant annotations.


## fetching rsids & HGVS protein annotations for variants

use scripts/hgvsp_annotate.py

### Installation
```
python3 -m pip install -r hgvs_requirements.txt
```

### Usage

```
python3 scripts/hgsvp_annotate.py FILE COLUMN --out OUTPUT_FILENAME
```