<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Pipelines for FINNGEN data munging](#pipelines-for-finngen-data-munging)
	- [Pre-requirements](#pre-requirements)
	- [Scrape info from release VCF files workflow](#scrape-info-from-release-vcf-files-workflow)
		- [Copying files to end bucket in Google cloud](#copying-files-to-end-bucket-in-google-cloud)
	- [Convert to bgen WORKFLOW](#convert-to-bgen-workflow)
	- [Split bgen/vcf to chunks WORKFLOW](#split-bgenvcf-to-chunks-workflow)
		- [splitting by parallelizing by chunk.](#splitting-by-parallelizing-by-chunk)
			- [Run splitting](#run-splitting)
		- [splitting by parallelizing by multiple chunks (e.g. all chunks in chromosomes).](#splitting-by-parallelizing-by-multiple-chunks-eg-all-chunks-in-chromosomes)
	- [Copying files to end bucket in Google cloud](#copying-files-to-end-bucket-in-google-cloud)

<!-- /TOC -->


# Pipelines for FINNGEN data munging

## Pre-requirements

See docker/Dockerfile for instructions on building required tools if not using Docker.

Check latest available pre-built images
eu.gcr.io/finngen-refinery-dev/conv_dock:0.008
eu.gcr.io/finngen-refinery-dev/conv_dock:0.008

To enable docker edit options data/workflow.options.json to include wanted docker image and provide the file as -options to cromwell. Remove/edit the output location as seen fit.

data/backends.conf gives reasonable defaults for running locally/SGE/google cloud. Modify default backend to use on or the other.


## Scrape info from release VCF files workflow

Pipeline for scraping all info fields from released VCF files and joining with external annotation.
Prepare configuration file by editing scripts/scrape.annot.wdl.json:
```
  "scrape_annots.vcfs": "gs://r1_data/r1_files_to_scrape" # vcf files in chromosomal order
  "scrape_annots.cpu": "2", how many cpus to use per node.
  "scrape_annots.memory": "8G", # how much memory to use per node
  "scrape_annots.external_annot": "gs://r1_data/R1_variant_annotation.tsv", ## optional external annotation. FIRST column MUST be called "variant" containing variant id in CHR:POS:REF:ALT format
  "scrape_annots.docker": "eu.gcr.io/finngen-refinery-dev/conv_dock:0.01" # which docker to use
```

**IMPORTANT: Give the files in chromosomal order as the chunks are joined to a single file in given order.**

Run scripts/scrape_annot.wdl with the prepared configuration file

### Copying files to end bucket in Google cloud

After successful execution of the pipeline copy all files to desired bucket. Locate the [cromwell job id] and [cromwell output bucket] from Cromwell server or standalone output.

```
  gsutil -m cp gs://[cromwell output bucket]/scrape_annot/[cromwell job id]/call-join_annot/**/*.gz* [destination bucket]
```



## Convert to bgen WORKFLOW
Prepare a configuration file where each line is a full path to a file to be converted.

Edit scripts/convert_to_bgen.conf.json
and add path to the configuration file to "convert_to_bgen.files_to_conv"

Run WDL scripts/convert_to_bgen.wdl with scripts/convert_to_bgen.conf.json as inputs

```
java -Dconfig.file=data/backend.conf -jar cromwell.jar  run scripts/convert_to_bgen.wdl --inputs scripts/convert_to_bgen.json --options data/workflow.options.json
```

## Split bgen/vcf to chunks WORKFLOW
Edit or prepare configuration file like data/files.conf to include full path to each chromosome bgen/vcf file. You can add arbitrary chr names given after "chr_" but they must match the chromosomes given in the split files.

### splitting by parallelizing by chunk WORKFLOW
Prepare a configuration file for bgen splitting using the previous file config and pre-calculated chunk points (github.com/FINNGEN/chrsplit). Provided helper script example writes configuration of 5k variant chunks to data/splitting.conf.

If making configuration file by hand make sure that the file is tab separated and not whitespace. By Default in finngen we want to have to sets of chunks. 1) chunks for gene based analysis (not splitting boundaries) and 2) chunks for single variants with more balanced chunk sizes. You run the configuration file generation for both files and concatenate them:

```
scripts/generate_chunk_conf.py data/files.conf data/chrpos_gene_split_5k_chunks.txt splitting_gene_based.conf
scripts/generate_chunk_conf.py data/files.conf data/chrpos_additional_split_5k_chunks.txt splitting_variant_based.conf

cat splitting_gene_based.conf splitting_variant_based.conf > all_chunk_points.txt

```
#### Run splitting

Edit configuration file scripts/split_to_chunks.conf.json and add full path to config file created in previous step to  "split_to_chunks.input_blocks". Memory requested can be adjusted in the same file.

Run wdl scripts/split_to_chunks.wdl using scripts/split_to_chunks.conf.json as an input and using optional configurations. Provided docker file in workflow options contains all the necessary tools.

```
java -Dconfig.file=data/backend.conf -jar cromwell.jar  run scripts/split_to_chunks.wdl --inputs scripts/split_to_chunks.conf.json --options data/workflow.options.json
```


### Splitting by parallelizing by multiple chunks (e.g. all chunks in chromosomes) WORKFLOW

More efficient chunking can be achieved by first parallelizing to high cpu compute nodes per chromosome and then parallelizing the individual chunks within the node. Add a --splitype switch to generate configuration for this type of analysis

```
scripts/generate_chunk_conf.py data/files.conf data/chrpos_additional_split_5k_chunks.txt splitting_variant_based.conf --splitype by_chrom
```

Run wdl scripts/split_to_chunks_single_chrom.wdl using scripts/split_to_chunks_single_chrom.wdl.json as an input and using optional configurations. Provided docker file in workflow options contains all the necessary tools.

```
java -Dconfig.file=data/backend.conf -jar cromwell.jar  run scripts/split_to_chunks_single_chrom.wdl --inputs scripts/split_to_chunks_single_chrom.wdl --options data/workflow.options.json
```

## Copying files to end bucket in Google cloud

After successful execution of the pipeline copy all files to desired bucket. Locate the [cromwell job id] and [cromwell output bucket] from Cromwell server or standalone output.  Use file clobbing to copy files_to_conv

```
  gsutil -m cp gs://[cromwell output bucket]/split_to_chunks/[cromwell job id]/call-split_chunk/**/*.bgen [destination bucket]
```
