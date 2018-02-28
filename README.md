# CONVERT_AND_SPLIT_VCF

Convert vcf to bgen and split to chunks

## Pre-requirements

These tools are needed for running the conversions
qctool_rc8 (tested on commit 76dbb1f of beta branch)
bgenix (Tested on 14d4b62 commit on master branch)

See docker/Dockerfile for instructions on building bgenix and qctool if not using Docker.
us.gcr.io/team-boston/finngen_qctool_small:0.001 docker available with tested builds.

To enable docker edit options data/workflow.options.json to include wanted docker image and provide the file as -options to cromwell. Remove/edit the output location as seen fit.

data/backends.conf gives reasonable defaults for running locally/SGE/google cloud. Modify default backend to use on or the other.

## Convert to bgen
Prepare a configuration file where each line is a full path to a file to be converted.

Edit scripts/convert_to_bgen.conf.json
and add path to the configuration file to "convert_to_bgen.files_to_conv"

Run WDL scripts/convert_to_bgen.wdl with scripts/convert_to_bgen.conf.json as inputs

```
java -Dconfig.file=data/backend.conf -jar cromwell.jar  run scripts/convert_to_bgen.wdl --inputs scripts/convert_to_bgen.json --options data/workflow.options.json
```

## Split bgen/vcf to chunks
Edit or prepare configuration file like data/files.conf to include full path to each chromosome bgen/vcf file. You can add arbitrary chr names given after "chr_" but they must match the chromosomes given in the split files.

### splitting by parallelizing by chunk.
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


### splitting by parallelizing by multiple chunks (e.g. all chunks in chromosomes).

More efficient chunking can be achieved by first parallelizing to high cpu compute nodes per chromosome and then parallelizing the individual chunks within the node. Add a --splitype switch to generate configuration for this type of analysis

```
scripts/generate_chunk_conf.py data/files.conf data/chrpos_additional_split_5k_chunks.txt splitting_variant_based.conf --splitype by_chrom
```

Run wdl scripts/split_to_chunks_single_chrom.wdl using scripts/split_to_chunks_single_chrom.wdl as an input and using optional configurations. Provided docker file in workflow options contains all the necessary tools.

```
java -Dconfig.file=data/backend.conf -jar cromwell.jar  run scripts/split_to_chunks_single_chrom.wdl --inputs scripts/split_to_chunks_single_chrom.wdl --options data/workflow.options.json
```

## Copying files to end bucket in Google cloud

After successful execution of the pipeline copy all files to desired bucket. Locate the [cromwell job id] and [cromwell output bucket] from Cromwell server or standalone output.  Use file clobbing to copy files_to_conv

```
  gsutil -m cp gs://[cromwell output bucket]/split_to_chunks/[cromwell job id]/call-split_chunk/**/*.bgen [destination bucket]
```
