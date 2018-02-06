# CONVERT_AND_SPLIT_VCF

Convert vcf to bgen and spli to chunks

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
java -jar cromwell.jar -Dconfig.file=data/backends.conf run scripts/convert_to_bgen.wdl --inputs scripts/convert_to_bgen.conf.json --options data/workflow.options.json
```

## Split bgen to chunks
Edit or prepare configuration file like data/files.conf to include full path to each chromosome bgen file. You can add arbitrary chr names given after "chr_" but they must match the chromosomes given in the split files.

Prepare a configuration file for bgen splitting using the previous file config and pre-calculated chunk points (github.com/FINNGEN/chrsplit). Provided helper script example writes configuration of 5k variant chunks to data/splitting.conf.

If making configuration file by hand make sure that the file is tab separated and not whitespace. By Default in finngen we want to have to sets of chunks. 1) chunks for gene based analysis (not splitting boundaries) and 2) chunks for single variants with more balanced chunk sizes. You run the configuration file generation for both files and concatenate them:

```
scripts/generate_chunk_conf.py data/files.conf data/chrpos_variant_split_5k_chunks.txt splitting_gene_based.conf
scripts/generate_chunk_conf.py data/files.conf data/chrpos_variant_split_5k_chunks.txt splitting_variant_based.conf

cat splitting_gene_based.conf splitting_variant_based.conf > all_chunk_points.txt

```

Edit configuration file scripts/split_to_chunks.conf.json and add full path to config file created in previous step to  "split_to_chunks.input_blocks"

Memory requested can be adjusted in the same file. Look

Run wdl scripts/split_to_chunks.wdl using scripts/split_to_chunks.conf.json as an input and using optional configurations.

```
java -jar cromwell.jar -Dconfig.file=data/backends.conf run scripts/split_to_chunks.wdl --inputs scripts/split_to_chunks.conf.json --options data/workflow.options.json
```
