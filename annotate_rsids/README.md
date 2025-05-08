# Annotate a file with rsids from dbsnp reference
This script  is for annotating a file of variants with rsids. For example, annotating custom GWAS results with rsids. It uses a dbSNP rsid VCF, available in the green library, to get the correct rsid for each variant. The input file can be either uncompressed or (b)gzip-compressed. The output files are not compressed.  
The input file should be a tabular file separated by some specific character, for example tab or comma. The separator can be specified through the options. The output will have an additional column "rsid" added.  
Output can be piped to other command-line utilities, e.g. gzip or bgzip by leaving out the output file. In this case, output is written to standard output (to the terminal).
The input file can either have chromosome, position, reference and alternate allele columns, or it can have a variant id column.

Usage:
```
usage: python3 rsidify.py[-h] --reference REFERENCE [--out OUT] (--cpra CPRA CPRA CPRA CPRA | --varid VARID VARID) [-s SEP] file

Annotate a tsv/csv file with rsids from a dbSNP reference. Input file must be sorted by chromosome and position.

positional arguments:
  file

options:
  -h, --help            show this help message and exit
  --reference REFERENCE
                        dbSNP reference VCF
  --out OUT             Output file path (not compressed). If output is not specified, output will be written to stdout.
  --cpra CPRA CPRA CPRA CPRA
                        chromosome, position, reference & alternate allele column names
  --varid VARID VARID   variant id and separator inside variant id, e.g. 'variant' ':'
  -s SEP, --sep SEP     File separator in input file


#Example: using tab as separator, writing output to uncompressed file
python3 rsidify.py FILE --reference DBSNP_REFERENCE_FILE -cpra CHROMOSOME_COLUMN  POSITION_COLUMN  REF_ALLELE_COLUMN ALT_ALLELE_COLUMN --out OUTPUT_FILE

#Example: comma as separator, piping the output to bgzip to get a bgzip-compressed output file
python3 rsidify.py FILE --reference DBSNP_REFERENCE_FILE -cpra CHROMOSOME_COLUMN  POSITION_COLUMN  REF_ALLELE_COLUMN ALT_ALLELE_COLUMN -s , |bgzip > OUTPUT_FILE.gz

#Example: variant id of CHR:POS:REF:ALT instead of chromosome,position,ref and alt allele columns
python3 rsidify.py FILE --reference DBSNP_REFERENCE_FILE --varid variant : --out OUTPUT_FILE
```
