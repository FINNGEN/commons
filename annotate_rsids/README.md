# Annotate a file with rsids from dbsnp reference
This script  is for annotating a file of variants with rsids. For example, annotating custom GWAS results with rsids. It uses a dbSNP rsid VCF, available in the green library, to get the correct rsid for each variant. The input file can be either uncompressed or (b)gzip-compressed. The output files are not compressed.  
The input file should be a tabular file separated by some specific character, for example tab or comma. The separator can be specified through the options. The output will have an additional column "rsid" added.  
Output can be piped to other command-line utilities, e.g. gzip or bgzip by leaving out the output file. In this case, output is written to standard output (to the terminal).
Usage:
```
usage: python3 rsidify.py [-h] --reference REFERENCE [--out OUT] -c CHROM -p POS -r REF
               -a ALT [-s SEP]
               file

Annotate a tsv/csv file with rsids from a dbSNP reference. Input file must be
sorted by chromosome and position.

positional arguments:
  file

optional arguments:
  -h, --help            show this help message and exit
  --reference REFERENCE
                        dbSNP reference VCF
  --out OUT             Output file path (not compressed). If output is not
                        specified, output will be written to stdout.
  -c CHROM, --chrom CHROM
                        Input file chromosome column
  -p POS, --pos POS     Input file position column
  -r REF, --ref REF     Input file reference allele column
  -a ALT, --alt ALT     Input file alternate allele column
  -s SEP, --sep SEP     File separator in input file


#Example: using tab as separator, writing output to uncompressed file
python3 rsidify.py FILE --reference DBSNP_REFERENCE_FILE -c CHROMOSOME_COLUMN -p POSITION_COLUMN -r REF_ALLELE_COLUMN -a ALT_ALLELE_COLUMN --out OUTPUT_FILE

#Example: comma as separator, piping the output to bgzip to get a bgzip-compressed output file
python3 rsidify.py FILE --reference DBSNP_REFERENCE_FILE -c CHROMOSOME_COLUMN -p POSITION_COLUMN -r REF_ALLELE_COLUMN -a ALT_ALLELE_COLUMN -s , |bgzip > OUTPUT_FILE.gz
```
