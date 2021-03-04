# Process dbSNp vcf to numeric chromosomes

This wdl pipeline 1) maps a dbSNP variation VCF file from NCBI Assembly identifiers to numeric chromosomes and 2) limits sequences to only chromosomes 1-25 (i.e. 1-22, X(23), Y(24), MT(25)).

## How to use

### Choose the correct db SNP build & check chromosome mapping
dbSNP variation release information can be found ont he dbSNP site.
For example, b154 release announcement: https://www.ncbi.nlm.nih.gov/mailman/pipermail/dbsnp-announce/2020q2/000220.html  
the VCF file can be found in https://ftp.ncbi.nih.gov/snp/redesign/latest_release/
The file is https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.38.gz, i.e. the NCBI Assembly version (or similar) is GCF_000001405.38: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.38/  

### Create chromosome mapping file
Create a chromosome mapping file. The file should have one chromosome in each line, with the original name first, followed with a tab, and the new name last. E.g. `old_name\tnew_name\n`.
Check the chromosome mapping from the NCBI Assembly page for the assemblies used in that build.

### Run WDL
upload the chromosome map into a bucket. Add that path to the wdl json file, under map_vcf.vcf_task.chrom_map.
change the vcf file path to a correct one for map_vcf.vcf_task.vcf_dl_link. This is NOT a bucket link: This is a ftp link.  
For an example input json, see e.g. [wdl/dbsnp_hg38b153.json](wdl/dbsnp_hg38b153.json)
