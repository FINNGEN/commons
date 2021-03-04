# Process dbSNp vcf to numeric chromosomes

This wdl pipeline 1) maps a dbSNP variation VCF file from NCBI Assembly identifiers to numeric chromosomes and 2) limits sequences to only chromosomes 1-25 (i.e. 1-22, X(23), Y(24), MT(25)).

## How to use

### Choose the correct db SNP build & check chromosome mapping
dbSNP variation release information can be found on the dbSNP site.
For example, b154 release announcement: https://www.ncbi.nlm.nih.gov/mailman/pipermail/dbsnp-announce/2020q2/000220.html  
the VCF file can be found in https://ftp.ncbi.nih.gov/snp/redesign/latest_release/ or in the archive: https://ftp.ncbi.nih.gov/snp/redesign/archive/  
The Assembly information can be found on the filename: if the file is https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.38.gz, the assembly version is the filename. The NCBI Assembly version for that file is GCF_000001405.38: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.38/  
Note that you can use either files from hg38 or hg19 - this pipeline works for both, just make sure to use the correct chromosome mapping.
### Create chromosome mapping file
Create a chromosome mapping file from NCBI Assembly identifiers to chromosomes. The file should have one chromosome in each line, with the original name first, followed with a tab, and the new name last. E.g. `old_name\tnew_name\n`. The old names are the ncbi assembly identiers, new names are numbers from 1-25. Check the correct chromosome mapping from the NCBI Assembly page for the assemblies used in that build.

### Run WDL
upload the chromosome map into a bucket. Add that path to the wdl json file, under map_vcf.vcf_task.chrom_map.
change the vcf file path to a correct one for map_vcf.vcf_task.vcf_dl_link. This is NOT a bucket link: This is a ftp link to the dbsnp ftp server (or anywhere else the dbsnp vcf file is stored).  
For an example input json, see e.g. [wdl/dbsnp_hg38b153.json](wdl/dbsnp_hg38b153.json) or [wdl/dbsnp_hg38b154.json](wdl/dbsnp_hg38b154.json).
