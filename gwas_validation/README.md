# Validation

## GWAS Validation

This wdl is for visual inspection of GWAS results, compared to last release.  
It produces:
- list of variants with pval < 1e-5 in new release, and matching variants in old release
- list of variants that were GWS in old release
- list of variants that are "new loci", i.e. not GWS in old release and more than 1.5MB away from them.
For each of those lists, comparison plots of new release vs previous release beta, mlogp are produced. A regression of beta_new ~ 0 + beta_old and mlogp_new ~ 0 + mlogp_old is done, and it's R² is shown.

The wdl uses bioinformatics docker, since it has the necessary plotting, tabix and other utilities.

## Finemap validation

This wdl is for visual inspection of finemapping results. It concats susie credset summaries and calculates the number of good, bad and total cs per phenotype. Plots are drawn:
- Avg CS size between previous and current release
- Current & previous release size vs avg r2
- Current & previous release size vs min r2

### Usage
Fill in the following:
- `finemap__validation.current`: A tuple of current release tag, file with all credset susie summaries listed
- `finemap__validation.previous`: A tuple of previous release tag, file with all credset susie summaries listed
- `finemap_validation.docker`: docker, bioinformatics docker works fine.
- `finemap_validation.zones`: zones, europe-west1-b,c,d work fine.
Then submit the WDL.