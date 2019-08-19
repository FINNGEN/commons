require(data.table)
require(tidyr)
require(dplyr)
require(tidylog)
require(tibble)
require(R.utils)

COHORT_FILE <- "../finngen_R4_cohort_information"
BATCH_FILE <- "../../geno/fgfactory_R4_passSamples_10-07-2019.txt"
PCA_FILE <- "../../geno/PCA/R4_final.eigenvec"
PHENO_FILE <- "finngen_R4_endpoint_pheno_only.gz"
MINIMUM_FILE <- "finngen_R4_minimum.gz"
SEX_FILE <- "../R4_seximputed.sexcheck"
OUT_FILE <- "R4_COV_PHENO_V1.txt"

# endpoints with less cases will be excluded
MIN_N <- 100
# first endpoint column in PHENO_FILE
FIRST_PHENO <- "DEATH"

inlier_samples <- fread(PCA_FILE) %>%
 rename(FINNGENID=`#FID`) %>%
 filter(grepl("^FG", FINNGENID)) %>%
 select(FINNGENID)

cohorts <- fread(COHORT_FILE, header=F) %>%
  rename(FINNGENID=V1, cohort=V2)

chips <- fread(BATCH_FILE, header=F) %>%
  separate(V1, into=c("batch", "n_var", "chip", "FINNGENID"), sep=":") %>%
  mutate(IS_AFFY = ifelse(grepl("^Axiom", chip), 1, 0))

if (endsWith(MINIMUM_FILE, ".gz")) {
  mindata <- fread(paste0("gunzip -c ", MINIMUM_FILE))
} else {
  mindata <- fread(MINIMUM_FILE)
}
mindata$BMI <- mindata$WEIGHT/mindata$HEIGHT/mindata$HEIGHT*10000
# mindata$IS_SMOKE <- NA
# mindata$IS_SMOKE[which(mindata$SMOKE2 == "no")] = 0
# mindata$IS_SMOKE[which(mindata$SMOKE2 == "yes")] = 1
# mindata$IS_SMOKE_OR_FORMER <- NA
# mindata$IS_SMOKE_OR_FORMER[which(mindata$SMOKE3 == "never")] = 0
# mindata$IS_SMOKE_OR_FORMER[which(mindata$SMOKE3 == "current" | mindata$SMOKE3 == "former")] = 1

sex <- fread(SEX_FILE) %>%
  rename(FINNGENID=IID) %>%
  mutate(SEX_IMPUTED = (SNPSEX - 1)) %>%
  select(FINNGENID, SEX_IMPUTED)

pca <- fread(PCA_FILE) %>%
  rename(FINNGENID=`#FID`) %>%
  select(-IID)

if (endsWith(PHENO_FILE, ".gz")) {
  pheno <- fread(paste0("gunzip -c ", PHENO_FILE))
} else {
  pheno <- fread(PHENO_FILE)
}

age <- data.table(cbind(FINNGENID=pheno$FINNGENID,
                        AGE_AT_DEATH_OR_NOW=ifelse(pheno$DEATH == 1, pheno$DEATH_AGE, round(pheno$BL_AGE + as.numeric(format(Sys.Date(), "%Y")) - pheno$BL_YEAR, 1))))

cov_pheno <- inlier_samples %>%
  left_join(age, by="FINNGENID") %>%
  left_join(chips, by="FINNGENID") %>%
  filter(!is.na(AGE_AT_DEATH_OR_NOW)) %>%
  mutate(AGE_AT_DEATH_OR_NOW = as.numeric(AGE_AT_DEATH_OR_NOW))
  
# assign batches removing one batch for non-collinearity
uniq_batch <- sort(unique(cov_pheno$batch))
for (batch in uniq_batch[1:length(uniq_batch) - 1]) {
  cov_pheno[[paste0("BATCH_", batch)]] <- ifelse(cov_pheno$batch == batch, 1, 0)
}

cov_pheno <- cov_pheno %>%
  left_join(mindata, by="FINNGENID") %>%
  left_join(cohorts, by="FINNGENID") %>%
  left_join(sex, by="FINNGENID") %>%
  left_join(pca, by="FINNGENID") %>%
  left_join(pheno %>% select(-BL_YEAR, -BL_AGE), by="FINNGENID") %>%
  mutate(BMI = as.numeric(BMI))

# restrict to endpoints with enough cases
first_pheno_index <- match(FIRST_PHENO, names(cov_pheno))[1]
cols <- c(rep(T, (first_pheno_index - 1)), colSums(cov_pheno[,first_pheno_index:length(cov_pheno)], na.rm=T) >= MIN_N)
cov_pheno <- cov_pheno[,cols]

# exclude individuals with different registry/imputed sex
cov_pheno <- cov_pheno %>% filter((SEX == "male" & SEX_IMPUTED == 0) | (SEX == "female" & SEX_IMPUTED == 1))

stopifnot(length(which(is.na(cov_pheno$batch))) == 0)
stopifnot(length(which(is.na(cov_pheno$SEX_IMPUTED))) == 0)
stopifnot(length(which(is.na(cov_pheno$PC1))) == 0)

fwrite(cov_pheno, sep="\t", quote=F, file=OUT_FILE, na="NA")
gzip(OUT_FILE)

### batch gwas

# old_illumina <- as.numeric(grepl("^Illumina|^lllumina|^Human|^Broad|^Psych", cov_pheno$chip))
# cov_pheno$GSA_VS_OTHER_ILLUMINA <- as.numeric(grepl("^GSA", cov_pheno$chip))
# cov_pheno$GSA_VS_OTHER_ILLUMINA <- ifelse(old_illumina == 0 & cov_pheno$GSA_VS_OTHER_ILLUMINA == 0, NA, cov_pheno$GSA_VS_OTHER_ILLUMINA)
# cov_pheno$SUPER1_VS_SUPER2 <- ifelse(grepl("^super_1", cov_pheno$batch), 1, ifelse(grepl("^super_2", cov_pheno$batch), 0, NA))
# cov_pheno$SUPER1_VS_OTHER_ILLUMINA <- ifelse(grepl("^super_1", cov_pheno$batch), 1, ifelse(old_illumina == 1, 0, NA))
# cov_pheno$SUPER2_VS_OTHER_ILLUMINA <- ifelse(grepl("^super_2", cov_pheno$batch), 1, ifelse(old_illumina == 1, 0, NA))
# fwrite(cov_pheno, sep="\t", quote=F, file=paste0(OUT_FILE, ".SUPER_BATCH.txt"))
# gzip(paste0(OUT_FILE, ".SUPER_BATCH.txt"))
