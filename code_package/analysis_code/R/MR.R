########### TwosampleMR ###########
######### version 0.6.15 ##########
# we used MR for retina eQTL and POAG MTAG for example code
# more detailed instruction could be found in website:"https://mrcieu.github.io/TwoSampleMR/"
# the eQTL data could be downloaded from supplementary table of paper:"DOI: 10.1038/s41588-019-0351-9"

# install packages
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)
library(data.table)
### 1.extract exposure
cis_retina <- "/MR/retina_sig_eqtl.txt"

MTAG <- "POAG_MTAG.txt"

bfile <- "g1000"
plink <- "/miniconda3/envs/plink/bin/plink"

# cis_retina_eQTL data 
exposure_dat <- read_exposure_data(
  filename = cis_retina,
  sep = "\t",
  snp_col = "SNP",
  chr_col = "CHR",
  pos_col = "BP",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  eaf_col = "FRQ",
  beta_col = "BETA",
  se_col = "SE",
  pval_col = "P",
  samplesize_col = "N",
  gene_col = "gene",
  phenotype_col = "gene"
)

# clump instrument SNPs
exposure_dat_clump <- clump_data(
  exposure_dat,
  clump_kb = 1000,
  clump_r2 = 0.01,
  clump_p1 = 5e-5,
  clump_p2 = 1,
  pop = "EUR",
  bfile = bfile,
  plink_bin = plink
)

### 2.extract outcome SNPs
# ao <- available_outcomes()
# read outcome GWAS
outcome_dat <- read_outcome_data(
  snps = exposure_dat_clump$SNP,
  filename = MTAG,
  sep = "\t",
  snp_col = "SNP",
  chr_col = "CHR",
  pos_col = "BP",
  effect_allele_col = "A2",
  other_allele_col = "A1",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "FRQ",
  pval_col = "P",
  samplesize_col = "N"
)
outcome_dat$outcome <- "POAG"

### 3.hormonized
harmonise_dat <- harmonise_data(
  exposure_dat = exposure_dat_clump,
  outcome_dat = outcome_dat
)

save(harmonise_dat, file = "/home/wangh/project1_glaucoma/analysis2/MR/harmonise_dat.RData")

### 4.perform MR and sensitive analysis
# mr_method_list()

wald <- mr_singlesnp(harmonise_dat, single_method = "mr_wald_ratio")

fixed <- mr_singlesnp(harmonise_dat, single_method = "mr_meta_fixed")

ivw <- mr(harmonise_dat, method_list = "mr_ivw")

median <- mr(harmonise_dat, method_list = "mr_weighted_median")

het <- mr_heterogeneity(harmonise_dat)
pleio <- mr_pleiotropy_test(harmonise_dat)

dflist <- list(ivw=ivw,median=median,wald=wald,fixed=fixed,het=het,pleio=pleio)

save_directory <- "/home/wangh/project1_glaucoma/analysis2/MR/result/retina_"

for (name in names(dflist)) {
  df <- dflist[[name]]
  file_path <- paste0(save_directory, name, ".txt")
  fwrite(df, file = file_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat("Saved:", file_path, "\n")
}

### result process ###

dir_path <- "/data01/user/wanghan/proj1/finemap/MR/result/TBC"
file_list <- list.files(dir_path, pattern = "\\.txt$", full.names = TRUE)

# read files
for (file in file_list) {
    
    fname <- basename(file)
    
    parts <- strsplit(fname, "_")[[1]]
    method <- gsub("\\.txt$", "", parts[length(parts)])

    df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    assign(method, df)

    cat("Loaded:", method, "\n")
}

### munge
pleio <- pleio %>% 
    mutate(model = "MR Egger") %>% 
    rename(pleio_p = pval)

fixed <- fixed[!is.na(fixed$p),] %>% 
    mutate(model = case_when(grepl("^All",SNP) ~ SNP, T ~ "wald"),
    SNP = case_when(grepl("^All",SNP) ~ "All", T ~ SNP)) %>% 
    select(exposure,outcome,model,SNP,b,se,p) %>% 
    group_by(model) %>% 
    mutate(FDR = p.adjust(p, method = "BH")) %>% 
    ungroup() %>% 
    mutate(model = gsub("^All - ", "", model)) %>% 
    rename(BETA = b, SE = se, P = p) %>% 
    left_join(het[,c(4,5,7,8)],by=c("exposure"="exposure","model"="method")) %>%
    left_join(pleio[,c(4,7,8)],by=c("exposure","model")) 

### save summary result, no FDR filter
fwrite(fixed, file = file.path(dir_path,"summary.txt"),
              sep = "\t", row.names = F, col.names = T)