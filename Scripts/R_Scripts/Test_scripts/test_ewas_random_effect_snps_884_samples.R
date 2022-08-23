#!/opt/gensoft/exe/R/3.6.2/scripts/Rscript
library(splines)
library(optparse)
library(dplyr)
library(parallel)
library(purrr)
library(lme4)
library(pbkrtest)
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")

# Parse options

cov <- opt$cov
n_cores <- opt$n_cores
chunk_path <- opt$chunk_path
methylation_path <- opt$methylation_path

cov <- "Age"
n_cores <- 12
methylation_path <- "./Data/RData/Methylation/MIMETH.minfi.final.MMatrix.autosomes.no_outliers.rds"
null_formula <- get_control_fm_15_cells()

p <- Sys.time()

fm <- null_formula %>%
  update(paste0("~ ", cov, " + .")) %>%
  rm_duplicate_age_terms()

# Load data
meth <- get_m_values()
snp_mat_list <- get_snp_mat_per_probe()

covariates <- get_covs_884() %>% 
  mutate(Age = Age / 50)

# Align IDs to new covariate IDs
snp_mat_list <- mclapply(snp_mat_list, function(x) x[covariates$SUBJID, , drop = FALSE], mc.cores = n_cores)
meth <- meth[match(covariates$SUBJID, meth$SUBJID),]
meth$SUBJID <- NULL
covariates <- covariates[c("Day_of_sampling", all.vars(fm))]

# Logic for special cases
if (cov == "Years_smoking") {
  
  fm <- update(fm, y ~ . - Smoking_status)
  ids <- !is.na(covariates$Smoking_status) & covariates$Smoking_status == "Smoker"
  covariates <- covariates[ids, ]
  covariates$Years_smoking <- covariates$Years_smoking / 20
  meth <- meth[ids, ]
  snp_mat_list <- map(snp_mat_list, ~ .[ids, , drop = FALSE])
  
} else if (cov == "Years_since_last_smoke") {
  
  fm <- update(fm, y ~ . - Smoking_status)
  ids <- !is.na(covariates$Smoking_status) & covariates$Smoking_status == "Ex_Smoker"
  covariates <- covariates[ids, ]
  covariates$Years_since_last_smoke <- covariates$Years_since_last_smoke / 20
  meth <- meth[ids, ]
  snp_mat_list <- map(snp_mat_list, ~ .[ids, , drop = FALSE])
  
} else if (cov == "CMV_serology") {
  
  fm <- update(fm, y ~ . - CMV_serostatus)
  ids <- !is.na(covariates$CMV_serostatus) & covariates$CMV_serostatus == "Positive"
  covariates <- covariates[ids, ]
  covariates$CMV_serology <- covariates$CMV_serology / 500
  meth <- meth[ids, ]
  snp_mat_list <- map(snp_mat_list, ~ .[ids, , drop = FALSE])
   
} else {
  
  fm <- update(fm, y ~ .)
  
}

terms_tib <- get_terms_tib(cov, covariates)



snp_mat_4_debug <- snp_mat_list[["cg05738196"]]
fit_ewas_debug_random_effect_snp(fm_add_snps_and_random_effect(fm, snp_mat_4_debug),
                                 na.omit(cbind(y = meth[["cg05738196"]], covariates, snp_mat_4_debug)),
                                 terms_tib)

gc()

res <- mcmapply(infer_ewas_random_effect_snps,
                meth[1:100],
                snp_mat_list[1:100],
                MoreArgs = list(fm = fm, covariates = covariates, terms_tib = terms_tib),
                SIMPLIFY = FALSE,
                mc.cores = n_cores) %>%
  bind_rows(.id = "Probe")

res <- compute_adjusted_p_values(res)
saveRDS(res, paste0(chunk_path, cov, ".rds"))
