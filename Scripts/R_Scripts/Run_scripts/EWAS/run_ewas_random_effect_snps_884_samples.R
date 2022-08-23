#!/opt/gensoft/exe/R/3.6.2/scripts/Rscript
library(splines)
library(optparse)
library(dplyr)
library(parallel)
library(purrr)
library(lme4)
library(pbkrtest)
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")

# Parse options
option_list = list(
  make_option("--cov", action = "store", type = "character"),
  make_option("--n_cores", action = "store", type = "integer"),
  make_option("--null_formula", action = "store", type = "character"),
  make_option("--chunk_path", action = "store", type = "character")
)

opt <- parse_args(OptionParser(option_list=option_list))
cov <- opt$cov
n_cores <- opt$n_cores
chunk_path <- opt$chunk_path
methylation_path <- opt$methylation_path

p <- Sys.time()

fm <- as.formula(opt$null_formula) %>%
  update(paste0("~ ", cov, " + .")) %>%
  rm_duplicate_age_terms()

# Load data
meth <- get_m_values()
snp_mat_list <- get_snp_mat_per_probe()

covariates <- get_covs_884() %>% 
  mutate(Age = Age / 50)

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

cat("Formula:")
cat("\n")
print(fm)
cat("Variables in dataframe:")
cat("\n")
print(names(covariates))
cat("Numbers of cores in use:")
cat("\n")
print(n_cores)
cat("\n")
cat("Save results at:")
cat("\n")
cat(chunk_path)


snp_mat_4_debug <- snp_mat_list[["cg05738196"]]
fit_ewas_debug_random_effect_snp(fm_add_snps_and_random_effect(fm, snp_mat_4_debug),
                                 na.omit(cbind(y = meth[["cg05738196"]], covariates, snp_mat_4_debug)),
                                 terms_tib)

gc()

res <- mcmapply(infer_ewas_random_effect_snps,
                meth,
                snp_mat_list,
                MoreArgs = list(fm = fm, covariates = covariates, terms_tib = terms_tib),
                SIMPLIFY = FALSE,
                mc.cores = n_cores) %>%
  bind_rows(.id = "Probe")

res <- compute_adjusted_p_values(res)
saveRDS(res, paste0(chunk_path, cov, ".rds"))
