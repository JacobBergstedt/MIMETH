library(splines)
library(dplyr)
library(parallel)
library(purrr)
library(lme4)
library(pbkrtest)
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
p <- Sys.time()
chunk_path <- "./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_extended_IDOL/"

covs_of_interest <- c("Age", 
                      "CMV_serostatus", 
                      "Sex", 
                      "Heart_rate", 
                      "Smoking_status", 
                      "Temperature_ear", 
                      "CMV_serology", 
                      "Log_CRP_levels", 
                      "Hour_of_sampling")

covs_of_interest <- covs_of_interest[c(8:9)]

n_cores <- 8
covariates <- get_covs_884() %>% mutate(Age = Age / 50)
str <- names(covariates)
str <- str[grepl("IDOL_extended", str) & !grepl("Neu", str)]
fm_IDOL_extended <- paste0(str, collapse = " + ")
null_fm <- as.formula(paste0("~ Sex + Smoking_status + CMV_serostatus + ns(Age, df = 3) + Ancestry_PC1 + Ancestry_PC2 + ", fm_IDOL_extended))

# Load data
meth <- get_m_values()
snp_mat_list <- get_snp_mat_per_probe()

snp_mat_list <- mclapply(snp_mat_list, function(x) x[covariates$SUBJID, , drop = FALSE], mc.cores = n_cores)
meth <- meth[match(covariates$SUBJID, meth$SUBJID),]

for (cov in covs_of_interest) {
  
  res <- fit_EWAS_full(cov = cov, 
                       null_fm = null_fm, 
                       covariates = covariates, 
                       meth = meth, 
                       snp_mat_list = snp_mat_list, 
                       n_cores = n_cores)  
  
  saveRDS(res, paste0(chunk_path, cov, ".rds"))
  
}


print(Sys.time() - p)

