library(splines)
library(dplyr)
library(parallel)
library(purrr)
library(lme4)
library(pbkrtest)
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")

p <- Sys.time()
chunk_path <- "./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_6_props/"

covs_of_interest <- c("Age",
                      "CMV_serostatus",
                      "Sex",
                      "Heart_rate",
                      "Log_WBC_hematology",
                      "Smoking_status",
                      "Temperature_ear",
                      "CMV_serology",
                      "Log_CRP_levels",
                      "Hour_of_sampling")

n_cores <- 10
null_fm <- ~ Sex + Smoking_status + CMV_serostatus + ns(Age, df = 3) + Ancestry_PC1 + Ancestry_PC2 + 
  X_NK_of_total.panel5 + X_mono_of_total.panel5 + X_CD4_of_total.panel5 + X_CD8_of_total.panel5 + X_CD19pos_of_total.panel5

# Load data
panel5_cells <- get_panel5_cells() %>% select(-X_CD4negCD8neg_of_total.panel5)

rowsums <- rowSums(select(panel5_cells, -SUBJID))
panel5_cells <- mutate(panel5_cells, across(-SUBJID, ~ . / rowsums))
meth <- get_m_values()
snp_mat_list <- get_snp_mat_per_probe()

covariates <- get_covs_884() %>% 
  select(-X_mono_of_total.panel5) %>% 
  mutate(Age = Age / 50) %>% 
  inner_join(panel5_cells)

snp_mat_list <- mclapply(snp_mat_list, 
                         function(x) x[covariates$SUBJID, , drop = FALSE], mc.cores = n_cores)
meth <- meth[match(covariates$SUBJID, meth$SUBJID),]


for (cov in covs_of_interest) {

  res <- fit_EWAS_full(cov = cov,
                       null_fm = null_fm,
                       covariates = covariates,
                       meth = meth,
                       snp_mat_list = snp_mat_list,
                       n_cores = n_cores)
  
  print(cov)

  saveRDS(res, paste0(chunk_path, cov, ".rds"))

}


print(Sys.time() - p)

