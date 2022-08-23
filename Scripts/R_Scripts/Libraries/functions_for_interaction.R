library(dplyr)
library(parallel)
library(purrr)
library(sandwich)
library(lmtest)
source("./Scripts/R_functions/functions4inference.R")

age_interaction <- function(methylation_type) {
  
  if (methylation_type == "beta") {
    var_fun <- vcovHC
    methylation_path <- "./Data/RData/meth_beta_values.rds"
    save_path <- "./Data/RData/Interactions/Beta_values/interactions_beta_values_age_sex.rds"
  } else if (methylation_type =="M_value") {
    var_fun <- NULL
    methylation_path <- "./Data/RData/meth.rds"
    save_path <- "./Data/RData/Interactions/M_values/interactions_M_values_age_sex.rds"
  }
  
  # Load methylation
  meth <- readRDS(methylation_path)
  
  covariates <- readRDS("./Data/RData/covariates.rds") %>% 
    right_join(meth["ID"], covariates, by = c("SUBJID" = "ID"))
  
  covariates[["Age"]] <- covariates[["Age"]] - mean(covariates[["Age"]], na.rm = TRUE)
  
  fm <- as.formula(paste0("y ~ ", paste0("facs_PC", 1:10, collapse = " + "), 
                          " + Age * Sex + Smoking + CMVPositiveSerology + MCHC"))
  
  
  interaction_terms <- paste0("Age:SexFemale")
  
  terms_tib <- tibble(Terms = interaction_terms, 
                      Exposure = interaction_terms,
                      Levels = "SexFemale")
  
  meth$ID <- NULL  
  res <- mclapply(meth,
         infer_ewas,
         fm = fm,
         covariates = covariates,
         var_fun = var_fun,
         terms_tib = terms_tib,
         mc.cores = 12) %>% 
    bind_rows(.id = "Probe")
  
  res$P_FDR <- p.adjust(res$P_value, "fdr")
  saveRDS(res, save_path)
}
