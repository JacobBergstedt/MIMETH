## Packages

library(MASS)
library(tidyverse)
library(parallel)
library(splines)
library(sandwich)
source("./Scripts/R_scripts/Libraries/functions_for_mediation2.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
select <- dplyr::select
n_cores <- 10
n_sim <- 5000

## Data
cell_list_15_cells <- scan("./Data/prop_controls_15.txt", what = character())


meth <- get_m_values()
meth <- meth[match(covs$SUBJID, meth$SUBJID),]
meth$SUBJID <- NULL
covs$SUBJID <- NULL
varlist <- c("Age", "CMV", "Sex", "Heart_rate", "Smoking", "Years_smoking", "Years_since_last_smoke")
varlist <- varlist[6]
for (var in varlist) {
  
  treatment <- var
  
  if (var == "Age") {
    
    mf <- get_all_vars(get_control_fm_15_cells(), covs)
    
    fm_cell_responses <- map(cell_list_15_cells, 
                             ~ paste0(., " ~ Age + Sex + CMV_serostatus + Smoking_status + Ancestry_PC1 + Ancestry_PC2"))
    
    fm_total <- y ~ Age + Sex + CMV_serostatus + Smoking_status + Ancestry_PC1 + Ancestry_PC2
    fm_direct <- get_control_fm_15_cells() %>% 
      update(y ~ Age + .) %>% 
      rm_duplicate_age_terms()
    
    
  } else if (var == "Years_since_last_smoke") {
    
    mf <- get_all_vars(update(get_control_fm_15_cells(), . ~ . + Years_since_last_smoke), covs)
    keep <- mf$Smoking_status %in% c("Ex_Smoker", "Non_Smoker")
    mf <- mf[keep, ]
    meth_sub <- meth[keep, ]
    
    fm_cell_responses <- map(cell_list_15_cells, 
                             ~ paste0(., " ~ Years_since_last_smoke + Sex + CMV_serostatus + ns(Age, df = 3) + Ancestry_PC1 + Ancestry_PC2"))
    
    fm_total <- y ~ Years_since_last_smoke + Sex + CMV_serostatus + ns(Age, df = 3) + Ancestry_PC1 + Ancestry_PC2
    
    fm_direct <- get_control_fm_15_cells() %>% 
      update(y ~ Years_since_last_smoke + . - Smoking_status)
    
  } else if (var == "Years_smoking") {
    
    
    mf <- get_all_vars(update(get_control_fm_15_cells(), . ~ . + Years_smoking), covs)
    
    keep <- mf$Smoking_status == "Smoker"
    mf <- mf[keep, ]
    meth_sub <- meth[keep, ]
    
    fm_cell_responses <- map(cell_list_15_cells, 
                             ~ paste0(., " ~ Years_smoking + Sex + CMV_serostatus + ns(Age, df = 3) + Ancestry_PC1 + Ancestry_PC2"))
    fm_total <- y ~ Years_smoking + Sex + CMV_serostatus + ns(Age, df = 3) + Ancestry_PC1 + Ancestry_PC2
    
    fm_direct <- get_control_fm_15_cells() %>% 
      update(y ~ Years_smoking + . - Smoking_status)
    
  } else if (var == "Smoking") {
    
    mf <- covs %>%
      make_smoking_binary()

    mf <- get_all_vars(update(get_control_fm_15_cells(), ~ . - Smoking_status + Smoker), mf)

    fm_cell_responses <- map(cell_list_15_cells, ~ paste0(., " ~ ns(Age, df = 3) + Sex + CMV_serostatus + Smoker + Ancestry_PC1 + Ancestry_PC2"))
    fm_total <- y ~ ns(Age, df = 3) + Sex + CMV_serostatus + Smoker + Ancestry_PC1 + Ancestry_PC2
    fm_direct <- get_control_fm_15_cells() %>%
      update(y ~ . - Smoking_status + Smoker)

    treatment <- "SmokerYes"
    
  } else if (var == "Heart_rate") {
    
    mf <- get_all_vars(update(get_control_fm_15_cells(), ~ Heart_rate + .), covs)
    
    fm_cell_responses <- map(cell_list_15_cells, ~ paste0(., " ~ ns(Age, df = 3) + Heart_rate + Sex + CMV_serostatus + Smoking_status + Ancestry_PC1 + Ancestry_PC2"))
    fm_total <- y ~ ns(Age, df = 3) + Sex + Heart_rate + CMV_serostatus + Smoking_status + Ancestry_PC1 + Ancestry_PC2
    fm_direct <- get_control_fm_15_cells() %>% 
      update(y ~ Heart_rate + .)
    
    treatment <- "Heart_rate"
    
  } else {
    
    mf <- get_all_vars(get_control_fm_15_cells(), covs)
    
    fm_cell_responses <- map(cell_list_15_cells, ~ paste0(., " ~ ns(Age, df = 3) + Sex + CMV_serostatus + Smoking_status + Ancestry_PC1 + Ancestry_PC2"))
    fm_total <- y ~ ns(Age, df = 3) + Sex + CMV_serostatus + Smoking_status + Ancestry_PC1 + Ancestry_PC2
    fm_direct <- get_control_fm_15_cells() %>% update(y ~ .)
    treatment <- switch(var,
                        Sex = "SexFemale",
                        CMV = "CMV_serostatusPositive")
  }
  
  fm_list <- list(fm_cell_responses = fm_cell_responses, fm_total = fm_total, fm_direct = fm_direct)
  
  print(var)
  
  res <- run_mediation(meth_sub, 
                       fm_list = fm_list, 
                       mf = mf, 
                       treatment = treatment, 
                       n_sim = n_sim, 
                       cell_list = cell_list_15_cells,
                       n_cores = n_cores) %>%
    bind_rows(.id = "Probe") %>%
    add_column(Variable = var)
  
  saveRDS(res, paste0("./Data/Chunk_data/Results/Cell_mediation/m_values_884/", var, ".rds"))
}


