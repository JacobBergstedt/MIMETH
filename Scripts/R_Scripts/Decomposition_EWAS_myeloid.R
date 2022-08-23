
library(sandwich)
library(tidyverse)
library(parallel)
library(lmtest)
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")

fit_example_int <- function(probe_beta, fm, mf) {
  
  mf$y <- probe_beta
  lm(fm, mf)
  
}

fit_int_myeloid <- function(probe_beta, fm, mf) {
  
  mf$y <- probe_beta
  m <- lm(fm, mf)
  res <- coeftest(m, vcov. = vcovHC)
  
  res[paste0(c("Age", "SexFemale", "CMV_serostatusPositive", "SmokerYes"), ":Myeloid"), c("Estimate", "Std. Error", "Pr(>|t|)")] %>% 
    as_tibble(rownames = "Term") %>% 
    mutate(Cell = "Myeloid",
           Variable = c("Age", "Sex", "CMV_serostatus", "Smoker")) %>% 
    select(-Term) %>% 
    rename(Standard_error = `Std. Error`,
           P_value = `Pr(>|t|)`)
  
}

fit_int <- function(probe_beta, fm, mf) {
  
  mf$y <- probe_beta
  m <- lm(fm, mf)
  res <- coeftest(m, vcov. = vcovHC)
  
  res[paste0(c("Age", "SexFemale", "CMV_serostatusPositive", "SmokerYes"), ":Log_lymph_mye_ratio"), c("Estimate", "Std. Error", "Pr(>|t|)")] %>% 
    as_tibble(rownames = "Term") %>% 
    mutate(Cell = "Log_lymph_mye_ratio",
           Variable = c("Age", "Sex", "CMV_serostatus", "Smoker")) %>% 
    select(-Term) %>% 
    rename(Standard_error = `Std. Error`,
           P_value = `Pr(>|t|)`)
  
}


n_cores <- 12
lineage_cells <- get_panel5_cells()

covs <- get_covs_884() %>% 
  mutate(Age = Age / 50) %>%
  dplyr::select(SUBJID, Age, Sex, CMV_serostatus, Smoking_status) %>% 
  inner_join(lineage_cells, by = "SUBJID") %>% 
  make_smoking_binary() %>% 
  mutate(Lymphoid = X_CD19pos_of_total.panel5 + X_NK_of_total.panel5 + X_CD8_of_total.panel5 + X_CD4_of_total.panel5 + X_CD4negCD8neg_of_total.panel5,
         Myeloid = X_PMN_of_total.panel5 + X_mono_of_total.panel5)

beta_values <- get_beta_values()
beta_values <- beta_values[match(covs$SUBJID, beta_values$SUBJID),]

fm <- y ~ (Age + Sex + CMV_serostatus + Smoker) * Myeloid
mf <- select(covs, Age, Sex, Smoker, CMV_serostatus, Myeloid)

res <- mclapply(beta_values[-1], fit_int_myeloid, fm = fm, mf = mf, mc.cores = n_cores) %>%
  bind_rows(.id = "Probe")

res <- full_join(res, get_anno_meth_location())
res <- res %>% 
  group_by(Variable) %>% 
  mutate(P_FDR = p.adjust(P_value, "BH"), P_bonferroni = p.adjust(P_value)) %>% 
  ungroup()

saveRDS(res, "./Data/RData/Results/EWAS/Beta_values/Environment/Sandwich/cell_decomposition_EWAS_Myeloid.rds")
