
library(tidyverse)
library(sandwich)
library(broom)
library(readxl)
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_compositions.R")
source("./Scripts/MiscFunctions.R")
library(robCompositions)
library(lmtest)
library(dirmult)

fit_compositional_model <- function(probe, snps, fm, cell_balances, covariates = NULL) {
  
  db <- cbind(y = probe, cell_balances)
  
  if (!is.null(covariates)) {
    
    fm <- update(fm, . ~ ns(Age, df = 3) + Sex + CMV_serostatus + Smoking_status + Ancestry_PC1 + Ancestry_PC2 + .)
    db <- cbind(db, covariates)
    
  }
  
  if (!is.null(snps)) {
    
    fm <- fm_add_snps(fm, snps)
    db <- cbind(db, snps)
    
  }
  
  
  m <- lm(fm, db)
  res <- coeftest(m, vcov. = vcovHC)
  
  
  tidy(res) %>% 
    filter(term %in% names(cell_balances)) %>% 
    select(Balance = term, Estimate = estimate, Standard_error = std.error, P_value = p.value) %>% 
    filter(Balance %in% c("Mye_Lym", "T_BNK", "CD4_CD8", "CD4Diff_CD4Naive", "CD8Diff_CD8Naive"))
  
}



to_m <- function(beta) {
  
  log2(beta / (1 - beta))
  
}

run_simulation <- function(X = NULL, change) {

  
  alpha <- sort(exp(rnorm(16, 0, 1)))
  alpha[cell_order] <- alpha
  
  
  if (is_null(X)) X <- rdirichlet(N, alpha = alpha)
  
  
  colnames(X) <- cell_list
  cell_balances <- balance_preds(X, sbp = read_xlsx("./Data/Cell_subset_partitions_3_with_DC.xlsx"))
  
  
  betas <- replicate(16, rbeta(N, 10, 1))

  colnames(betas) <- colnames(X)
  betas[, change] <- replicate(length(change), rbeta(N, 5, 1))
  y <- rowSums(X * betas)
  y_m <- to_m(y)
  
  
  fm <- as.formula(paste0("y ~ ", paste0(names(cell_balances), collapse = " + ")))

  res <- fit_compositional_model(y_m, NULL, fm = fm, cell_balances = cell_balances, covariates = NULL) %>% 
    rename(term = Balance)
  
  list(X = X, res_tib = res, change = change, betas = betas)
  

}

covs <- get_covs_884()
cell_list <- get_cell_list()
mu <-  covs[cell_list] %>% colMeans()


cell_order <- covs[cell_list] %>% colMeans() %>% order()


N <- nrow(covs)
# X <- covs[cell_list]


myeloid <- c("X_VIABLE_NEUTROPHILS_OF_TOTAL.panel7", "X_mono_of_total.panel5", "X_VIABLE_BASOPHILS_OF_TOTAL.panel7", "X_VIABLE_EOSINOPHILS_OF_TOTAL.panel7", "X_dendritic_cells.panel8")
lymphoid <- cell_list[!cell_list %in% myeloid]
CD4_diff <- c("X_CD4_EM_of_total.panel1", "X_CD4_EMRA_of_total.panel1", "X_CD4_CM_of_total.panel1")
CD4 <- c("X_CD4_EM_of_total.panel1", "X_CD4_EMRA_of_total.panel1", "X_CD4_naive_of_total.panel1", "X_CD4_CM_of_total.panel1")
CD8 <- c("X_CD8_EM_of_total.panel1", "X_CD8_EMRA_of_total.panel1", "X_CD8_naive_of_total.panel1", "X_CD8_CM_of_total.panel1")


si_mye <- map(1:100, ~ run_simulation(change = myeloid), .id = "Rep")
si_mye_tib <- map_dfr(si_mye, "res_tib")
si_CD4 <- map(1:100, ~ run_simulation(change = CD4), .id = "Rep")
si_CD4_tib <- map_dfr(si_CD4, "res_tib")
si_CD4_diff <- map(1:100, ~ run_simulation(change = CD4_diff), .id = "Rep")
si_CD4_diff_tib <- map_dfr(si_CD4_diff, "res_tib")

save(si_mye_tib, si_CD4_tib, si_CD4_diff_tib, file = "./Data/RData/composition_simulation_objects.RData")

