
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(glue)
library(parallel)
library(tidyverse)
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
# source("./Scripts/R_scripts/Libraries/functions_for_meQTL_mapping.R")
options(stringsAsFactors = FALSE)
options(width = 200)

# Functions ---------------------------------------------------------------

generate_cv_indices <- function() {
  ids <- get_meth_ids()
  sample.int(5, length(ids), replace = TRUE)
}

compute_R2 <- function(obs, pred) {
  cor(obs, pred) ^ 2
}

estimate_genetic_interaction_P_value <- function(meth_site, snps_site, mf, cv_ids) {
  
  db <- cbind(y = meth_site, mf, snps_site)
  fm <- paste0("y ~ Age + Sex + CMV_serostatus + Smoking_status + Log_CRP_levels + ", paste0(colnames(snps_site), collapse = " + "))
  fm_interaction <- paste0("y ~ (Age + Sex + CMV_serostatus + Smoking_status + Log_CRP_levels) * (", paste0(colnames(snps_site), collapse = " + "), ")")
  m <- lm(fm, db)
  m_interaction <- lm(fm_interaction, db)
  as.data.frame(anova(m, m_interaction))[[2, "Pr(>F)"]]
  
}

estimate_genetic_interaction_R2 <- function(meth_site, snps_site, mf, cv_ids) {
  
  db <- cbind(y = meth_site, mf, snps_site)
  fm <- paste0("y ~ Age + Sex + CMV_serostatus + Smoking_status + Log_CRP_levels + ", paste0(colnames(snps_site), collapse = " + "))
  fm_interaction <- paste0("y ~ (Age + Sex + CMV_serostatus + Smoking_status + Log_CRP_levels) * (", paste0(colnames(snps_site), collapse = " + "), ")")
  res <- double(5)
  
  for (fold in 1:5) {
    
    db_train <- na.omit(db[cv_ids != fold, ])
    db_test <- na.omit(db[cv_ids == fold, ])
    m <- lm(fm, db_train)
    m_interaction <- lm(fm_interaction, db_train)
    pred_m <- predict(m, db_test)
    pred_m_interaction <- predict(m_interaction, db_test)
    R2_m <- compute_R2(pred_m, db_test$y)
    R2_m_interaction <- compute_R2(pred_m_interaction, db_test$y) 
    res[fold] <- R2_m_interaction - R2_m
  }
  
  tibble(R2 = res, fold = 1:5)
  
}

# Main --------------------------------------------------------------------
n_cores <- 6
cell_list_15_cells <- get_15_props()
mf <- get_covs_884() %>% 
  dplyr::select(SUBJID, all_of(cell_list_15_cells), Sex, Smoking_status, Log_CRP_levels, Age, CMV_serostatus) %>% 
  mutate(Age = Age / 50, Smoking_status = fct_collapse(Smoking_status, Non_Smoker = c("Non_Smoker", "Ex_Smoker"))) %>%
  as.data.frame()

snps <- readRDS("./Data/RData/Genotypes/snp_mat_per_probe.rds")
snps <- purrr::discard(snps, is_null)

snps <- mclapply(snps, 
                 function(x) x[mf$SUBJID, , drop = FALSE], mc.cores = n_cores)

meth <- get_beta_values()
meth <- meth[match(mf$SUBJID, meth$SUBJID),]
meth <- meth[names(snps)]

cv_ids <- generate_cv_indices()
R2_1 <- mcmapply(estimate_genetic_interaction_R2, meth, snps, MoreArgs = list(mf = mf, cv_ids = cv_ids), SIMPLIFY = FALSE, mc.cores = 1) %>% 
  bind_rows(.id = "Probe")

gc()

cv_ids <- generate_cv_indices()
R2_2 <- mcmapply(estimate_genetic_interaction_R2, meth, snps, MoreArgs = list(mf = mf, cv_ids = cv_ids), SIMPLIFY = FALSE, mc.cores = n_cores) %>% 
  bind_rows(.id = "Probe") %>% 
  mutate(fold = fold + 5)

gc()

cv_ids <- generate_cv_indices()
R2_3 <- mcmapply(estimate_genetic_interaction_R2, meth, snps, MoreArgs = list(mf = mf, cv_ids = cv_ids), SIMPLIFY = FALSE, mc.cores = 12) %>% 
  bind_rows(.id = "Probe") %>% 
  mutate(fold = fold + 10)

gc()

cv_ids <- generate_cv_indices()
R2_4 <- mcmapply(estimate_genetic_interaction_R2, meth, snps, MoreArgs = list(mf = mf, cv_ids = cv_ids), SIMPLIFY = FALSE, mc.cores = 12) %>% 
  bind_rows(.id = "Probe") %>% 
  mutate(fold = fold + 15)


R2 <- bind_rows(R2_1, R2_2, R2_3, R2_4)
saveRDS(R2, "./Data/RData/Results/Proportion_of_variance_explained/prop_var_explained_genetic_interactions_884.rds")

p_values <- mcmapply(estimate_genetic_interaction_P_value, meth, snps, 
                     MoreArgs = list(mf = mf), SIMPLIFY = FALSE, mc.cores = 12)

p_values <- tibble(Probe = names(p_values), P = unlist(p_values), P_FDR = p.adjust(unlist(p_values), "BH"))
saveRDS(p_values, "./Data/RData/Results/Proportion_of_variance_explained/prop_var_explained_genetic_interactions_P_value_884.rds")


