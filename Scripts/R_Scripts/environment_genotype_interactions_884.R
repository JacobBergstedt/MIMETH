

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(glue)
library(parallel)
library(corrr)
library(tidyverse)
library(sandwich)
library(lmtest)
library(broom)
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
options(stringsAsFactors = FALSE)
options(width = 200)

# Functions ---------------------------------------------------------------
fit_example_main_effect_model <- function(meth_site, snps_site, mf) {
  
  db <- cbind(y = meth_site, mf, snps_site)
  
  fm_snps <- paste0(colnames(snps_site), collapse = " + ")
  fm_terms <- glue('Age + Sex + CMV_serostatus + Smoking_status + {fm_snps}')
  fm <- glue('y ~ Ancestry_PC1 + Ancestry_PC2 + {fm_cells} + {fm_terms}')
  
  lm(fm, db)
}

fit_example_interaction_model <- function(meth_site, snps_site, mf) {
  
  db <- cbind(y = meth_site, mf, snps_site)
  
  fm_snps <- paste0(colnames(snps_site), collapse = " + ")
  fm_interaction_terms <- glue('(Age + Sex + CMV_serostatus + Smoking_status) * ({fm_snps})')
  fm_interaction <- glue('y ~ Ancestry_PC1 + Ancestry_PC2 + {fm_cells} + {fm_interaction_terms}')
  
  lm(fm_interaction, db)
}

fit_example_CRP_interaction_model <- function(meth_site, snps_site, mf) {
  
  db <- cbind(y = meth_site, mf, snps_site)
  
  fm_snps <- paste0(colnames(snps_site), collapse = " + ")
  fm_interaction_terms <- glue('(Age + Sex + CMV_serostatus + Smoking_status + Log_CRP_levels) * ({fm_snps})')
  fm_interaction <- glue('y ~ Ancestry_PC1 + Ancestry_PC2 + {fm_cells} + {fm_interaction_terms}')
  
  lm(fm_interaction, db)
}


run_interaction_model <- function(meth_site, snps_site, mf) {
  
  db <- cbind(y = meth_site, mf, snps_site)
  
  fm_snps <- paste0(colnames(snps_site), collapse = " + ")
  fm_interaction_terms <- glue('(Age + Sex + CMV_serostatus + Smoking_status) * ({fm_snps})')
  fm_interaction <- glue('y ~ Ancestry_PC1 + Ancestry_PC2 + {fm_cells} + {fm_interaction_terms}')
  
  m <- lm(fm_interaction, db)
  
  terms <- expand.grid(c("Age", "SexFemale", "Smoking_statusSmoker", "CMV_serostatusPositive"), 
                       colnames(snps_site), stringsAsFactors = FALSE)
  terms_int <- paste0(terms$Var1, ":", terms$Var2)
  
  sum_m <- coeftest(m, vcov. = vcovHC)
  keep_ident <- terms_int %in% rownames(sum_m)
  terms_int <- terms_int[keep_ident]
  terms <- terms[keep_ident, ]
  
  tidy(sum_m) %>% 
    filter(term %in% terms_int) %>%
    rename(Estimate = estimate, Standard_error = std.error, P_value = p.value) %>% 
    add_column(Variable = terms$Var1, SNP = terms$Var2, .before = 1)
}

run_CRP_interaction_model <- function(meth_site, snps_site, mf) {
  
  db <- cbind(y = meth_site, mf, snps_site)
  
  fm_snps <- paste0(colnames(snps_site), collapse = " + ")
  fm_interaction_terms <- glue('(Age + Sex + CMV_serostatus + Smoking_status + Log_CRP_levels) * ({fm_snps})')
  fm_interaction <- glue('y ~ Ancestry_PC1 + Ancestry_PC2 + {fm_cells} + {fm_interaction_terms}')
  
  m <- lm(fm_interaction, db)
  
  terms <- expand.grid("Log_CRP_levels", 
                       colnames(snps_site), 
                       stringsAsFactors = FALSE)
  
  terms_int <- paste0(terms$Var1, ":", terms$Var2)
  
  sum_m <- coeftest(m, vcov. = vcovHC)
  keep_ident <- terms_int %in% rownames(sum_m)
  terms_int <- terms_int[keep_ident]
  terms <- terms[keep_ident, ]
  
  tidy(sum_m) %>% 
    filter(term %in% terms_int) %>%
    rename(Estimate = estimate, Standard_error = std.error, P_value = p.value) %>% 
    add_column(Variable = terms$Var1, SNP = terms$Var2, .before = 1)
}


# Main --------------------------------------------------------------------

n_cores <- 8
cell_list_15_cells <- get_15_props()
fm_cells <- paste0(cell_list_15_cells, collapse = " + ")

mf <- get_covs_884() %>% 
  dplyr::select(SUBJID, all_of(cell_list_15_cells), Sex, Smoking_status, Log_CRP_levels, Age, CMV_serostatus, Ancestry_PC1, Ancestry_PC2) %>% 
  mutate(Age = Age / 50, Smoking_status = fct_collapse(Smoking_status, Non_Smoker = c("Non_Smoker", "Ex_Smoker"))) %>%
  as.data.frame()

snps <- readRDS("./Data/RData/Genotypes/significant_snp_mat_per_probe.rds")
snps <- purrr::discard(snps, is_null)
  
snps <- mclapply(snps, function(x) x[mf$SUBJID, , drop = FALSE], mc.cores = n_cores)

meth <- get_m_values()
meth <- meth[match(mf$SUBJID, meth$SUBJID),]
meth <- meth[names(snps)]

res <- mcmapply(run_interaction_model, meth, snps, MoreArgs = list(mf = mf), SIMPLIFY = FALSE, mc.cores = n_cores) %>%
  bind_rows(.id = "Probe") %>%
  select(Probe, Variable, SNP, Estimate, Standard_error, P_value) %>%
  mutate(P_FDR = p.adjust(P_value, "BH"), P_bonferroni = p.adjust(P_value))

saveRDS(res, "./Data/RData/Results/gene_environment_interaction_sign_SNPs_884.rds")

res <- mcmapply(run_CRP_interaction_model, meth, snps, MoreArgs = list(mf = mf), SIMPLIFY = FALSE, mc.cores = n_cores) %>%
  bind_rows(.id = "Probe") %>%
  select(Probe, Variable, SNP, Estimate, Standard_error, P_value) %>%
  mutate(P_FDR = p.adjust(P_value, "BH"), P_bonferroni = p.adjust(P_value))

saveRDS(res, "./Data/RData/Results/gene_environment_interaction_sign_SNPs_CRP_884.rds")






  