library(tidyverse)
library(broom)
library(magrittr)
library(parallel)
library(splines)
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")

compute_unexplained_variance <- function(probe_name) {
  if (is.na(snp_4_probe[[probe_name]])) {
    fm <- update(fm, . ~ . - snp)
    mf <- add_column(mf, y = meth[[probe_name]])
  } else {
    mf <- add_column(mf, snp = snp_matrix[, snp_4_probe[[probe_name]]], y = meth[[probe_name]])
  }
  m <- lm(fm, mf)
  R2 <- summary(m)$r.squared
  var_total <- var(mf$y, na.rm = TRUE)
  snp_coef <- tidy(m) %>% filter(term == "snp")
  p_value <- ifelse(is_empty(snp_coef$p.value), NA, snp_coef$p.value)
  mean <- ifelse(is_empty(snp_coef$estimate), NA, snp_coef$estimate)
  tibble(Probe = probe_name, 
         SNP = snp_4_probe[[probe_name]], 
         Var = var_total, 
         R2 = R2, 
         Unexplained_var = Var * (1 - R2), 
         Pvalue = p_value,
         Mean = mean)
}

meth <- get_m_values()
covs <- get_covs()
cis_res <- readRDS("./Data/RData/Results/MeQTL/cis.rds")
snp_matrix <- readRDS("./Data/RData/Genotypes/top_cis_snps.rds")[meth$SUBJID, ]

best_snp_per_probe <- cis_res %>%
  group_by(Probe) %>%
  summarize(best_snp = SNP[which.min(Pvalue)]) %>% 
  ungroup()

snp_4_probe <- setNames(best_snp_per_probe$best_snp, best_snp_per_probe$Probe)

control_fm <- get_control_fm()
mf <- get_model_frame(covs)
fm <- update(control_fm, y ~ . + snp)

unexplained_variance <- mclapply(best_snp_per_probe$Probe, 
                                 compute_unexplained_variance, 
                                 mc.cores = 12) %>% bind_rows()

probe_selection_20K <- unexplained_variance %>% 
  top_n(2e4, Unexplained_var) %>% 
  pull(Probe)

probe_selection_30K <- unexplained_variance %>% 
  top_n(3e4, Unexplained_var) %>% 
  pull(Probe)

probe_selection_50K <- unexplained_variance %>% 
  top_n(5e4, Unexplained_var) %>% 
  pull(Probe)

probe_selection_100K <- unexplained_variance %>% 
  top_n(1e5, Unexplained_var) %>% 
  pull(Probe)

saveRDS(unexplained_variance, "./Data/RData/Results/MeQTL/best_snp_per_probe.rds")
saveRDS(probe_selection_20K, "./Data/RData/Methylation/probes_trans_20K.rds")
saveRDS(probe_selection_30K, "./Data/RData/Methylation/probes_trans_30K.rds")
saveRDS(probe_selection_50K, "./Data/RData/Methylation/probes_trans_50K.rds")
saveRDS(probe_selection_100K, "./Data/RData/Methylation/probes_trans_100K.rds")
