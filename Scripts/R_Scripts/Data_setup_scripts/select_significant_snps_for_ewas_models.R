#!/local/gensoft2/exe/R/3.6.0/scripts/Rscript

library(tidyverse)
library(corrr)
library(parallel)
library(dtplyr)
meth <- get_m_values()[-1]
cis <- readRDS("./Data/RData/Results/MeQTL/cis_adjusted_probes.rds") %>% 
  filter(Significant_family)

# Old results
# trans <- readRDS("./Data/RData/Results/MeQTL/trans_adjust_on_probes_significant_and_independent.rds")

trans <- readRDS("./Data/RData/Results/MeQTL/trans_snp_bottom_layer_sign_independent.rds")
snps <- readRDS("./Data/RData/Genotypes/LabExMI_imputation_1000x5699237_snp_matrix.rds")[get_meth_ids(), ]

snp_list_tib <- trans %>% 
  select(Probe, SNP) %>% 
  bind_rows(select(cis, Probe, SNP)) %>% 
  group_by(Probe) %>% 
  summarize(snp_list = list(SNP))

snps_per_probe <- snp_list_tib$snp_list
names(snps_per_probe) <- snp_list_tib$Probe
snp_mat_list <- mclapply(snps_per_probe, function(snp_name) snps[, snp_name, drop = FALSE], mc.cores = 2)
n <- names(meth)[!names(meth) %in% names(snp_mat_list)]
names(snp_mat_list[n]) <- n
snp_mat_list <- snp_mat_list[names(meth)]
saveRDS(snp_mat_list, "./Data/RData/Genotypes/significant_snp_mat_per_probe.rds")
