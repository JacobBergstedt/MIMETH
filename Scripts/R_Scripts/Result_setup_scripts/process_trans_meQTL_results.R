library(tidyverse)
library(parallel)
library(broom)
source("./Scripts/R_scripts/Libraries/functions_for_collections.R")
source("./Scripts/R_scripts/Libraries/functions_for_meQTL_mapping2.R")

n_cores <- 11
meth <- get_m_values()
geno <- get_snp_matrix()[meth$SUBJID, ]


paths <- paste0("./Data/Chunk_data/Results/MeQTL/Trans_meQTL/M_values/Normalized_cells_884/chromosome_", 1:22, ".rds")
trans <- collect_meQTL_results(paths, n_cores, threshold = 1, add_annotation = FALSE)




# BH ----------------------------------------------------------------------


trans_snp_bottom_layer <- two_stage_adjust(trans, 
                                           top_level = "Probe", 
                                           level = 0.05, 
                                           bottom_level = "SNP", 
                                           Bonferroni_bottom_adjustment = TRUE, 
                                           manual_local_adjustment_factor = 5699237)

trans_snp_bottom_layer_sign <- trans_snp_bottom_layer %>% 
  filter(Locally_significant)

# trans_probe_bottom_layer <- two_stage_adjust(trans, top_level = "SNP", level = 0.05, bottom_level = "Probe", manual_top_adjustment_factor = 5699237,
#                                              manual_local_adjustment_factor = 50000, Bonferroni_bottom_adjustment = FALSE)


meth_sign <- meth[trans_snp_bottom_layer_sign$Probe]
cis_sign <- readRDS("./Data/RData/Results/MeQTL/cis_adjusted_probes_884.rds") %>% 
  filter(Significant_family) %>% 
  select(Probe, Cis_SNP = SNP)

trans_with_cis <- left_join(trans_snp_bottom_layer_sign, cis_sign)
geno_sign <- geno[, unique(c(cis_sign$Cis_SNP, trans_with_cis$SNP))]
trans_snp_bottom_layer_sign_independent <- filter_snps_conditionally(trans_with_cis, meth_sign, geno_sign, p_tresh = 1e-6)

saveRDS(trans_snp_bottom_layer, "./Data/RData/Results/MeQTL/trans_snp_bottom_layer_884.rds")
saveRDS(trans_snp_bottom_layer_sign, "./Data/RData/Results/MeQTL/trans_snp_bottom_layer_sign_884.rds")
saveRDS(trans_snp_bottom_layer_sign_independent, "./Data/RData/Results/MeQTL/trans_snp_bottom_layer_sign_independent_884.rds")



# Bonferroni --------------------------------------------------------------



# trans_snp_bottom_layer <- two_stage_adjust(trans, 
#                                            top_level = "Probe", 
#                                            level = 0.05, 
#                                            bottom_level = "SNP", 
#                                            Bonferroni_bottom_adjustment = TRUE, 
#                                            manual_local_adjustment_factor = 5699237,
#                                            manual_top_adjustment_factor = 50000)


trans$P_bonferroni <- trans$Pvalue * 50000 * 5699237
trans_sign <- trans_snp_bottom_layer %>% 
  filter(P_bonferroni < 0.05)

# trans_probe_bottom_layer <- two_stage_adjust(trans, top_level = "SNP", level = 0.05, bottom_level = "Probe", manual_top_adjustment_factor = 5699237,
#                                              manual_local_adjustment_factor = 50000, Bonferroni_bottom_adjustment = FALSE)


meth_sign <- meth[trans_sign$Probe]
cis_sign <- readRDS("./Data/RData/Results/MeQTL/cis_adjusted_probes_884_bonferroni.rds") %>% 
  filter(P_bonferroni < 0.05) %>% 
  select(Probe, Cis_SNP = SNP)

trans_with_cis <- left_join(trans_sign, cis_sign)
geno_sign <- geno[, unique(c(cis_sign$Cis_SNP, trans_with_cis$SNP))]
trans_sign_independent <- filter_snps_conditionally(trans_with_cis, meth_sign, geno_sign, p_tresh = 1e-6)

saveRDS(trans, "./Data/RData/Results/MeQTL/trans_884_bonferroni.rds")
saveRDS(trans_sign, "./Data/RData/Results/MeQTL/trans_884_bonferroni_sign.rds")
saveRDS(trans_sign_independent, "./Data/RData/Results/MeQTL/trans_884_bonferroni_sign_independent.rds")