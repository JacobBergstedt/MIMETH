
library(tidyverse)
library(parallel)

# Load annotation files ---------------------------------------------------
  
roadmap_desc <- read_tsv("./Data/Roadmap_15_states/roadmap_descriptions.tsv")
roadmap_label_keys <- set_names(paste0(roadmap_desc$STATE_NU, "_", roadmap_desc$MNEMONIC), 
                                nm = roadmap_desc$SHORT_DESC)

meth_anno_roadmap <- get_anno_meth_roadmap()
meth_anno_geography <- get_anno_meth_geography()
meth_location <- get_anno_meth_location()
# snp_loc <- get_anno_snp_location() %>% 
#   dplyr::rename(SNP = SNP_ID)

# geno <- get_snp_matrix(NULL)
genotype_map <- readRDS("./Data/RData/Genotypes/LabExMI_imputation_1000x5699237_annotated_map_with_ancestral.rds") %>% 
  select(SNP = SNP_ID, SNP_gene, Minor_allele, Major_allele)

# get_maf <- function(x) {
#   (2 * sum(x == 0, na.rm = TRUE) + sum(x == 1, na.rm = TRUE)) / sum(!is.na(x)) / 2
# }
# maf <- apply(geno, 2, get_maf)
# keep_snps <- names(maf)[maf > 0.10]

keep_snps <- readRDS("./Data/RData/Genotypes/SNPs_to_keep_after_MAF_filter.rds")
  
  
# Load meQTL mapping results ----------------------------------------------

# trans <- readRDS("./Data/RData/Results/MeQTL/trans_snp_bottom_layer_sign_independent.rds") %>% 
#   select(Probe, SNP)
# 
# cis <- readRDS("./Data/RData/Results/MeQTL/cis_adjusted_probes.rds") %>% 
#   select(Probe, SNP)

# Load interaction results ------------------------------------------------

res_env_crp <- readRDS("./Data/RData/Results/gene_environment_interaction_sign_SNPs_CRP_884.rds") %>% 
  filter(SNP %in% keep_snps, !is.na(Standard_error)) %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  left_join(meth_location) %>% 
  left_join(meth_anno_geography) %>% 
  left_join(meth_anno_roadmap) %>% 
  left_join(genotype_map)
  
res_env <- readRDS("./Data/RData/Results/gene_environment_interaction_sign_SNPs_884.rds") %>%
  filter(SNP %in% keep_snps, !is.na(Standard_error)) %>%
  group_by(Variable) %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  ungroup() %>% 
  left_join(meth_location) %>% 
  left_join(meth_anno_geography) %>% 
  left_join(meth_anno_roadmap) %>% 
  left_join(genotype_map)

res_cell <- readRDS("./Data/RData/Results/decomposition_myeloid_genotype_sign_SNPs_beta_values_884.rds") %>%
  filter(SNP %in% keep_snps, !is.na(Standard_error)) %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  left_join(meth_location) %>% 
  left_join(meth_anno_geography) %>% 
  left_join(meth_anno_roadmap) %>% 
  left_join(genotype_map)

res <- bind_rows(res_env, res_cell, res_env_crp)
res <- res %>% 
  mutate(Variable = fct_recode(factor(Variable),
                               CRP = "Log_CRP_levels",
                               CMV = "CMV_serostatusPositive", 
                               Smoking = "Smoking_statusSmoker", 
                               Sex = "SexFemale",
                               Age = "Age")) %>% 
  mutate(Variable = fct_relevel(Variable, "Myeloid", "CRP", "CMV", "Smoking", "Sex", "Age")) %>%
  filter(P_bonferroni < 0.05)

saveRDS(res, "./Data/RData/Results/genotype_interactions_sign_SNPs_884.rds")


# EWAS x ratio ------------------------------------------------------------

res_EWAS_mye <- readRDS("./Data/RData/Results/EWAS/Beta_values/Environment/Sandwich/cell_decomposition_EWAS_Myeloid.rds") %>% 
  left_join(meth_location) %>% 
  left_join(meth_anno_geography) %>% 
  left_join(meth_anno_roadmap) %>% 
  filter(P_bonferroni < 0.05)

saveRDS(res_EWAS_mye, "./Tables/EWAS_cell_decomposition.rds")


# Construct interaction table ---------------------------------------------
res_age_int <- collect_ewas_results("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Age_interactions/", 3)
res_age_int_sign <- res_age_int %>% 
  filter(P_bonferroni < 0.05) %>% 
  rename(Variable1 = Interacting_term1, Variable2 = Interacting_term2) %>% 
  select(-Exposure, -Terms) %>% 
  left_join(meth_location)

res_int <- res %>% 
  select(-contains("Chromatin"), -contains("allele")) %>% 
  rename(Variable1 = Variable, Variable2 = SNP)

res_EWAS_decomposition <- res_EWAS_mye %>% 
  rename(Variable1 = Variable, Variable2 = Cell)

  
res_int <- bind_rows(res_int, res_age_int_sign, res_EWAS_decomposition) %>% select(-Geography)
saveRDS(res_int, "./Tables/interaction_table.rds")
write_tsv(res_int, "./Tables/interaction_table.tsv")