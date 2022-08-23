
library(tidyverse)

# Load annotation files ---------------------------------------------------
  
roadmap_desc <- read_tsv("./Data/Roadmap_15_states/roadmap_descriptions.tsv")
roadmap_label_keys <- set_names(paste0(roadmap_desc$STATE_NU, "_", roadmap_desc$MNEMONIC), 
                                nm = roadmap_desc$SHORT_DESC)

meth_anno_roadmap <- get_anno_meth_roadmap()
meth_anno_geography <- get_anno_meth_geography()
meth_location <- get_anno_meth_location()
# snp_loc <- get_anno_snp_location() %>% 
#   dplyr::rename(SNP = SNP_ID)

geno <- get_snp_matrix(NULL)
genotype_map <- readRDS("./Data/RData/Genotypes/LabExMI_imputation_1000x5699237_annotated_map_with_ancestral.rds") %>% 
  select(SNP = SNP_ID, SNP_gene, Minor_allele, Major_allele)

get_maf <- function(x) {
  (2 * sum(x == 0, na.rm = TRUE) + sum(x == 1, na.rm = TRUE)) / sum(!is.na(x)) / 2
}
maf <- apply(geno, 2, get_maf)
keep_snps <- names(maf)[maf > 0.10]

# Load meQTL mapping results ----------------------------------------------

# trans <- readRDS("./Data/RData/Results/MeQTL/trans_snp_bottom_layer_sign_independent.rds") %>% 
#   select(Probe, SNP)
# 
# cis <- readRDS("./Data/RData/Results/MeQTL/cis_adjusted_probes.rds") %>% 
#   select(Probe, SNP)

# Load interaction results ------------------------------------------------

# res_env_crp <- readRDS("./Data/RData/Results/gene_environment_interaction_filtered_CRP.rds")
res_env <- readRDS("./Data/RData/Results/gene_environment_interaction_filtered_884.rds") %>%
  filter(SNP %in% keep_snps, !is.na(Standard_error)) %>% 
  left_join(meth_location) %>% 
  left_join(meth_anno_geography) %>% 
  left_join(meth_anno_roadmap) %>% 
  left_join(genotype_map)

res_cell <- readRDS("./Data/RData/Results/decomposition_lymphoid_genotype.rds") %>% 
  filter(SNP %in% keep_snps) %>% 
  rename(Variable = Cell) %>% 
  left_join(meth_location) %>% 
  left_join(meth_anno_geography) %>% 
  left_join(meth_anno_roadmap) %>% 
  left_join(genotype_map)

res <- bind_rows(res_env, res_cell)

res <- res %>% 
  mutate(Variable = fct_recode(factor(Variable),
                               CRP = "Log_CRP_levels",
                               CMV = "CMV_serostatusPositive", 
                               Smoking = "Smoking_statusSmoker", 
                               Sex = "SexFemale",
                               Age = "Age",
                               `B cells` = "X_CD19pos_of_CD45pos.panel5",
                               `CD4-CD8-` = "X_CD4negCD8neg_of_CD45pos.panel5",
                               `CD4` = "X_CD4pos_of_CD45pos.panel5",
                               `CD8` = "X_CD8bpos_of_CD45pos.panel5",
                               Monocytes = "X_mono_of_CD45pos.panel5",
                               `NK cells` = "X_NK_of_CD45pos.panel5")) %>% 
  mutate(Variable = fct_relevel(Variable, "CD4-CD8-", "CD8", "CD4", "B cells", "NK cells", "Monocytes", "CRP", "CMV", "Smoking", "Sex", "Age")) %>%
  filter(P_FDR < 0.05)

saveRDS(res, "./Data/RData/Results/genotype_interactions_of_CD45.rds")


# Construct interaction table ---------------------------------------------

res_age_int <- collect_ewas_results("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect/Age_interactions/", 3)
res_age_int_sign <- res_age_int %>% 
  filter(P_FDR < 0.05) %>% 
  rename(Variable1 = Interacting_term1, Variable2 = Interacting_term2) %>% 
  select(-Exposure, -Terms, -P_bonferroni) %>% 
  left_join(meth_location)

res_int <- res %>% 
  select(-contains("Chromatin"), -contains("Geography"), -contains("allele"), -SNP_distance_to_gene_bp) %>% 
  rename(Variable1 = Variable, Variable2 = SNP)

res_int <- bind_rows(res_int, res_age_int_sign)
saveRDS(res_int, "./Tables/interaction_table.rds")