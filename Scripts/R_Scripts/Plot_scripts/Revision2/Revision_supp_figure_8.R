

# Load packages -----------------------------------------------------------

library(tidyverse)
library(scales)
library(latex2exp)
library(patchwork)
source("./Scripts/R_scripts/Libraries/functions_for_plotting.R")

# Load data ---------------------------------------------------------------

genotype_map <- readRDS("./Data/RData/Genotypes/LabExMI_imputation_1000x5699237_annotated_map_with_ancestral.rds") %>%
  mutate(Which_allele_is_derived = if_else(Which_allele_is_ancestral == "Minor", "Major", "Minor"))

trans_adjust_on_snps_independent <- readRDS("./Data/RData/Results/MeQTL/trans_884_bonferroni_sign_independent.rds") %>%
  left_join(select(genotype_map, SNP_ID, Minor_allele, Major_allele, SNP_distance_to_gene_bp, Which_allele_is_derived), by = c("SNP" = "SNP_ID")) %>%
  filter(P_bonferroni < 0.05) %>%
  mutate(Mean = case_when(
    Which_allele_is_derived == "Minor" ~ - 2 * Mean,
    Which_allele_is_derived == "Major" ~ 2 * Mean,
    is.na(Which_allele_is_derived) ~ -2 * Mean
  ))

# roadmap_desc <- read_tsv("./Data/Roadmap_15_states/roadmap_descriptions.tsv")
# original_roadmap_levels <- paste0(roadmap_desc$STATE_NU, "_", roadmap_desc$MNEMONIC)

meth_anno_roadmap <- get_anno_meth_roadmap()
meth_anno_loc <- get_anno_meth_location()
meth_anno_geo <- get_anno_meth_geography()

cis_adjusted_probes <- readRDS("./Data/RData/Results/MeQTL/cis_adjusted_probes_884_bonferroni.rds") %>%
  filter(P_bonferroni < 0.05) %>%
  left_join(select(genotype_map, SNP_ID, Which_allele_is_derived), by = c("SNP" = "SNP_ID")) %>%
  mutate(Mean = case_when(
    Which_allele_is_derived == "Minor" ~ - 2 * Mean,
    Which_allele_is_derived == "Major" ~ 2 * Mean,
    is.na(Which_allele_is_derived) ~ -2 * Mean
  ))


# EWAS
ewas_age <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Age.rds") %>% filter(P_bonferroni < 0.05)
ewas_sex <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Sex.rds") %>% filter(P_bonferroni < 0.05)
ewas_smoking <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Smoking_status.rds") %>% filter(P_bonferroni < 0.05)
ewas_CMV <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/CMV_serostatus.rds") %>% filter(P_bonferroni < 0.05)
ewas_CRP <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Log_CRP_levels.rds") %>% filter(P_bonferroni < 0.05)

# Age dispersion
disp <- readRDS("./Data/RData/Results/Dispersion/res_dispersion_gamlss_884.rds") %>%
  filter(P_bonferroni < 0.05)

# A -----------------------------------------------------------------------

pltA <- plot_m_effects_without_counts(cis_adjusted_probes, 
                                      "Mean", 
                                      "P_bonferroni", 
                                      color = col_cis, 
                                      title = "Local meQTL",
                                      base_size = font_size_volcano, 
                                      point_size = point_size_volcano) +
  labs(tag = "A")

pltA


# B -----------------------------------------------------------------------

pltB <- plot_m_effects_without_counts(trans_adjust_on_snps_independent, 
                                      "Mean", 
                                      "P_bonferroni", 
                                      color = col_trans, 
                                      title = "Long range meQTL",
                                      base_size = font_size_volcano, 
                                      point_size = point_size_volcano) +
  labs(tag = "B")
pltB


# C -----------------------------------------------------------------------
pltC <- plot_m_effects_without_counts(ewas_age, 
                                      "Estimate", 
                                      "P_bonferroni", 
                                      color = col_age, 
                                      title = "Age",
                                      base_size = font_size_volcano, 
                                      point_size = point_size_volcano) +
  labs(tag = "C")

pltC


# D -----------------------------------------------------------------------

disp <- disp %>% 
  mutate(Estimate = Estimate / 50)
pltD <- plot_m_effects_without_counts(disp, 
                                      "Estimate", 
                                      "P_bonferroni", 
                                      color = col_age, 
                                      title = "Age dispersion",
                                      base_size = font_size_volcano, 
                                      point_size = point_size_volcano) +
  labs(tag = "D")

pltD

# E -----------------------------------------------------------------------

pltE <- plot_m_effects_without_counts(ewas_sex, 
                                      "Estimate", 
                                      "P_bonferroni", 
                                      color = col_sex, 
                                      title = "Sex",
                                      base_size = font_size_volcano, 
                                      point_size = point_size_volcano) +
  labs(tag = "E")
pltE


# F -----------------------------------------------------------------------

pltF <- plot_m_effects_without_counts(ewas_CMV, 
                                      "Estimate", 
                                      "P_bonferroni", 
                                      color = col_CMV, 
                                      title = "CMV",
                                      base_size = font_size_volcano, 
                                      point_size = point_size_volcano) +
  labs(tag = "F")
pltF

# G -----------------------------------------------------------------------

pltG <- plot_m_effects_without_counts(ewas_smoking, 
                                      "Estimate", 
                                      "P_bonferroni", 
                                      color = col_smoking, 
                                      title = "Smoking",
                                      base_size = font_size_volcano, 
                                      point_size = point_size_volcano) +
  labs(tag = "G")
pltG



# H -----------------------------------------------------------------------

pltH <- plot_m_effects_without_counts(ewas_CRP, 
                                      "Estimate", 
                                      "P_bonferroni", 
                                      color = col_CRP, 
                                      title = "CRP levels",
                                      base_size = font_size_volcano, 
                                      point_size = point_size_volcano) +
  labs(tag = "H")
pltH

# Assemble ----------------------------------------------------------------

plt <- (pltA | pltB | pltC | pltD) / (pltE | pltF | pltG | pltH)
ggsave("./Plots/Revision_plots/Supp_fig8/Sup_figure_8.pdf", plt, width = 200, height = 147, units = "mm")
ggsave("./Plots/Revision_plots/Supp_fig8/Sup_figure_8.png", plt, width = 200, height = 147, units = "mm", dpi = 500)


