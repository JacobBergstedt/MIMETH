
library(GenomicRanges)
library(rtracklayer)
library(fs)
library(tidyverse)
library(missMethyl)
library(vroom)
library(broom)
library(goseq)
select <- dplyr::select
library(scales)
source("./Scripts/R_scripts/Libraries/functions_for_enrichments.R")
source("./Scripts/R_scripts/Libraries/functions_for_plotting.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")

TF_meth_dep <- scan("./Data/meth_binding_TF.txt", what = character())
TF_CD8 <- c("TBX21", "ID3", "EOMES", "BCL6", "FOXO1", "STAT3", "ZEB1", "BACH1")

# Setup data structures ---------------------------------------------------


# path_to_chain = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
# lift_over_chain = import.chain(path_to_chain)
# 
# meth_location <- get_meth_location(fix = FALSE)
# meth_anno_roadmap <- get_meth_roadmap_annotation()
# 
# TFs <- list.files("./Data/TFBS_2019/pwm_tfbs_per_tf")
# TFBS <- map(TF. ~ add_TFBS(., meth_anno = meth_location, lift_over_chain = lift_over_chain)) %>%
#   reduce(inner_join)
# 
# TFBS <- bind_cols(TFBS["Probe"], discard(TFBS[-1], ~ min(table(.)) < 100)) %>% 
#   rename_with(.fn = ~ paste0("TF_", .), .cols = -Probe)
# 

# Get probe lists ---------------------------------------------------------

# dice <- get_dice() %>% 
#   pivot_longer(cols = -Gene, names_to = "Cell", values_to = "Expression")

meth_location <- get_anno_meth_location(fix = FALSE)
meth_anno_roadmap <- get_anno_meth_roadmap_all_cells() %>% 
  select(Probe, State = Chromatin_states_Mononuclear_cells)

meth_anno_loc_roadmap <- meth_location %>% 
  inner_join(meth_anno_roadmap)

meth_anno_TFBS <- readRDS("./Data/RData/Methylation/Annotation/meth_TFBS_annotation.rds") %>% 
  left_join(meth_anno_roadmap)


# EWAS probes -------------------------------------------------------------

ewas_age <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Age.rds") %>% 
  left_join(meth_anno_loc_roadmap) %>% 
  mutate(Probe_gene = simplify_probe_genes(Probe_gene))

ewas_sex <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Sex.rds") %>% 
  left_join(meth_anno_loc_roadmap) %>% 
  mutate(Probe_gene = simplify_probe_genes(Probe_gene))

ewas_CMV <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/CMV_serostatus.rds") %>% 
  left_join(meth_anno_loc_roadmap) %>% 
  mutate(Probe_gene = simplify_probe_genes(Probe_gene))

ewas_smoking <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Smoking_status.rds") %>% 
  left_join(meth_anno_loc_roadmap) %>% 
  filter(Levels != "Ex_Smoker") %>% 
  mutate(Probe_gene = simplify_probe_genes(Probe_gene))

thresh <- 0.05
age_pos_probes <- ewas_age %>% filter(P_bonferroni < thresh, Estimate > 0) %>% pull(Probe)
sex_pos_probes <- ewas_sex %>% filter(P_bonferroni < thresh, Estimate > 0) %>% pull(Probe)
CMV_pos_probes <- ewas_CMV %>% filter(P_bonferroni < thresh, Estimate > 0) %>% pull(Probe)
smoking_pos_probes <- ewas_smoking %>% filter(P_bonferroni < thresh, Estimate > 0) %>% pull(Probe)

age_neg_probes <- ewas_age %>% filter(P_bonferroni < thresh, Estimate < 0) %>% pull(Probe)
sex_neg_probes <- ewas_sex %>% filter(P_bonferroni < thresh, Estimate < 0) %>% pull(Probe)
CMV_neg_probes <- ewas_CMV %>% filter(P_bonferroni < thresh, Estimate < 0) %>% pull(Probe)
smoking_neg_probes <- ewas_smoking %>% filter(P_bonferroni < thresh, Estimate < 0) %>% pull(Probe)


# Mediation probes --------------------------------------------------------

ewas_med_age <- readRDS("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/Age.rds") %>% 
  filter(Effect == "Mediation") %>% 
  mutate(P_bonferroni = p.adjust(P_value))

ewas_med_sex <- readRDS("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/Sex.rds") %>% 
  filter(Effect == "Mediation") %>% 
  mutate(P_bonferroni = p.adjust(P_value))

ewas_med_smoking <- readRDS("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/Smoking.rds") %>% 
  filter(Effect == "Mediation") %>% 
  mutate(P_bonferroni = p.adjust(P_value))

ewas_med_CMV <- readRDS("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/CMV.rds") %>% 
  filter(Effect == "Mediation") %>% 
  mutate(P_bonferroni = p.adjust(P_value))

ewas_med_age_sign <- ewas_med_age %>% filter(P_bonferroni < thresh)
ewas_med_sex_sign <- ewas_med_sex %>% filter(P_bonferroni < thresh)
ewas_med_smoking_sign <- ewas_med_smoking %>% filter(P_bonferroni < thresh)
ewas_med_CMV_sign <- ewas_med_CMV %>% filter(P_bonferroni < thresh)

age_med_pos_probes <- ewas_med_age_sign %>% filter(Estimate > 0) %>% pull(Probe)
age_med_neg_probes <- ewas_med_age_sign %>% filter(Estimate < 0) %>% pull(Probe)
CMV_med_pos_probes <- ewas_med_CMV_sign %>% filter(Estimate > 0) %>% pull(Probe)
CMV_med_neg_probes <- ewas_med_CMV_sign %>% filter(Estimate < 0) %>% pull(Probe)

# Overall -----------------------------------------------------------------
# active_states <- c("TSS", "Fl. TSS")
# meth_anno_sub <- meth_anno_TFBS %>% filter(State %in% active_states)


TF_list <- select(meth_anno_TFBS, -Probe) %>% names()
TF_meth_dep_list <- TF_meth_dep[TF_meth_dep %in% TF_list]

res_direct     <- list(Age_Pos     = age_pos_probes,
                       Age_Neg     = age_neg_probes,
                       CMV_Pos     = CMV_pos_probes,
                       CMV_Neg     = CMV_neg_probes,
                       Sex_Pos     = sex_pos_probes,
                       Sex_Neg     = sex_neg_probes,
                       Smoking_Pos = smoking_pos_probes,
                       Smoking_Neg = smoking_neg_probes) %>% 
  mclapply(get_TFBS_enrichment, anno = meth_anno_TFBS, mc.cores = 8) %>% 
  bind_rows(.id = "Variable") %>% 
  separate(col = "Variable", into = c("Variable", "Direction"), sep = "_")

# Enrichments positive direction
hits <- res_direct  %>%
  filter(Variable == "Age", Direction == "Pos") %>%
  arrange(desc(conf.low)) %>%
  filter(P_FDR < 1e-5, estimate > 1) %>%
  pull(TF)

go_enrichment(hits, TF_list)

#
hits <- res_direct %>%
  filter(Variable == "CMV", Direction == "Pos") %>%
  arrange(desc(conf.low)) %>%
  filter(P_FDR < 1e-5, estimate > 1) %>%
  pull(TF)

p <- go_enrichment(hits, TF_list)

# Negative direction
hits <- res_direct %>%
  filter(Variable == "CMV", Direction == "Neg") %>%
  arrange(desc(conf.low)) %>%
  filter(P_FDR < 0.05, estimate > 1) %>%
  pull(TF)

p <- go_enrichment(hits, TF_list)


# Mediation ---------------------------------------------------------------
res_mediated     <- list(Age_Pos     = age_med_pos_probes,
                         Age_Neg     = age_med_neg_probes,
                         CMV_Pos     = CMV_med_pos_probes,
                         CMV_Neg     = CMV_med_neg_probes) %>% 
  map_dfr(get_TFBS_enrichment, anno = meth_anno_TFBS, .id = "Variable") %>% 
  separate(col = "Variable", into = c("Variable", "Direction"), sep = "_")


hits <- res_mediated  %>%
  filter(Variable == "Age", Direction == "Pos") %>%
  arrange(desc(conf.low)) %>%
  filter(P_FDR < 1e-5, estimate > 1) %>%
  pull(TF)

pos_age_med_enrichments <- go_enrichment(hits, TF_list = TF_list)

hits <- res_mediated  %>%
  filter(Variable == "CMV", Direction == "Pos") %>%
  arrange(desc(conf.low)) %>%
  filter(P_FDR < 1e-100, estimate > 1) %>%
  pull(TF)

# dice %>%
#   filter(Gene %in% hits) %>%
#   group_by(Gene) %>%
#   arrange(desc(Expression), .by_group = TRUE) %>%
#   dplyr::slice(1) %>%
#   ungroup() %>%
#   dplyr::count(Cell)

pos_CMV_med_enrichments <- go_enrichment(hits, TF_list = TF_list)

hits <- res_mediated  %>%
  filter(Variable == "CMV", Direction == "Neg") %>%
  arrange(desc(conf.low)) %>%
  filter(P_FDR < 1e-5, estimate > 1) %>%
  pull(TF)

neg_CMV_med_enrichments <- go_enrichment(hits, TF_list = TF_list)
saveRDS(res_direct, file = "./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_dir_EWAS_884_bonferroni.rds")
saveRDS(res_mediated, file = "./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_med_EWAS_884_bonferroni.rds")



# Outside of CGI ----------------------------------------------------------
meth_anno_geo <- get_anno_meth_geography()
meth_anno_TFBS_non_CGI <- meth_anno_TFBS %>% 
  left_join(meth_anno_geo) %>% 
  filter(Geography != "Island")

ewas_age_non_CGI <- left_join(meth_anno_TFBS_non_CGI, ewas_age)

age_pos_probes <- ewas_age %>% filter(P_bonferroni < thresh, Estimate > 0) %>% pull(Probe)
age_neg_probes <- ewas_age %>% filter(P_bonferroni < thresh, Estimate < 0) %>% pull(Probe)

res_direct_age_non_cgi     <- list(Age_Pos     = age_pos_probes,
                                   Age_Neg     = age_neg_probes) %>% 
  mclapply(get_TFBS_enrichment, anno = meth_anno_TFBS_non_CGI, mc.cores = 2) %>% 
  bind_rows(.id = "Variable") %>% 
  separate(col = "Variable", into = c("Variable", "Direction"), sep = "_")

saveRDS(res_direct_age_non_cgi, "./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_dir_age_non_CGI_884_bonferroni.rds")


# Dispersion --------------------------------------------------------------
# disp <- readRDS("./Data/RData/Results/Dispersion/res_dispersion_gamlss.rds") %>%
#   mutate(P_FDR = p.adjust(P_value, "BH"))
# 
# disp_pos_sites <- disp %>% filter(Estimate > 0, P_FDR < 0.05) %>% pull(Probe)
# 
# get_TFBS_enrichment(disp_pos_sites, meth_anno_TFBS)


# Plot enrichments --------------------------------------------------------


# Direct
for (var in unique(res_direct$Variable)) {
  for (dir in unique(res_direct$Direction)) {
    plt <- res_direct %>% 
      filter(Variable == var, Direction == dir) %>%
      top_n(n = 15, conf.low) %>% 
      arrange(desc(conf.low)) %>% 
      mutate(Term = factor(TF, TF)) %>% 
      dplyr::rename(Estimate = estimate, Low = conf.low, High = conf.high) %>% 
      plot_odds_ratio(NULL,
                      bar_size = odds_ratio_errorbar_size, 
                      color = col_age, 
                      base_size = font_size_odds_ratio)
    ggsave(paste0("./Plots/Odds_ratios_bonferroni_884/TFBS/EWAS/Direct/", var, "_", dir, ".pdf"), plt, width = 5, height = 5, units = "cm")      
  }
}



# Mediation
for (var in unique(res_mediated$Variable)) {
  for (dir in unique(res_mediated$Direction)) {
    plt <- res_mediated %>% 
      filter(Variable == var, Direction == dir) %>%
      top_n(n = 15, conf.low) %>% 
      arrange(desc(conf.low)) %>% 
      mutate(Term = factor(TF, TF)) %>% 
      dplyr::rename(Estimate = estimate, Low = conf.low, High = conf.high) %>% 
      plot_odds_ratio(NULL,
                      bar_size = odds_ratio_errorbar_size, 
                      color = col_age, 
                      base_size = font_size_odds_ratio)
    ggsave(paste0("./Plots/Odds_ratios_bonferroni_884/TFBS/EWAS/Mediated/", var, "_", dir, ".pdf"), plt, width = 5, height = 5, units = "cm")      
  }
}
