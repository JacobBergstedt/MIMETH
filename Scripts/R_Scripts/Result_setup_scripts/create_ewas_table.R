
library(tidyverse)
library(parallel)
n_cores <- 10

# Load ewas tibbles -------------------------------------------------------
ewas <- collect_ewas_results("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/", n_cores)
ewas_total <- collect_ewas_results("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_no_cells/", n_cores)
fem_spec <- collect_ewas_results("./Data/Chunk_data/Results/EWAS/M_values/Female_specific/Random_effect_all_genotypes_884/Correct_for_16_props/", n_cores)
fem_spec_total <- collect_ewas_results("./Data/Chunk_data/Results/EWAS/M_values/Female_specific/Random_effect_all_genotypes_884/Correct_for_no_cells/", n_cores)
ewas_effect_decomposition <- collect_ewas_results("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/", n_cores)

mediation <- ewas_effect_decomposition %>%
  filter(Effect == "Mediation") %>%
  group_by(Variable) %>%
  mutate(P_bonferroni = p.adjust(P_value),
         P_FDR = p.adjust(P_value, "BH")) %>%
  ungroup()

ewas_IDOL <- collect_ewas_results("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_IDOL/", n_cores)
ewas_extended_IDOL <- collect_ewas_results("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_extended_IDOL/", n_cores)
ewas_6_cells <- collect_ewas_results("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_6_props/", n_cores)
ewas_Houseman <- collect_ewas_results("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_Houseman/", n_cores)

# Assign significance -----------------------------------------------------
nr_assoc_ewas <- ewas %>% 
  distinct(Exposure, Probe, .keep_all = TRUE) %>% 
  group_by(Exposure) %>% 
  summarize(Direct_Bonferroni = sum(P_bonferroni < 0.05)) %>% 
  arrange(desc(Direct_Bonferroni))

nr_assoc_ewas_total <- ewas_total %>%
  distinct(Exposure, Probe, .keep_all = TRUE) %>% 
  group_by(Exposure) %>% 
  summarize(Total_Bonferroni = sum(P_bonferroni < 0.05)) %>% 
  arrange(desc(Total_Bonferroni))

nr_assoc_fem <- fem_spec %>%
  distinct(Exposure, Probe, .keep_all = TRUE) %>%
  group_by(Exposure) %>%
  summarize(Direct_Bonferroni = sum(P_FDR < 0.05)) %>%
  arrange(desc(Direct_Bonferroni))

nr_assoc_fem_total <- fem_spec_total %>%
  distinct(Exposure, Probe, .keep_all = TRUE) %>%
  group_by(Exposure) %>%
  summarize(Total_Bonferroni = sum(P_bonferroni < 0.05)) %>%
  arrange(desc(Total_Bonferroni))

nr_assoc_med <- mediation %>% 
  filter(Effect == "Mediation") %>% 
  group_by(Variable) %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  summarize(Mediated_Bonferroni = sum(P_bonferroni < 0.05)) %>% 
  rename(Exposure = Variable) %>% 
  mutate(Exposure = recode(Exposure, CMV = "CMV_serostatus", Smoking = "Smoking_status"))

nr_assoc_IDOL <- ewas_IDOL %>% 
  distinct(Exposure, Probe, .keep_all = TRUE) %>%
  group_by(Exposure) %>%
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  summarize(IDOL_Bonferroni = sum(P_bonferroni < 0.05))

nr_assoc_6_cells <- ewas_6_cells %>% 
  distinct(Exposure, Probe, .keep_all = TRUE) %>%
  group_by(Exposure) %>%
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  summarize(Cells_6_Bonferroni = sum(P_bonferroni < 0.05))


nr_assoc_extended_IDOL <- ewas_extended_IDOL %>% 
  distinct(Exposure, Probe, .keep_all = TRUE) %>%
  group_by(Exposure) %>%
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  summarize(IDOL12_Bonferroni = sum(P_bonferroni < 0.05))

nr_assoc_Houseman <- ewas_Houseman %>%
  distinct(Exposure, Probe, .keep_all = TRUE) %>%
  group_by(Exposure) %>%
  mutate(P_bonferroni = p.adjust(P_value)) %>%
  summarize( Houseman_Bonferroni = sum(P_bonferroni < 0.05))

# Create table ------------------------------------------------------------

ewas_table <- bind_rows(nr_assoc_ewas_total, nr_assoc_fem_total) %>% 
  left_join(bind_rows(nr_assoc_ewas, nr_assoc_fem)) %>%
  left_join(nr_assoc_extended_IDOL) %>% 
  left_join(nr_assoc_6_cells) %>%
  left_join(nr_assoc_IDOL) %>%
  left_join(nr_assoc_Houseman) %>%
  left_join(nr_assoc_med) %>%
  filter(!Exposure %in% c("Year_of_last_pregnancy",
                          "Log WBC hematology",
                          "MCHC", 
                          "Log_lymphocytes_hematology", 
                          "Log1p_eosinophils_hematology", 
                          "Log_mean_corpuscular_volume", 
                          "Log_monocytes_hematology", 
                          "Log_neutrophils_hematology", 
                          "Log_RBC",
                          "Log1p_basophils_hematology")) %>%
  arrange(desc(Total_Bonferroni), desc(Direct_Bonferroni), desc(Cells_6_Bonferroni), desc(IDOL12_Bonferroni), desc(IDOL_Bonferroni), desc(Mediated_Bonferroni)) %>% 
  map_dfc(as.character)

names(ewas_table) <- c("Exposure", "None", "16 cells", "IDOL 12 cells", "6 cells", "IDOL 6 cells", "Houseman 6 cells", "Mediation")
ewas_table[is.na(ewas_table)] <- ""
ewas_table <- ewas_table %>% 
  mutate(Exposure = str_replace_all(Exposure, "_", " "), across(!Exposure, format, big.mark = ","))



saveRDS(ewas_table, "./Tables/ewas_table_884_Bonferroni.rds")
write_tsv(ewas_table, "./Tables/ewas_table_884_Bonferroni.tsv")

ewas_table_small <- ewas_table %>% 
  filter(Exposure %in% c("Age", "CMV serostatus", "Sex", "Heart rate", "Log WBC hematology", "Smoking status", "Temperature ear", "CMV serology", "Log CRP levels", "Hour of sampling"))



ewas_table_small_long <- pivot_longer(ewas_table_small, cols = -Exposure, names_to = "Adjustment", values_to = "NR_sign") %>% mutate(NR_sign = as.integer(NR_sign))  
order_tib <- tibble(Exposure = c("Age", "CMV serostatus", "Sex", "Smoking status", "Log CRP levels", "Heart rate", "Temperature ear", "Hour of sampling"))
ewas_table_small_long <- left_join(order_tib, ewas_table_small_long)
saveRDS(ewas_table_small_long, "./Tables/EWAS_884_NR_sign_small_table.rds")  



# Balances ----------------------------------------------------------------

balances <- readRDS("./Data/RData/Results/EWAS/M_values/Environment/EWAS_balance_effects_compositional_analysis_partition3.rds") %>%
  group_by(Balance) %>%
  mutate(P_bonferroni = p.adjust(P_value)) %>%
  ungroup() %>% 
  left_join(get_anno_meth_location()) %>% 
  left_join(get_anno_meth_roadmap())


keys <- unique(balances$Balance)
names(keys) <- c("Myeloid v Lymphoid", "Neutrophils vother myeloid cells", "Monocytes v non neutrophil myeloid cells",
                 "Basophils v eosinophils", "T cell v B and NK cells", "B cell v NK cells", "CD8negCD4neg v T cells", 
                 "CD4 v CD8", "CD4 diff v CD4 naive", "CD4 EM and EMRA v CD4 CM", "CD4 EMRA v CD4EM", "CD8 diff v CD8 naive",
                 "CD8 EM and EMRA v CD8 CM", "CD8 EMRA v CD8 EM")


balances_summary <- balances %>% 
  group_by(Balance) %>% 
  summarize(NR_sign = sum(P_bonferroni < 0.05))

balances_summary <- left_join(tibble(Balance_labels = names(keys), Balance = keys), balances_summary)
saveRDS(balances_summary, "./Tables/Balances_NR_sign_table.rds")




