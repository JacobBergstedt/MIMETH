library(tidyverse)
CMV_16_props <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/CMV_serostatus.rds") %>% 
  select(Probe, CMV_16_props = Estimate)

CMV_IDOL <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_IDOL/CMV_serostatus.rds") %>% 
  select(Probe, CMV_IDOL = Estimate)
res_EWAS <- readRDS("./Data/RData/Results/EWAS/M_values/Environment/EWAS_balance_effects_compositional_analysis_partition3_with_DC.rds")
CD4Diff <- res_EWAS %>% 
  filter(Balance == "CD4Diff_CD4Naive") %>% 
  select(Probe, CD4 = Estimate) 

CD8Diff <- res_EWAS %>% 
  filter(Balance == "CD8Diff_CD8Naive") %>% 
  select(Probe, CD8 = Estimate) 


db <- inner_join(CMV_16_props, CMV_IDOL) %>% 
  inner_join(CD4Diff) %>% 
  inner_join(CD8Diff)