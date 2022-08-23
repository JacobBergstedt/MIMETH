library(tidyverse)
library(bacon)
library(parallel)


ewas <- collect_ewas_results("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/", 10)
fem_spec <- collect_ewas_results("./Data/Chunk_data/Results/EWAS/M_values/Female_specific/Random_effect_all_genotypes_884/Correct_for_16_props/", 10)



ewas_bacon <- bind_rows(ewas, fem_spec) %>% 
  group_by(Exposure)


groups <- unlist(group_keys(ewas_bacon))

ewas_bacon <- ewas_bacon %>% 
  group_split() %>% 
  set_names(groups) %>% 
  mclapply(function(x) bacon(effectsizes = x$Estimate, standarderrors = x$Standard_error), mc.cores = 12)

res <- ewas_bacon %>%
  map_dfr(~ tibble(Inflation = inflation(.), Bias = bias(.)), .id = "Exposure")

saveRDS(res, "./Data/RData/Results/EWAS/EWAS_bacon_inflation_884_with_fem.rds")