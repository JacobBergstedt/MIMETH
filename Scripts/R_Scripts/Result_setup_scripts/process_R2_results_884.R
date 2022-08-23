source("./Scripts/R_scripts/Libraries/functions_for_collections.R")
library(tidyverse)
library(broom)
library(parallel)
library(missMethyl)

res_R2 <- collect_variance_explained_results_884()

res_R2_averaged <- res_R2 %>% 
  group_by(Probe) %>% 
  summarize(across(contains("effect"), list(mean, sd)))

res_R2_mean <- res_R2_averaged %>% 
  select(Probe, contains("1"))

names(res_R2_mean) <- gsub("_1", "", names(res_R2_mean))

res_R2_sd <- res_R2_averaged %>% select(Probe, contains("2"))
names(res_R2_sd) <- gsub("_1", "", names(res_R2_sd))

res_R2_gen_int <- readRDS("./Data/RData/Results/Proportion_of_variance_explained/prop_var_explained_genetic_interactions_884.rds") %>% 
  mutate(R2 = 100 * R2) %>% 
  group_by(Probe) %>% 
  summarize(R2_mean = mean(R2))

meth_anno_roadmap <- get_anno_meth_roadmap()

vars <- c("Total_effect_intrinsic", "Total_effect_exposures", "Total_effect_cells", "Total_effect_genetic")

res_R2_mean <- res_R2_mean %>%
  left_join(select(res_R2_gen_int, Conditional_effect_genetic_interaction = R2_mean, Probe)) %>% 
  mutate(Conditional_effect_genetic_interaction = ifelse(is.na(Conditional_effect_genetic_interaction), 0, Conditional_effect_genetic_interaction)) %>% 
  rowwise() %>% 
  mutate(Top_pred = vars[which.max(c_across(all_of(vars)))]) %>% 
  ungroup() %>% 
  mutate(Top_pred = as.character(fct_recode(Top_pred, 
                                            `Cell composition` = "Total_effect_cells", 
                                            Intrinsic = "Total_effect_intrinsic", 
                                            Genetic = "Total_effect_genetic", 
                                            Exposures = "Total_effect_exposures")))

saveRDS(res_R2_mean, "./Data/RData/Results/Proportion_of_variance_explained/prop_var_explained_genetic_interactions_averaged_884.rds")



res_R2_mean %>% filter(Full_effect > 25) %>% dplyr::count(Top_pred) %>% mutate(prop = n / sum(n))
res_R2_mean %>% filter(Full_effect > 75) %>% dplyr::count(Top_pred) %>% mutate(Prop = n / sum(n))
sum(res_R2_mean$Total_effect_genetic > 25)


# Some analysis for sites explained by cell composition
target_probes <- res_R2_mean$Probe[res_R2_mean$Total_effect_cells > 50]
is_case <- meth_anno_roadmap$Probe %in% target_probes
states <- names(get_anno_roadmap_translation())
res_enrich <- map(set_names(states, states), ~ meth_anno_roadmap$State_labels == .) %>%
  map(function(is_state) xtabs(~ is_state + is_case)) %>%
  map(fisher.test) %>%
  map_dfr(tidy, .id = "States")

go <- gometh(target_probes, res_R2_mean$Probe, array.type = "EPIC", collection = "GO")
as_tibble(go) %>% arrange(FDR)