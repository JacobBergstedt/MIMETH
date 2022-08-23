roadmap_desc <- read_tsv("./Data/Roadmap_15_states/roadmap_descriptions.tsv")
roadmap_label_keys <- set_names(paste0(roadmap_desc$STATE_NU, "_", roadmap_desc$MNEMONIC), 
                                nm = roadmap_desc$SHORT_DESC)
meth_anno_roadmap <- readRDS("./Data/RData/Methylation/Annotation/meth_roadmap_consensus_annotation.rds") %>% 
  mutate(State_labels = fct_recode(factor(State, roadmap_label_keys), !!!roadmap_label_keys))
meth_location <- get_meth_location()
age_ewas <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect/Correct_for_16_props/Age.rds")
age_ewas <- left_join(age_ewas, meth_location) %>%
  left_join(select(meth_anno_roadmap, Probe, State, State_labels)) %>%
  left_join(meth_anno_geography) %>% 
  select(Probe, Estimate, P_value, P_FDR, Probe_gene, State, State_labels, Geography) %>%
  arrange(P_value)