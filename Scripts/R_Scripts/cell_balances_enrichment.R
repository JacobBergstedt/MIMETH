library(tidyverse)
library(missMethyl)
library(parallel)
library(broom)
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
source("./Scripts/R_scripts/Libraries/functions_for_enrichments.R")

meth_anno_roadmap <- get_anno_meth_roadmap_all_cells()
rename <- dplyr::rename

keep <- c("Mye_Lym", "T_BNK", "CD4_CD8", "CD4_Tcells", "CD4Diff_CD4Naive", "CD8Diff_CD8Naive")

res <- readRDS("./Data/RData/Results/EWAS/M_values/Environment/EWAS_balance_effects_compositional_analysis_partition3_with_DC.rds") %>% 
  filter(Balance %in% keep) %>% 
  group_by(Balance) %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  ungroup() %>% 
  left_join(meth_anno_roadmap)

res_list <- res %>% group_by(Balance)
keys <- res_list %>% group_keys()
res_list <- res_list %>% 
  group_split() %>% 
  set_names(keys$Balance)

res_chromatin_enrichments <- mclapply(res_list, run_EWAS_enrichment_pipeline_bonferroni, mc.cores = 6) %>% 
  bind_rows(.id = "Balance")

res_go <- mclapply(res_list, run_gometh, mc.cores = 6) %>% 
  bind_rows(.id = "Balance") %>%
  filter(ONTOLOGY == "BP") %>% 
  separate(col = Terms, into = c("Direction", "Thresh")) %>% 
  group_by(Balance, Direction, Thresh) %>% 
  arrange(FDR, .by_group = TRUE) %>% 
  ungroup()

saveRDS(res_chromatin_enrichments, "./Data/RData/Results/Cell_16_balances_chromatin_state_enrichments.rds")
saveRDS(res_go, "./Data/RData/Results/Cell_16_balances_go_terms.rds")


# For TSS only ------------------------------------------------------------
keep <- c("TSS", "Fl. TSS")
res_TSS <- res %>% 
  filter(Chromatin_states_Mononuclear_cells %in% keep)

keys <- factor(res_TSS$Balance)
res_list_TSS <- res_TSS %>% 
  split(f = keys)

res_go_TSS <- mclapply(res_list_TSS, gometh_enrichment_EWAS, mc.cores = 6) %>% 
  bind_rows(.id = "Balance") %>%
  filter(ONTOLOGY == "BP") %>% 
  separate(col = Terms, into = c("Direction", "Thresh")) %>% 
  group_by(Balance, Direction, Thresh) %>% 
  arrange(FDR, .by_group = TRUE) %>% 
  ungroup()


saveRDS(res_go_TSS, "./Data/RData/Results/Cell_16_balances_go_terms_TSS.rds")
