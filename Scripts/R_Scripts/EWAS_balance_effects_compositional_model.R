## Packages

library(MASS)
library(lme4)
library(tidyverse)
library(parallel)
library(splines)
library(pbkrtest)
library(sandwich)
library(broom)
library(readxl)
source("./Scripts/R_scripts/Libraries/functions_for_mediation2.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_compositions.R")
library(robCompositions)

fit_compositional_model <- function(probe, snps, fm, cell_balances, covariates = NULL) {
  
  db <- cbind(y = probe, cell_balances)
  
  if (!is.null(covariates)) {
    
    fm <- update(fm, . ~ ns(Age, df = 3) + Sex + CMV_serostatus + Smoking_status + Ancestry_PC1 + Ancestry_PC2 + .)
    db <- cbind(db, covariates)
    
  }
  
  if (!is.null(snps)) {
    
    fm <- fm_add_snps(fm, snps)
    db <- cbind(db, snps)
    
  }
  
  
  m <- lm(fm, db)
  tidy(m) %>% 
    filter(term %in% names(cell_balances)) %>% 
    select(Balance = term, Estimate = estimate, Standard_error = std.error, P_value = p.value)
  
}

meth_anno_loc <- get_anno_meth_location() 


# Run analysis ------------------------------------------------------------

n_cores <- 12

## Data
cell_list_16_cells <- get_cell_list()

meth <- get_m_values()
snp_mat_list <- get_snp_mat_per_probe()

covariates <- get_covs_884() %>% 
  mutate(Age = Age / 50) %>% 
  select(SUBJID, Sex, Smoking_status, CMV_serostatus, Age, Ancestry_PC1, Ancestry_PC2, all_of(cell_list_16_cells)) %>% 
  filter(!if_any(all_of(cell_list_16_cells), .fns = ~ . <= 0))

snp_mat_list <- mclapply(snp_mat_list, function(x) x[covariates$SUBJID, , drop = FALSE], mc.cores = n_cores)
meth <- meth[match(covariates$SUBJID, meth$SUBJID),]
meth$SUBJID <- NULL


cells <- covariates %>% 
  select(all_of(cell_list_16_cells))

cell_balances <- balance_preds(cells, sbp = read_xlsx("./Data/Cell_subset_partitions_3_with_DC.xlsx"))
covariates <- covariates %>% select(Age, Sex, CMV_serostatus, Smoking_status, Ancestry_PC1, Ancestry_PC2)

fm <- as.formula(paste0("y ~ ", paste0(names(cell_balances), collapse = " + ")))
fit_compositional_model(meth[[2]], snp_mat_list[[2]], fm = fm, cell_balances = cell_balances, covariates = covariates)
fit_compositional_model(meth[[2]], snp_mat_list[[2]], fm = fm, cell_balances = cell_balances, covariates = NULL)

res <- mcmapply(fit_compositional_model, meth, snp_mat_list, MoreArgs = list(fm = fm, cell_balances = cell_balances, covariates = NULL),  
                SIMPLIFY = FALSE, mc.cores = 12)
res <- res %>% 
  bind_rows(.id = "Probe") %>% 
  left_join(meth_anno_loc) %>% 
  group_by(Balance) %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  ungroup()

saveRDS(res, "./Data/RData/Results/EWAS/M_values/Environment/EWAS_balance_effects_compositional_analysis_partition3_with_DC.rds")


# 
# p %>%
#   filter(Cell %in% c("X_VIABLE_NEUTROPHILS_OF_TOTAL.panel7", "X_CD4_naive_of_total.panel1")) %>%
#   filter(p.adjust(P_value) < 0.05) %>% 
#   select(Cell, Estimate, Probe) %>%
#   pivot_wider(names_from = Cell, values_from = Estimate) %>%
#   ggplot(aes(x = X_VIABLE_NEUTROPHILS_OF_TOTAL.panel7, y = X_CD4_naive_of_total.panel1)) +
#   geom_point()
# 
# 
# 
# p <- o %>% 
#   rename(Cell = term, Estimate = estimate, Standard_error = std.error, P_value = p.value)
# 
# nr_sign <- p %>%
#   group_by(Cell) %>%
#   mutate(P_bonferroni = p.adjust(P_value)) %>%
#   summarize(Nr_sign = sum(P_bonferroni < 0.05))
# 
# cells <- get_covs_884()[get_cell_list()]
# cell_props <- colMeans(cells)
# cell_props <- tibble(Cell = names(cell_props), Props = 100 * cell_props)
# 
# keys <- tibble(
#   Cell = c("X_VIABLE_NEUTROPHILS_OF_TOTAL.panel7",
#            "X_VIABLE_BASOPHILS_OF_TOTAL.panel7",
#            "X_VIABLE_EOSINOPHILS_OF_TOTAL.panel7",
#            "X_mono_of_total.panel5",
#            "X_CD4_naive_of_total.panel1",
#            "X_CD4_CM_of_total.panel1",
#            "X_CD4_EM_of_total.panel1",
#            "X_CD4_EMRA_of_total.panel1",
#            "X_CD8_naive_of_total.panel1",
#            "X_CD8_CM_of_total.panel1",
#            "X_CD8_EM_of_total.panel1",
#            "X_CD8_EMRA_of_total.panel1",
#            "X_CD8bnegCD4neg_of_total.panel1",
#            "X_Bcells_ofTotal.panel6",
#            "X_NK_of_total.panel4",
#            "X_dendritic_cells.panel8"),
#   Name = c("Neu",
#            "Basophils",
#            "Eosinophils",
#            "Monocytes",
#            "CD4 naive",
#            "CD4 CM",
#            "CD4 EM",
#            "CD4 EMRA",
#            "CD8 naive",
#            "CD8 CM",
#            "CD8 EM",
#            "CD8 EMRA",
#            "CD8bnegCD4neg",
#            "B_cells",
#            "NK",
#            "DC")
# )
# 
# plt_frame <- inner_join(nr_sign, cell_props) %>%
#   inner_join(keys)
# 
# plt <- ggplot(plt_frame, aes(y = Props, x = Nr_sign, label = Name)) +
#   geom_point() +
#   geom_text(hjust = 0, nudge_x = 300, nudge_y = 0.05) +
#   ylab("Cell proportion") +
#   xlab("Nr sign. sites") +
#   theme_bw()
# 
# plt


# 
# ggsave("./Plots/EWAS_cell_types.pdf", plt, width = 20, height = 15, units = "cm")