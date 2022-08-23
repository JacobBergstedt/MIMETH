library(tidyverse)
library(parallel)
library(fs)
library(patchwork)
library(scales)
library(readxl)
library(robCompositions)
source("./Scripts/R_scripts/Libraries/functions_for_plotting.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
source("./Scripts/R_scripts/Libraries/functions_for_compositions.R")
select <- dplyr::select

# balances <- readRDS("./Data/RData/Results/EWAS/M_values/Environment/EWAS_balance_effects_compositional_analysis_partition2.rds") %>% 
#   group_by(Balance) %>% 
#   mutate(P_bonferroni = p.adjust(P_value)) %>% 
#   ungroup()

# balances_adjusted <- readRDS("./Data/RData/Results/EWAS/M_values/Environment/EWAS_balance_effects_compositional_analysis_partition2_adjusted.rds")

# CD8CMEMEMRA_CD8Naive <- balances %>% filter(Balance == "CD8CMEMEMRA_CD8Naive")
# CD8EMEMRA_CD8CM <- balances %>% filter(Balance == "CD8EMEMRA_CD8CM")
# CD4EMEMRA_CD4CM <- balances %>% filter(Balance == "CD4EMEMRA_CD4CM")
# CD4CMEMEMRA_CD4Naive <- balances %>% filter(Balance == "CD4CMEMEMRA_CD4Naive")

age_IDOL <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_IDOL/Age.rds") %>% 
  add_column(Adjustment = "IDOL")

# age_panel5 <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_6_props/Age.rds") %>% 
#   add_column(Adjustment = "Panel_5")

age_16_cells <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Age.rds") %>% 
  add_column(Adjustment = "16_cells")

age <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_no_cells/Age.rds") %>% 
  add_column(Adjustment = "None")

CMV_IDOL <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_IDOL/CMV_serostatus.rds") %>% 
  add_column(Adjustment = "IDOL")

CMV_16_cells <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/CMV_serostatus.rds") %>% 
  add_column(Adjustment = "16_cells")

# CMV_panel5 <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_6_props/CMV_serostatus.rds") %>% 
#   add_column(Adjustment = "Panel_5")

nr_sign_table <- readRDS("./Tables/EWAS__884_NR_sign_small_table.rds")
res_PCA <- readRDS("./Data/RData/Results/meth_PCA_regression_results.rds") %>% 
  mutate(PC = factor(PC, unique(PC)))
meth_PC <- readRDS("./Data/RData/Methylation/Meth_PCA.rds")
covariates <- get_covs_884()


cell_list_16_cells <- get_cell_list()
cell_list <- cell_list_16_cells[!cell_list_16_cells == "X_dendritic_cells.panel8"]

cells <- covariates %>% 
  select(SUBJID, all_of(cell_list))

cell_balances <- add_column(balance_preds(cells[-1], sbp = read_xlsx("./Data/Cell_subset_partitions_2.xlsx")), 
                            SUBJID = cells[["SUBJID"]], .before = 1)


# compute_correlation <- function(x, y) {
#   
#   cor(x$Estimate, y$Estimate)
#   
# }
# 
# f <- factor(balances$Balance)
# IDOL <- balances %>% 
#   split(f = f) %>%
#   map_dbl( ~ compute_correlation(x = ., y = CMV_IDOL))
# 
# 
# f <- factor(balances$Balance)
# Cells_16 <- balances %>% 
#   split(f = f) %>%
#   map_dbl( ~ compute_correlation(x = ., y = CMV_16_cells))
# 
# 
# f <- factor(balances$Balance)
# IDOL <- balances %>% 
#   split(f = f) %>%
#   map_dbl( ~ compute_correlation(x = ., y = age_IDOL))
# 
# 
# f <- factor(balances$Balance)
# Cells_16 <- balances %>% 
#   split(f = f) %>%
#   map_dbl( ~ compute_correlation(x = ., y = age_16_cells))



# Heatmap -----------------------------------------------------------------


plt <- nr_sign_table %>%
  ggplot(aes(x = Adjustment_labels, y = fct_rev(Exposure_labels), fill = NR_sign)) +
  geom_raster() +
  geom_text(aes(label = format(NR_sign, big.mark = ",")), color = "white", size = 4) +
  scale_fill_gradient(low = blues[2], high = blues[9], trans = "sqrt") +
  scale_x_discrete(position = "top") + 
  ylab(label = NULL) +
  xlab(label = NULL) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        panel.grid = element_blank())
  
ggsave("./Plots/NR_sign_heatmap_sqrtscale.pdf", plt, height = 20, width = 12, units = "cm")



plt <- nr_sign_table %>%
  ggplot(aes(x = Adjustment_labels, y = fct_rev(Exposure_labels), fill = NR_sign)) +
  geom_raster() +
  geom_text(aes(label = format(NR_sign, big.mark = ",")), color = "white", size = 4) +
  scale_fill_gradient(low = blues[4], high = blues[9]) +
  ylab(label = NULL) +
  xlab(label = NULL) + 
  theme_classic() +
  theme(legend.position = "none")

ggsave("./Plots/NR_sign_heatmap.pdf", plt, height = 20, width = 12, units = "cm")

pltA <- plt



# IDOL v DIRECT -----------------------------------------------------------

# Age ---------------------------------------------------------------------

plt_frame <- age_16_cells %>% 
  select(Probe, Estimate_16_cells = Estimate, P_16_cells = P_bonferroni) %>% 
  left_join(select(age_IDOL, Probe, Estimate_IDOL = Estimate, P_IDOL = P_bonferroni))


plt <- plt_frame %>% 
  filter(P_16_cells < 0.05 | P_IDOL < 0.05) %>% 
  ggplot(aes(x = Estimate_IDOL, y = Estimate_16_cells)) +
  geom_point(color = col_age) +
  xlim(-1, 1) +
  ylim(-1, 1) + 
  geom_abline(size = 1, color = "orange") +
  theme_classic()

ggsave("./Plots/Revision_plots/scatter_Age_IDOL_v_16_cells.pdf", plt, width = 30, height = 20, units = "cm")

pltB <- plt


# CMV ---------------------------------------------------------------------

plt_frame <- CMV_16_cells %>% 
  select(Probe, Estimate_16_cells = Estimate, P_16_cells = P_bonferroni) %>% 
  left_join(select(CMV_IDOL, Probe, Estimate_IDOL = Estimate, P_IDOL = P_bonferroni)) %>%
  left_join(select(CD8CMEMEMRA_CD8Naive, Probe, Estimate_CD8CMEMEMRA = Estimate, P_CD8CMEMEMRA_Naive = P_bonferroni)) %>% 
  left_join(select(CD8EMEMRA_CD8CM, Probe, Estimate_CD8CM = Estimate, P_CD8EMEMRA_CM = P_bonferroni)) %>% 
  left_join(select(CD4CMEMEMRA_CD4Naive, Probe, Estimate_CD4CMEMEMRA = Estimate, P_CD4CMEMEMRA_Naive = P_bonferroni)) %>% 
  left_join(select(CD4EMEMRA_CD4CM, Probe, Estimate_CD4EMRA = Estimate, P_CD4EMEMRA_CM = P_bonferroni))


plt <- plt_frame %>% 
  filter(P_16_cells < 0.05 | P_IDOL < 0.05) %>% 
  ggplot(aes(x = Estimate_IDOL, y = Estimate_16_cells)) +
  geom_point(color = col_CMV) +
  xlim(-1, 1) +
  ylim(-1, 1) + 
  geom_abline(size = 1, color = "orange") +
  theme_classic()

ggsave("./Plots/Revision_plots/scatter_CMV_IDOL_v_16_cells.pdf", plt, width = 30, height = 20, units = "cm")

pltC <- plt


plt1 <- plt_frame %>% 
  filter(P_IDOL < 0.05) %>% 
  ggplot(aes(x = Estimate_IDOL, y = Estimate_CD4CMEMEMRA)) +
  geom_point(color = col_CMV) +
  xlab("IDOL estimate") +
  ylab("CD4 Naive v CM * EM * EMRA") +
  xlim(-1, 1) +
  ylim(-1, 1) +
  theme_bw()

ggsave("./Plots/Revision_plots/scatter_CMV_IDOL_v_CD4Diff_balance.pdf", plt1, width = 30, height = 20, units = "cm")

plt2 <- plt_frame %>% 
  filter(P_IDOL < 0.05) %>% 
  ggplot(aes(x = Estimate_16_cells, y = Estimate_CD4CMEMEMRA)) +
  geom_point(color = col_CMV) +
  xlab("16 cells estimate") +
  ylab(NULL) +
  xlim(-1, 1) +
  ylim(-1, 1) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank())

ggsave("./Plots/Revision_plots/scatter_CMV_16_cells_v_CD4Diff_balance.pdf", plt2, width = 30, height = 20, units = "cm")

plt <- plt1 + plt2
ggsave("./Plots/Revision_plots/scatter_CMV_est_v_CD4Diff_balance.pdf", plt, width = 30, height = 20, units = "cm")

plt <- plt_frame %>% 
  filter(P_IDOL < 0.05) %>% 
  ggplot(aes(x = Estimate_IDOL, y = Estimate_CD8CM)) +
  geom_point() +
  theme_bw() +
  xlab("16 cells estimate") +
  ylab("CD4 Naive v CM * EM * EMRA") +
  xlim(-1, 1) +
  ylim(-0.5, 0.5)

# Cell type estimates ---------------------------------------------------------
cells_panel5 <- get_panel5_cells()
norm_const <- rowSums(cells_panel5[-1])
cells_panel5 <- mutate_at(cells_panel5, .vars = vars(!SUBJID), ~ . / norm_const) %>% 
  select(-X_CD4negCD8neg_of_total.panel5)
keys <- names(cells_panel5)[-1]
names(keys) <- c("IDOL_NK", "IDOL_Mono", "IDOL_CD4T", "IDOL_CD8T", "IDOL_Bcell", "IDOL_Neu")

plt_frame_panel5 <- pivot_longer(cells_panel5, -SUBJID, names_to = "Panel5_name", values_to = "Panel5_value") %>% 
  mutate(Cell = fct_recode(factor(Panel5_name, levels = keys), !!!keys))

plt_frame_IDOL <- pivot_longer(covariates[c("SUBJID", names(keys))], -SUBJID, names_to = "IDOL_name", values_to = "IDOL_value") %>% 
  mutate(Cell = IDOL_name)

plt_frame <- inner_join(plt_frame_IDOL, plt_frame_panel5)

plt <-  plt_frame %>% 
  ggplot(aes(x = IDOL_value, y = Panel5_value)) +
  geom_point() +
  geom_abline(size = 1) +
  facet_wrap(vars(Cell)) + 
  theme_bw()

pltD <- plt

ggsave("./Plots/Revision_plots/cell_estimates_comparison.pdf", plt, height = 10, width = 16, units = "cm")

# PC plots ----------------------------------------------------------------

plt_frame <- meth_PC[c("SUBJID", "PC1", "PC2")] %>% 
  inner_join(select(covariates, SUBJID, CMV_serostatus, Age, all_of(cell_list))) %>% 
  inner_join(cell_balances)

plt <- plt_frame %>% 
  ggplot(aes(x = PC1, y = PC2, color = CMV_serostatus)) +
  geom_point(size = 2) +
  scale_color_manual(values = c(reds[3], col_CMV)) +
  theme_bw()

plt

ggsave("./Plots/Revision_plots/PC1_PC2_CMV.pdf", plt, height = 20, width = 30, units = "cm")

plt <- plt_frame %>% 
  ggplot() +
  geom_point(aes(x = CD4CMEMEMRA_CD4Naive, y = PC1, color = CMV_serostatus)) +
  scale_color_manual(values = c(reds[3], col_CMV)) +
  xlab("f(CD4 CM,EM,EMRA / CD4 Naive)") +
  ylab("Meth PC1") +
  geom_smooth(aes(x = CD4CMEMEMRA_CD4Naive, y = PC1), method = "lm", se = FALSE, color = "black") +
  guides(color = guide_legend(title = "CMV")) + 
  theme_bw() +
  theme(legend.position = c(0.85, 0.85))

ggsave("./Plots/Revision_plots/PC1_CD4_CMV.pdf", plt, width = 30, height = 20, units = "cm")

pltE <- plt


# Combine plot ------------------------------------------------------------


design <- "
  AABC
  AADE
"

plt <- pltA + pltB + pltC + pltD + pltE + plot_layout(design = design) &
  theme(plot.tag = element_text(size = 6, face = "bold"))

ggsave("./Plots/Revision_plots/Figure_1.png", plt, width = 172, height = 147, units = "mm", dpi = 500)
ggsave("./Plots/Revision_plots/Figure_1.pdf", plt, width = 172, height = 147, units = "mm")
