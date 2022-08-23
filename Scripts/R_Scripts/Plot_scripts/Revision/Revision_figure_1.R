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

balances <- readRDS("./Data/RData/Results/EWAS/M_values/Environment/EWAS_balance_effects_compositional_analysis_partition3.rds") %>%
  group_by(Balance) %>%
  mutate(P_bonferroni = p.adjust(P_value)) %>%
  ungroup() %>% 
  left_join(get_anno_meth_location()) %>% 
  left_join(get_anno_meth_roadmap())

# balances_adjusted <- readRDS("./Data/RData/Results/EWAS/M_values/Environment/EWAS_balance_effects_compositional_analysis_partition2_adjusted.rds")

CD8Diff_CD8Naive <- balances %>% filter(Balance == "CD8Diff_CD8Naive")
CD4Diff_CD4Naive <- balances %>% filter(Balance == "CD4Diff_CD4Naive")

age_IDOL <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_IDOL/Age.rds") %>% 
  add_column(Adjustment = "IDOL")

age_16_cells <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Age.rds") %>% 
  add_column(Adjustment = "16_cells")

age <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_no_cells/Age.rds") %>% 
  add_column(Adjustment = "None")

age_6_cells <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_6_props/Age.rds") %>% 
  add_column(Adjustment = "6_cells")

CMV_IDOL <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_IDOL/CMV_serostatus.rds") %>% 
  add_column(Adjustment = "IDOL")

CMV_16_cells <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/CMV_serostatus.rds") %>% 
  add_column(Adjustment = "16_cells")

CMV_total <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_no_cells/CMV_serostatus.rds") %>% 
  add_column(Adjustment = "None")

CMV_6_cells <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_6_props/CMV_serostatus.rds") %>% 
  add_column(Adjustment = "6_cells")

# CMV_panel5 <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_6_props/CMV_serostatus.rds") %>% 
#   add_column(Adjustment = "Panel_5")

nr_sign_table <- readRDS("./Tables/EWAS_884_NR_sign_small_table.rds")
res_PCA <- readRDS("./Data/RData/Results/meth_PCA_regression_results_3.rds") %>% 
  mutate(PC = factor(PC, unique(PC)))
meth_PC <- readRDS("./Data/RData/Methylation/Meth_PCA.rds")
covariates <- get_covs_884()


cell_list_16_cells <- get_cell_list()
cell_list <- cell_list_16_cells[!cell_list_16_cells == "X_dendritic_cells.panel8"]

cells <- covariates %>% 
  select(SUBJID, all_of(cell_list))

cell_balances <- add_column(balance_preds(cells[-1], sbp = read_xlsx("./Data/Cell_subset_partitions_3.xlsx")), 
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

nr_sign_table <- readRDS("./Tables/EWAS_884_NR_sign_small_table.rds") %>% 
  filter(Adjustment != "Mediation")
keys <- c("16 cells", "IDOL 12 cells", "6 cells", "IDOL 6 cells", "Houseman 6 cells", "None")
names(keys) <- c("16 cells", "IDOL 12", "6 cells", "IDOL 6", "Houseman 6", "None")
nr_sign_table <- nr_sign_table %>%
  mutate(Adjustment = fct_recode(factor(Adjustment, keys), !!!keys))

plt <- nr_sign_table %>%
  ggplot(aes(x = Adjustment, y = fct_rev(factor(Exposure, unique(Exposure))), fill = NR_sign)) +
  geom_raster() +
  geom_text(aes(label = format(NR_sign, big.mark = ",")), color = "white", size = font_size_label) +
  scale_fill_gradient(low = blues[4], high = blues[9], trans = "sqrt") +
  scale_x_discrete(position = "top") + 
  ylab(label = NULL) +
  xlab(label = NULL) + 
  theme_bw(base_size = font_size_scatter) +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  labs(tag = "A")
  
ggsave("./Plots/Revision_plots/Figure1/NR_sign_heatmap_sqrtscale.pdf", plt, width = 10, height = 10, units = "cm")



# plt <- nr_sign_table %>%
#   ggplot(aes(x = Adjustment_labels, y = fct_rev(Exposure_labels), fill = NR_sign)) +
#   geom_raster() +
#   geom_text(aes(label = format(NR_sign, big.mark = ",")), color = "white", size = 4) +
#   scale_fill_gradient(low = blues[4], high = blues[9]) +
#   ylab(label = NULL) +
#   xlab(label = NULL) + 
#   theme_classic() +
#   theme(legend.position = "none")
# 
# ggsave("./Plots/NR_sign_heatmap.pdf", plt, height = 20, width = 12, units = "cm")

pltA <- plt



# IDOL v DIRECT -----------------------------------------------------------

# Age ---------------------------------------------------------------------

plt_frame_age <- age_16_cells %>% 
  select(Probe, Estimate_16_cells = Estimate, P_16_cells = P_bonferroni) %>% 
  left_join(select(age, Probe, Estimate_total = Estimate, P_total = P_bonferroni)) %>%
  left_join(select(age_6_cells, Probe, Estimate_6_cells = Estimate, P_6_cells = P_bonferroni)) %>%
  left_join(select(age_IDOL, Probe, Estimate_IDOL = Estimate, P_IDOL = P_bonferroni)) %>%
  left_join(select(CD8Diff_CD8Naive, Probe, Estimate_CD8Diff_CD8Naive = Estimate, P_CD8Diff_CD8Naive = P_bonferroni)) %>% 
  left_join(select(CD4Diff_CD4Naive, Probe, Estimate_CD4Diff_CD4Naive = Estimate, P_CD4Diff_CD4Naive = P_bonferroni))

plt_frame <- age_16_cells %>% 
  select(Probe, Estimate_16_cells = Estimate, P_16_cells = P_bonferroni) %>% 
  left_join(select(age_IDOL, Probe, Estimate_IDOL = Estimate, P_IDOL = P_bonferroni))


plt <- plt_frame %>% 
  filter(P_16_cells < 0.05 | P_IDOL < 0.05) %>% 
  ggplot(aes(x = Estimate_IDOL, y = Estimate_16_cells)) +
  geom_point(size = point_size_scatter, color = col_age) +
  geom_abline(size = line_size_blacks, color = "black") +
  theme_bw(base_size = font_size_scatter) +
  labs(tag = "B")

ggsave("./Plots/Revision_plots/Figure1/scatter_Age_IDOL_v_16_cells.pdf", plt, width = 30, height = 20, units = "cm")

pltB <- plt


# CMV ---------------------------------------------------------------------

plt_frame <- CMV_16_cells %>% 
  select(Probe, Estimate_16_cells = Estimate, P_16_cells = P_bonferroni) %>% 
  left_join(select(CMV_total, Probe, Estimate_total = Estimate, P_total = P_bonferroni)) %>%
  left_join(select(CMV_6_cells, Probe, Estimate_6_cells = Estimate, P_6_cells = P_bonferroni)) %>%
  left_join(select(CMV_IDOL, Probe, Estimate_IDOL = Estimate, P_IDOL = P_bonferroni)) %>%
  left_join(select(CD8Diff_CD8Naive, Probe, Estimate_CD8Diff_CD8Naive = Estimate, P_CD8Diff_CD8Naive = P_bonferroni)) %>% 
  left_join(select(CD4Diff_CD4Naive, Probe, Estimate_CD4Diff_CD4Naive = Estimate, P_CD4Diff_CD4Naive = P_bonferroni))


plt <- plt_frame %>% 
  filter(P_16_cells < 0.05 | P_IDOL < 0.05) %>% 
  ggplot(aes(x = Estimate_IDOL, y = Estimate_16_cells)) +
  geom_point(size = point_size_scatter, color = col_CMV) +
  geom_abline(size = line_size_blacks, color = "black") +
  theme_bw(base_size = font_size_scatter) +
  labs(tag = "C")

ggsave("./Plots/Revision_plots/Figure1/scatter_CMV_IDOL_v_16_cells.pdf", plt, width = 10, height = 10, units = "cm")

pltC <- plt


pltF <- plt_frame %>% 
  filter(P_IDOL < 0.05) %>% 
  ggplot(aes(x = Estimate_IDOL, y = Estimate_CD4Diff_CD4Naive)) +
  geom_point(size = point_size_scatter, color = col_CMV) +
  xlim(-1, 1) +
  ylim(-1, 1) +
  xlab("IDOL estimate") +
  ylab("CD4 Diff. effect") +
  theme_bw(base_size = font_size_scatter) +
  labs(tag = "F")

ggsave("./Plots/Revision_plots/Figure1/scatter_CMV_IDOL_v_CD4Diff_balance.pdf", pltF, width = 10, height = 10, units = "cm")

pltG <- plt_frame %>% 
  filter(P_IDOL < 0.05) %>% 
  ggplot(aes(x = Estimate_16_cells, y = Estimate_CD4Diff_CD4Naive)) +
  geom_point(size = point_size_scatter, color = col_CMV) +
  xlim(-1, 1) +
  ylim(-1, 1) +
  xlab("16 cells estimate") +
  ylab("CD4 Diff. effect") +
  theme_bw(base_size = font_size_scatter) +
  labs(tag = "G")

ggsave("./Plots/Revision_plots/Figure1/scatter_CMV_16_cells_v_CD4Diff_balance.pdf", pltG, width = 10, height = 10, units = "cm")

pltH <- plt_frame %>% 
  filter(P_total < 0.05) %>% 
  ggplot(aes(x = Estimate_total, y = Estimate_CD4Diff_CD4Naive)) +
  geom_point(size = point_size_scatter, color = col_CMV) +
  xlim(-1, 1) +
  ylim(-1, 1) +
  xlab("Total est.") +
  ylab("CD4 Diff. effect") +
  theme_bw(base_size = font_size_scatter) +
  labs(tag = "H")

ggsave("./Plots/Revision_plots/Figure1/scatter_CMV_total_v_CD4Diff_balance.pdf", pltH, width = 10, height = 10, units = "cm")



plt <- pltF | pltG | pltH
ggsave("./Plots/Revision_plots/Figure1/scatter_CMV_v_CD4Diff_balance.pdf", plt, width = 16, height = 16, units = "cm")

# PC plots ----------------------------------------------------------------

plt_frame <- meth_PC[c("SUBJID", "PC1", "PC2", "PC3")] %>% 
  inner_join(select(covariates, SUBJID, CMV_serostatus, Age, all_of(cell_list))) %>% 
  inner_join(cell_balances) %>% 
  mutate(CD8pheno = ifelse(CD8Diff_CD8Naive > median(CD8Diff_CD8Naive), "Diff.", "Naive"))

plt <- plt_frame %>% 
  ggplot(aes(x = PC1, y = PC2, color = CMV_serostatus)) +
  geom_point(size = point_size_scatter) +
  scale_color_manual(values = c(reds[3], col_CMV)) +
  theme_bw(base_size = font_size_scatter)

plt

ggsave("./Plots/Revision_plots/Figure1/PC1_PC2_CMV.pdf", plt, height = 20, width = 30, units = "cm")

plt <- plt_frame %>% 
  ggplot() +
  geom_point(aes(x = CD4Diff_CD4Naive, y = PC1, color = CMV_serostatus), size = point_size_scatter) +
  scale_color_manual(values = c(reds[3], col_CMV)) +
  xlab("CD4 Diff.") +
  ylab("Meth PC1") +
  geom_smooth(aes(x = CD4Diff_CD4Naive, y = PC1), method = "lm", se = FALSE, size = line_size_blacks, color = "black") +
  guides(color = guide_legend(title = "CMV")) + 
  theme_bw(base_size = font_size_scatter) +
  theme(legend.position = c(0.8, 0.9),
        legend.key.size = unit(0.1,"cm")) +
  labs(tag = "D")

ggsave("./Plots/Revision_plots/Figure1/PC1_CD4_CMV.pdf", plt, width = 10, height = 10, units = "cm")

pltD <- plt

plt <- plt_frame %>% 
  ggplot() +
  geom_point(aes(x = Age, y = PC3, color = fct_rev(CD8pheno)), size = point_size_scatter) +
  scale_color_manual(values = c(oranges[3], col_cell)) +
  xlab("Age") +
  ylab("Meth PC3") +
  geom_smooth(aes(x = Age, y = PC3), method = "lm", se = FALSE, size = line_size_blacks, color = "black") +
  guides(color = guide_legend(title = "CD8 Pheno")) + 
  theme_bw(base_size = font_size_scatter) +
  theme(legend.position = c(0.8, 0.9),
        legend.key.size = unit(0.1,"cm")) +
  labs(tag = "E")


pltE <- plt
ggsave("./Plots/Revision_plots/Figure1/PC3_Age_CD8diff.pdf", plt, width = 10, height = 10, units = "cm")


# Combine plot ------------------------------------------------------------


design <- "
AAAABBCC
AAAADDEE
AAAAFFGG
"

plt <- pltA + pltB + pltC + pltD + pltE + pltF + pltG + plot_layout(design = design) &
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave("./Plots/Revision_plots/Figure1/Figure_1.png", plt, width = 172, height = 172, units = "mm", dpi = 500)
ggsave("./Plots/Revision_plots/Figure1/Figure_1.pdf", plt, width = 172, height = 172, units = "mm")
