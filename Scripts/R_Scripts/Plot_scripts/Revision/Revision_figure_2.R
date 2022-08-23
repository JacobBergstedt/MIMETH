

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(parallel)
library(patchwork)
library(scales)
library(ggrepel)
library(splines)
library(broom)
library(robCompositions)
library(readxl)
source("./Scripts/MiscFunctions.R")
source("./Scripts/R_scripts/Libraries/functions_for_compositions.R")
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
source("./Scripts/R_scripts/Plot_scripts/locuszoom.R")
source("./Scripts/R_scripts/Libraries/functions_for_plotting.R")


# Load data ---------------------------------------------------------------


balances <- readRDS("./Data/RData/Results/EWAS/M_values/Environment/EWAS_balance_effects_compositional_analysis_partition3.rds") %>%
  group_by(Balance) %>%
  mutate(P_bonferroni = p.adjust(P_value)) %>%
  ungroup() %>%
  left_join(get_anno_meth_location()) %>%
  left_join(get_anno_meth_roadmap())

CD8Diff_CD8Naive <- balances %>% filter(Balance == "CD8Diff_CD8Naive")
CD4Diff_CD4Naive <- balances %>% filter(Balance == "CD4Diff_CD4Naive")

meth <- get_beta_values()
m_values <- get_m_values()
covs <- get_covs_884()
meth <- meth[match(covs$SUBJID, meth$SUBJID),]
m_values <- m_values[match(covs$SUBJID, m_values$SUBJID),]

meth_anno_roadmap <- get_anno_meth_roadmap_all_cells() %>%
  mutate(State_labels = Chromatin_states_Mononuclear_cells) %>%
  mutate(State_labels = factor(State_labels, names(get_anno_roadmap_translation()))) %>%
  select(!contains("Chromatin_states"))

meth_anno_geography <- get_anno_meth_geography()
meth_anno_location <- get_anno_meth_location()

ewas_med_CMV_beta <- readRDS("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/CMV.rds") %>%
  filter(Effect == "Mediation") %>%
  left_join(meth_anno_location, by = "Probe") %>%
  left_join(meth_anno_roadmap, by = "Probe") %>%
  left_join(meth_anno_geography, by = "Probe") %>%
  mutate(P_bonferroni = p.adjust(P_value))

ewas_med_CMV <- readRDS("./Data/Chunk_data/Results/Cell_mediation/m_values_884/CMV.rds") %>%
  filter(Effect == "Mediation") %>%
  left_join(meth_anno_location, by = "Probe") %>%
  left_join(meth_anno_roadmap, by = "Probe") %>%
  left_join(meth_anno_geography, by = "Probe") %>%
  mutate(P_bonferroni = p.adjust(P_value))

ewas_med_CMV_sign <- ewas_med_CMV %>%
  filter(P_bonferroni < 0.05) %>%
  mutate(Is_outlier = is_outlier(Estimate))

ewas_CMV <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/CMV_serostatus.rds")
ewas_CMV <- left_join(ewas_CMV, meth_anno_location) %>%
  left_join(select(meth_anno_roadmap, Probe, State_labels)) %>%
  left_join(meth_anno_geography) %>%
  select(Probe, Estimate, Standard_error, P_value, P_bonferroni, Probe_chr, Probe_position, Probe_gene, State_labels, Geography_labels) %>%
  arrange(P_value)

ewas_sign_CMV <- filter(ewas_CMV, P_bonferroni < 0.05) %>%
  mutate(Is_outlier = is_outlier(Estimate))

ewas_no_corr_CMV <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_no_cells/CMV_serostatus.rds")
ewas_no_corr_CMV <- left_join(ewas_no_corr_CMV, meth_anno_location) %>%
  left_join(select(meth_anno_roadmap, Probe, State_labels)) %>%
  left_join(meth_anno_geography) %>%
  select(Probe, Estimate, Standard_error, P_value, P_bonferroni, Probe_chr, Probe_position, Probe_gene, State_labels, Geography_labels) %>%
  arrange(P_value)

cell_list_16_cells <- get_cell_list()
cell_list <- cell_list_16_cells[!cell_list_16_cells == "X_dendritic_cells.panel8"]

cells <- covs %>%
  select(SUBJID, all_of(cell_list))

cell_balances <- add_column(balance_preds(cells[-1], sbp = read_xlsx("./Data/Cell_subset_partitions_3.xlsx")),
                            SUBJID = cells[["SUBJID"]], .before = 1)


state_enrichments <- readRDS("./Data/RData/Results/EWAS_chromatin_state_enrichments_bonferroni_884.rds")
TFBS_enrichments <- readRDS("./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_dir_EWAS_884_bonferroni.rds")
TFBS_med_enrichments <- readRDS("./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_med_EWAS_884_bonferroni.rds")


compute_correlation <- function(x, y) {

  cor(x$Estimate, y$Estimate)

}

f <- factor(balances$Balance)
balances %>%
  split(f = f) %>%
  map_dbl( ~ compute_correlation(x = ., y = ewas_med_CMV))




# A ----------------------------------------------------------------------


plt_frame <- ewas_med_CMV %>%
  select(Probe, Med = Estimate, P_med = P_bonferroni) %>%
  left_join(select(ewas_no_corr_CMV, Probe, Probe_gene, Total = Estimate, P_total = P_bonferroni)) %>%
  left_join(select(ewas_CMV, Probe, Probe_gene, Direct = Estimate, P_direct = P_bonferroni)) %>%
  mutate(Probe_gene = simplify_probe_genes(Probe_gene)) %>%
  arrange(desc(abs(Total)))

#
# cpgs <- c("cg04875128", "cg13108341", "cg16867657", "cg06419846", "cg12480547", "cg27111050", "cg01877366")
# tib_text <- plt_frame %>%
#   filter(Probe %in% cpgs) %>%
#   mutate(nudge_x = 1.4 * ifelse(Mediated > 0, 1, -1),
#          nudge_y = 0.2 * ifelse(Total > 0, 1, -1))
#
# font_size = font_size_table * (1 / ggplot2:::.pt)

pltA <-  plt_frame %>%
  filter(P_total < 0.05) %>% 
  ggplot(aes(x = Med, y = Total)) +
  geom_point(color = col_CMV, size = point_size_scatter) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_x_continuous(limits = c(-1, 1)) +
  xlab("Direct effect size") +
  ylab("Total effect size") +
  theme_bw(base_size = font_size_scatter) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(tag = "A")


pltB <-  plt_frame %>%
  filter(P_total < 0.05) %>% 
  ggplot(aes(x = Direct, y = Total)) +
  geom_point(color = col_CMV, size = point_size_scatter) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_x_continuous(limits = c(-1, 1)) +
  xlab("Direct effect size") +
  ylab("Total effect size") +
  theme_bw(base_size = font_size_scatter) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(tag = "A")

pltC <-  plt_frame %>%
  filter(P_direct < 0.05) %>% 
  ggplot(aes(x = Direct, y = Total)) +
  geom_point(color = col_CMV, size = point_size_scatter) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_x_continuous(limits = c(-1, 1)) +
  xlab("Direct effect size") +
  ylab("Total effect size") +
  theme_bw(base_size = font_size_scatter) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(tag = "A")



ggsave("./Plots/Revision_plots/Figure2/Fig2A_CMV_Total_v_Med.pdf", pltA, width = 6, height = 4, units = "cm")
ggsave("./Plots/Revision_plots/Figure2/Fig2A_CMV_Total_v_Dir.pdf", pltB, width = 6, height = 4, units = "cm")
ggsave("./Plots/Revision_plots/Figure2/Fig2A_CMV_Total_v_Dir_sign_Dir.pdf", pltC, width = 6, height = 4, units = "cm")


# B -----------------------------------------------------------------------

plt_frame <- ewas_med_CMV %>% 
  select(Probe, Estimate_med = Estimate, P_med = P_bonferroni) %>%
  left_join(select(CD8Diff_CD8Naive, Probe, Estimate_CD8Diff_CD8Naive = Estimate, P_CD8Diff_CD8Naive = P_bonferroni)) %>% 
  left_join(select(CD4Diff_CD4Naive, Probe, Estimate_CD4Diff_CD4Naive = Estimate, P_CD4Diff_CD4Naive = P_bonferroni))


pltB <- plt_frame %>%
  filter(P_med < 0.05) %>% 
  ggplot(aes(y = Estimate_med, x = Estimate_CD8Diff_CD8Naive)) +
  geom_point(size = point_size_scatter, color = col_CMV) +
  xlim(-0.5, 0.5) +
  ylim(-0.5, 0.5) +
  ylab("Mediated effect") +
  xlab("CD8 Diff. effect") +
  theme_bw(base_size = font_size_scatter) +
  labs(tag = "B")

ggsave("./Plots/Revision_plots/Figure2/scatter_CMV_med_v_CD8Diff_balance.pdf", pltB, width = 10, height = 10, units = "cm")

# C & D -----------------------------------------------------------------------

plt_frame <- inner_join(select(covs, SUBJID, CMV_serostatus), cell_balances)
plt_frame <- add_column(beta = meth[["cg15302376"]], plt_frame)


pltC1 <- plt_frame %>%
  ggplot(aes(x = CMV_serostatus, y = CD4Diff_CD4Naive, fill = CMV_serostatus)) +
  geom_violin(color = "black", size = line_size_scatter) +
  scale_fill_manual(values = c(reds[3], col_CMV)) +
  ylim(-4, 3) +
  ylab("CD4 differentiation skew") +
  theme_classic(base_size = font_size_scatter) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none")

ggsave("./Plots/Revision_plots/Figure2/violin_CMV_v_CD4Diff_m_values.pdf", pltC, width = 10, height = 10, units = "cm")


pltC2 <- plt_frame %>%
  ggplot(aes(x = beta, y = CD4Diff_CD4Naive)) +
  geom_point(size = point_size_scatter, col = col_cell) +
  geom_smooth(method = "lm", size = line_size_scatter, se = FALSE, color = "black") +    
  ylim(-4, 3) +
  ylab(NULL) +
  theme_classic(base_size = font_size_scatter) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

ggsave("./Plots/Revision_plots/Figure2/scatter_CD4Diff_m_values.pdf", pltC2, width = 10, height = 10, units = "cm")


pltC <- pltC1 + pltC2 + plot_layout(widths = c(1, 3))
pltC[[1]] <- pltC[[1]] + labs(tag = "C")

ggsave("./Plots/Revision_plots/Figure2/mediation_DNMT3A.pdf", pltC, width = 10, height = 10, units = "cm")



# D -----------------------------------------------------------------------



plt_frame <- inner_join(select(covs, SUBJID, CMV_serostatus), cell_balances)
plt_frame <- add_column(beta = meth[["cg15302376"]], plt_frame)


pltD1 <- plt_frame %>%
  ggplot(aes(x = CMV_serostatus, y = CD8Diff_CD8Naive, fill = CMV_serostatus)) +
  geom_violin(color = "black", size = line_size_scatter) +
  scale_fill_manual(values = c(reds[3], col_CMV)) +
  ylim(-4, 3) +
  ylab("CD8 differentiation skew") +
  theme_classic(base_size = font_size_scatter) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none")
ggsave("./Plots/Revision_plots/Figure2/violin_CMV_v_CD8Diff_m_values.pdf", pltD1, width = 10, height = 10, units = "cm")


pltD2 <- plt_frame %>%
  ggplot(aes(x = beta, y = CD8Diff_CD8Naive)) +
  geom_point(size = point_size_scatter, col = col_cell) +
  geom_smooth(method = "lm", size = line_size_scatter, se = FALSE, color = "black") +    
  ylim(-4, 3) +
  ylab(NULL) +
  theme_classic(base_size = font_size_scatter) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

ggsave("./Plots/Revision_plots/Figure2/scatter_CD8Diff_m_values.pdf", pltD2, width = 10, height = 10, units = "cm")

pltD <- pltD1 + pltD2 + plot_layout(widths = c(1, 3))
pltD[[1]] <- pltD[[1]] + labs(tag = "D")

ggsave("./Plots/Revision_plots/Figure2/mediation_DNMT3A.pdf", pltD, width = 10, height = 10, units = "cm")



# E -----------------------------------------------------------------------


thresh <- 5000
ref <- get_transcript_tib(chr = "chr11", start = 65302349, stop = 65357476, thresh = thresh, only_major_transcript = TRUE)

ref_sub <- ref %>% filter(Transcript_ID == "NM_001130144")

plt_locuszoom <- plot_locus_zoom(ref_sub, thresh = thresh, font_size_label = font_size_table)

ewas_CMV_sub <- ewas_CMV %>% 
  filter(Probe_chr == 11, Probe_position > plt_locuszoom$lim_min, Probe_position < plt_locuszoom$lim_max)

plt_locuszoom_scatter <- ggplot(ewas_CMV_sub, aes(x = Probe_position / 1e3, y = Estimate)) +
  geom_point(col = col_CMV, size = point_size_scatter) + 
  scale_x_continuous(breaks = breaks_pretty(n = 4), limits = c(plt_locuszoom$lim_min / 1e3, plt_locuszoom$lim_max / 1e3)) +   
  theme_classic(base_size = font_size_scatter) +
  theme(axis.text = element_text(color = "black")) +
  xlab("Position (KB)") +
  ylab("Direct effect size")

pltE <- plt_locuszoom$plt + plt_locuszoom_scatter + plot_layout(heights = c(1, 6))
pltE[[1]] <- pltE[[1]] + labs(tag = "E")
ggsave("./Plots/Revision_plots/Figure2/locus_zoomDNMT3A.pdf", pltE, width = 10, height = 10, units = "cm")



# F -----------------------------------------------------------------------

plt_frame <- TFBS_enrichments %>%
  filter(Variable == "CMV", Direction == "Pos") %>%
  top_n(n = 15, conf.low) %>%
  arrange(desc(estimate)) %>%
  mutate(Term = factor(TF, TF)) %>%
  dplyr::rename(Estimate = estimate, Low = conf.low, High = conf.high)



pltF <- plot_odds_ratio(plt_frame,
                        NULL,
                        bar_size = odds_ratio_errorbar_size,
                        color = col_CMV,
                        base_size = font_size_odds_ratio,
                        limits = c(0, 400),
                        breaks = seq(0, 400, 20)) +
  labs(tag = "F")

ggsave("./Plots/Revision_plots/Figure2/TFBS_pos.pdf", pltF, width = 10, height = 10, units = "cm")


# G -----------------------------------------------------------------------

plt_frame <- TFBS_enrichments %>% 
  filter(Variable == "CMV", Direction == "Neg") %>% 
  top_n(n = 15, conf.low) %>% 
  arrange(desc(estimate)) %>%
  mutate(Term = factor(TF, TF)) %>% 
  dplyr::rename(Estimate = estimate, Low = conf.low, High = conf.high)


pltG <- plot_odds_ratio(plt_frame, 
                        NULL, 
                        bar_size = odds_ratio_errorbar_size, 
                        color = col_CMV, 
                        base_size = font_size_odds_ratio,
                        limits = c(0, 400),
                        breaks = seq(0, 400, 20)) +
  labs(tag = "G")


ggsave("./Plots/Revision_plots/Figure2/TFBS_neg.pdf", pltG, width = 10, height = 10, units = "cm")

# Assemble ----------------------------------------------------------------

design <- "
  AB
  CD
  EE
  FG
"

plt <- pltA + pltB + pltC + pltD + pltE + pltF + pltG + plot_layout(design = design) &
  theme(plot.tag = element_text(size = 6, face = "bold"))

ggsave("./Plots/Revision_plots/Figure2/Figure_2.png", plt, width = 172, height = 200, units = "mm", dpi = 500)
ggsave("./Plots/Revision_plots/Figure2/Figure_2.pdf", plt, width = 172, height = 200, units = "mm")
