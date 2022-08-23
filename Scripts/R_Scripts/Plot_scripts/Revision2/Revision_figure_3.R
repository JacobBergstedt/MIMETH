

# Load libraries ----------------------------------------------------------


library(tidyverse)
library(patchwork)
library(scales)
library(ggrepel)
source("./Scripts/MiscFunctions.R")
source("./Scripts/R_scripts/Plot_scripts/locuszoom.R")
source("./Scripts/R_scripts/Libraries/functions_for_plotting.R")


# Load annotation ---------------------------------------------------------


meth_anno_roadmap <- get_anno_meth_roadmap_all_cells() %>%
  mutate(State_labels = Chromatin_states_Mononuclear_cells) %>%
  mutate(State_labels = factor(State_labels, names(get_anno_roadmap_translation()))) %>%
  select(!contains("Chromatin_states"))

meth_anno_geography <- get_anno_meth_geography()
meth_anno_location <- get_anno_meth_location()

beta_values <- get_beta_values()
m_values <- get_m_values()
covs <- get_covs_884()
CD4_cells <- readRDS("./Data/RData/Cells/proportions_CD4_subsets.rds")
CD8_cells <- readRDS("./Data/RData/Cells/proportions_CD8_subsets.rds")

# Load result data --------------------------------------------------------

ewas_med_age_beta <- readRDS("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/Age.rds") %>%
  filter(Effect == "Mediation") %>%
  left_join(meth_anno_location, by = "Probe") %>%
  left_join(meth_anno_roadmap, by = "Probe") %>%
  left_join(meth_anno_geography, by = "Probe") %>%
  mutate(P_bonferroni = p.adjust(P_value))

ewas_med_age <- readRDS("./Data/Chunk_data/Results/Cell_mediation/m_values_884/Age.rds") %>%
  filter(Effect == "Mediation") %>%
  left_join(meth_anno_location, by = "Probe") %>%
  left_join(meth_anno_roadmap, by = "Probe") %>%
  left_join(meth_anno_geography, by = "Probe") %>%
  mutate(P_bonferroni = p.adjust(P_value))

ewas_med_age_sign <- ewas_med_age %>%
  filter(P_bonferroni < 0.05) %>%
  mutate(Is_outlier = is_outlier(Estimate))

ewas_age <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Age.rds")
ewas_age <- left_join(ewas_age, meth_anno_location) %>%
  left_join(select(meth_anno_roadmap, Probe, State_labels)) %>%
  left_join(meth_anno_geography) %>%
  select(Probe, Estimate, Standard_error, P_value, P_bonferroni, Probe_chr, Probe_position, Probe_gene, State_labels, Geography_labels) %>%
  arrange(P_value)

ewas_sign_age <- filter(ewas_age, P_bonferroni < 0.05) %>%
  mutate(Is_outlier = is_outlier(Estimate))

ewas_no_corr_age <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_no_cells/Age.rds")
ewas_no_corr_age <- left_join(ewas_no_corr_age, meth_anno_location) %>%
  left_join(select(meth_anno_roadmap, Probe, State_labels)) %>%
  left_join(meth_anno_geography) %>%
  select(Probe, Estimate, Standard_error, P_value, P_bonferroni, Probe_chr, Probe_position, Probe_gene, State_labels, Geography_labels) %>%
  arrange(P_value)


state_enrichments <- readRDS("./Data/RData/Results/EWAS_chromatin_state_enrichments_bonferroni_884.rds")
TFBS_enrichments <- readRDS("./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_dir_EWAS_884_bonferroni.rds")
TFBS_med_enrichments <- readRDS("./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_med_EWAS_884_bonferroni.rds")


disp <- readRDS("./Data/RData/Results/Dispersion/res_dispersion_gamlss_884.rds") %>%
  mutate(P_bonferroni = p.adjust(P_value)) %>%
  left_join(meth_anno_roadmap) %>%
  mutate(Is_outlier = is_outlier(Estimate))

disp_no_cells <- readRDS("./Data/RData/Results/Dispersion/res_dispersion_gamlss_no_cells_884.rds") %>%
  mutate(P_bonferroni = p.adjust(P_value)) %>%
  left_join(meth_anno_roadmap) %>%
  mutate(Is_outlier = is_outlier(Estimate))


disp_adj_var <- readRDS("./Data/RData/Results/Dispersion/res_dispersion_gamlss_variance_adjusted_for_cells_884.rds") %>%
  mutate(P_bonferroni = p.adjust(P_value)) %>%
  left_join(meth_anno_roadmap) %>%
  mutate(Is_outlier = is_outlier(Estimate))



#  A ---------------------------------------------------------------------


plt_frame <- ewas_med_age %>% 
  select(Probe, Effect, Mediated = Estimate, P_value, P_med = P_bonferroni) %>% 
  left_join(select(ewas_no_corr_age, Probe, Probe_gene, Total = Estimate, P_total = P_bonferroni)) %>% 
  mutate(Probe_gene = simplify_probe_genes(Probe_gene)) %>% 
  filter(P_total < 0.05) %>%
  arrange(desc(abs(Total)))

cpgs <- c("cg04875128", "cg13108341", "cg16867657", "cg06419846", "cg12480547", "cg27111050", "cg01877366")
tib_text <- plt_frame %>% 
  filter(Probe %in% cpgs) %>%
  mutate(nudge_x = 1.4 * ifelse(Mediated > 0, 1, -1),
         nudge_y = 0.2 * ifelse(Total > 0, 1, -1))

font_size = font_size_table * (1 / ggplot2:::.pt)

pltA <-  plt_frame %>% 
  ggplot(aes(x = Mediated, y = Total)) +
  geom_point(color = col_age, size = point_size_scatter) +
  geom_label_repel(aes(x = Mediated, y = Total, label = Probe_gene), data = tib_text, 
                   label.padding = 0.15,
                   size = 1, 
                   nudge_x = tib_text$nudge_x,
                   nudge_y = tib_text$nudge_y) +
  xlim(c(-4, 4)) +
  ylim(c(-4, 4)) +
  xlab("Mediated effect size") +
  ylab("Total effect size") +
  theme_bw(base_size = font_size_scatter) +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  labs(tag = "A")

ggsave("./Plots/Revision2_plots/Figure3/Fig3A_total_v_med.pdf", pltA, width = 6, height = 4, units = "cm")

#  B ---------------------------------------------------------------------

plt_frame <- state_enrichments %>% 
  filter(Type == "Direct", Variable == "Age", Cell == "CD4 naive")

pltB <- plot_odds_ratio_split(plt_frame, 
                              NULL,
                              point_size = point_size_volcano,
                              bar_size = odds_ratio_errorbar_size,
                              colors = c(blues[4], col_age),
                              base_size = font_size_odds_ratio,
                              legend_position = c(0.78, 0.88)) +
  labs(tag = "B")


ggsave("./Plots/Revision2_plots/Figure3/Fig3B_effect_size_enrich_dir.pdf", pltB, width = 6, height = 4, units = "cm")

#  C ---------------------------------------------------------------------

pltC <- plot_effect_distribution_over_state(ewas_sign_age,
                                            title = NULL,
                                            states = "State_labels",
                                            base_size = font_size_violin,
                                            line_size = line_size_violin,
                                            outlier_point_size = odds_ratio_errorbar_size,
                                            fontsize_count =  font_size_table,
                                            color = col_age)

pltC[[1]] <- pltC[[1]] +
  labs(tag = "C") + 
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave("./Plots/Revision2_plots/Figure3/Fig3C_effect_size_dist_dir.pdf", pltC, width = 6, height = 4, units = "cm")


#  D ---------------------------------------------------------------------

plt_frame <- TFBS_enrichments %>% 
  filter(Variable == "Age", Direction == "Pos") %>% 
  top_n(n = 15, conf.low) %>% 
  arrange(desc(conf.low)) %>%
  mutate(TF = gsub("EZH2phosphoT487", "EZH2phos.", TF)) %>% 
  mutate(Term = factor(TF, TF)) %>% 
  dplyr::rename(Estimate = estimate, Low = conf.low, High = conf.high)



pltD <- plot_odds_ratio(plt_frame, 
                        NULL, 
                        bar_size = odds_ratio_errorbar_size, 
                        color = col_age, 
                        base_size = font_size_odds_ratio,
                        limits = c(0, 22),
                        breaks = seq(0, 22, 2)) +
  labs(tag = "D")

ggsave("./Plots/Revision2_plots/Figure3/Fig3D_TFBS_pos_dir.pdf", pltD, width = 6, height = 4, units = "cm")


#  E ---------------------------------------------------------------------
plt_frame <- covs %>% inner_join(CD4_cells)

pltE <- plt_frame %>% 
  ggplot(aes(x = Age, y = X_CD4_EMRA_of_CD4.panel1)) + 
  geom_point(size = point_size_scatter, color = col_cell) + 
  xlab("Age") +
  ylab("% EMRA cells in CD4 T cells") +
  theme_bw(base_size = font_size_scatter) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(tag = "E")
  

ggsave("./Plots/Revision2_plots/Figure3/Fig3E_Age_v_CD4EMRA.pdf", pltE, width = 6, height = 4, units = "cm")


#  F ---------------------------------------------------------------------

plt_frame <- disp %>% 
  filter(P_bonferroni < 0.05) %>% 
  mutate(Direction = ifelse(Estimate > 0, "Pos.", "Neg.")) %>% 
  mutate(Direction = factor(Direction, c("Pos.", "Neg.")))

pltF <- ggplot(plt_frame, aes(x = Direction)) + 
  geom_bar(fill = col_age, col = col_age) +
  ylab("# assoc. CpG sites") +
  xlab("Direction") + 
  theme_classic(base_size = font_size_scatter) +
  theme(axis.text = element_text(color = "black")) +
  labs(tag = "F")


ggsave("./Plots/Revision2_plots/Figure3/Fig3F_disp_bar_plot.pdf", pltF, width = 6, height = 4, units = "cm")

# 3 G ---------------------------------------------------------------------

plt_frame <- inner_join(select(ewas_age, Probe, Direct = Estimate, P_dir = P_bonferroni), 
                        select(disp, Probe, Disp = Estimate, P_disp = P_bonferroni))
pltG <- plt_frame %>% 
  filter(P_disp < 0.05) %>% 
  ggplot(aes(x = Disp, y = Direct)) +
  geom_point(color = col_age, size = point_size_scatter) +
  xlab("Dispersion") +
  ylab("Direct effect") +
  theme_bw(base_size = font_size_scatter) +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  labs(tag = "G")

ggsave("./Plots/Revision2_plots/Figure3/Fig3G_ewas_v_disp.pdf", pltI, width = 6, height = 4, units = "cm")



# 3 H ---------------------------------------------------------------------

plt_frame <- state_enrichments %>% 
  filter(Type == "Mediated", Variable == "Age", Cell == "CD4 naive")

pltH <- plot_odds_ratio_split(plt_frame, 
                              NULL,
                              point_size = point_size_volcano,
                              bar_size = odds_ratio_errorbar_size,
                              colors = c(blues[4], col_age),
                              base_size = font_size_odds_ratio,
                              legend_position = c(0.80, 0.10)) +
  labs(tag = "H")

ggsave("./Plots/Revision2_plots/Figure3/Fig3H_med_state_enrich.pdf", pltI, width = 6, height = 4, units = "cm")



# 3 I ---------------------------------------------------------------------

pltI <- plot_effect_distribution_over_state(ewas_med_age_sign,
                                            title = NULL,
                                            states = "State_labels",
                                            base_size = font_size_violin,
                                            line_size = line_size_violin,
                                            outlier_point_size = odds_ratio_errorbar_size,
                                            fontsize_count =  font_size_table,
                                            color = col_age)
pltI[[1]] <- pltI[[1]] + 
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(tag = "I")

ggsave("./Plots/Revision2_plots/Figure3/Fig3I_med_dist.pdf", pltI, width = 6, height = 4, units = "cm")


# 3 J ---------------------------------------------------------------------

plt_frame <- TFBS_med_enrichments %>% 
  filter(Variable == "Age", Direction == "Pos") %>% 
  top_n(n = 15, conf.low) %>% 
  arrange(desc(conf.low)) %>% 
  mutate(Term = factor(TF, TF)) %>% 
  dplyr::rename(Estimate = estimate, Low = conf.low, High = conf.high)

pltJ <- plot_odds_ratio(plt_frame,
                        NULL,
                        bar_size = odds_ratio_errorbar_size, 
                        color = col_age, 
                        base_size = font_size_odds_ratio,
                        limits = c(0, 22),
                        breaks = seq(0, 22, 2)) +
  labs(tag = "J")

ggsave("./Plots/Revision2_plots/Figure3/Fig3I_med_TFBS_enrich.pdf", pltI, width = 6, height = 4, units = "cm")

# Assemble figure ---------------------------------------------------------


design <- "
  AAABCCC
  DDEEFGG
  HIIIJJJ
"

plt <- (pltA + pltB + pltC + pltD + pltE + pltF + pltG + pltH + pltI + pltJ +
          plot_layout(design = design)) &
  theme(plot.tag = element_text(size = 6, face = "bold"))

ggsave("./Plots/Revision2_plots/Figure3/Figure_3.png", plt, width = 172, height = 147, units = "mm", dpi = 500)
ggsave("./Plots/Revision2_plots/Figure3/Figure_3.pdf", plt, width = 172, height = 147, units = "mm")
