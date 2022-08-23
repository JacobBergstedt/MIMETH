

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(patchwork)
library(scales)
library(ggrepel)
library(latex2exp)
source("./Scripts/R_scripts/Plot_scripts/locuszoom.R")
source("./Scripts/R_scripts/Libraries/functions_for_plotting.R")


# Load annotation ---------------------------------------------------------

meth_anno_roadmap <- get_anno_meth_roadmap()
meth_anno_geography <- get_anno_meth_geography()
meth_anno_location <- get_anno_meth_location()
covs <- get_covs_884()
beta_values <- get_beta_values()
beta_values <- beta_values[match(covs$SUBJID, beta_values$SUBJID),]

# Load result data --------------------------------------------------------

ewas_effect_decomposition <- readRDS("./Data/RData/Results/Cell_mediation/mediation.rds")

ewas_med_sex <- readRDS("./Data/Chunk_data/Results/Cell_mediation/m_values_884/Sex.rds") %>%
  filter(Effect == "Mediation") %>%
  mutate(P_bonferroni = p.adjust(P_value)) %>%
  left_join(meth_anno_location) %>%
  left_join(select(meth_anno_roadmap, Probe, State_labels)) %>%
  left_join(meth_anno_geography) %>%
  select(Probe, Estimate, Standard_error, P_value, P_bonferroni, Probe_chr, Probe_position, Probe_gene, State_labels, Geography_labels) %>%
  arrange(P_value)

ewas_med_sex_sign <- ewas_med_sex %>%
  filter(P_bonferroni < 0.05) %>%
  mutate(Is_outlier = is_outlier(Estimate))


ewas_sex <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Sex.rds")
ewas_sex <- left_join(ewas_sex, meth_anno_location) %>%
  left_join(select(meth_anno_roadmap, Probe, State_labels)) %>%
  left_join(meth_anno_geography) %>%
  select(Probe, Estimate, Standard_error, P_value, P_bonferroni, Probe_chr, Probe_position, Probe_gene, State_labels, Geography_labels) %>%
  arrange(P_value)

ewas_sign_sex <- filter(ewas_sex, P_bonferroni < 0.05) %>%
  mutate(Is_outlier = is_outlier(Estimate))

ewas_no_corr_sex <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_no_cells/Sex.rds") %>%
  left_join(meth_anno_location) %>%
  left_join(meth_anno_roadmap) %>%
  left_join(meth_anno_geography) %>%
  select(Probe, Estimate, Standard_error, P_value, P_bonferroni, Probe_chr, Probe_position, Probe_gene, State_labels, Geography_labels) %>%
  arrange(P_value)

state_enrichments <- readRDS("./Data/RData/Results/EWAS_chromatin_state_enrichments_bonferroni_884.rds")
geography_enrichments <- readRDS("./Data/RData/Results/EWAS_geographic_region_enrichments_bonferroni_884.rds") %>%
  dplyr::rename(Term = Geographic_region_labels)
TFBS_enrichments <- readRDS("./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_dir_EWAS_884_bonferroni.rds")
TFBS_med_enrichments <- readRDS("./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_med_EWAS_884_bonferroni.rds")


# Plots 

# A -----------------------------------------------------------------------

plt_frame <- ewas_med_sex %>% 
  select(Probe, Mediated = Estimate, P_value, P_med = P_bonferroni) %>% 
  left_join(select(ewas_no_corr_sex, Probe, Probe_gene, Total = Estimate, P_total = P_bonferroni)) %>% 
  mutate(Probe_gene = simplify_probe_genes(Probe_gene)) %>% 
  filter(P_total < 0.05 | P_med < 0.05) %>%
  arrange(desc(abs(Total)))


pltA <-  plt_frame %>% 
  ggplot(aes(x = Mediated, y = Total)) +
  geom_point(color = col_age, size = point_size_scatter) +
  xlim(c(-4, 4)) +
  ylim(c(-4, 4)) +
  xlab("Mediated effect size") +
  ylab("Total effect size") +
  theme_bw(base_size = font_size_scatter) +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  labs(tag = "A")

pltA

# B -----------------------------------------------------------------------

plt_frame <- tibble(Sex = covs$Sex, CpG = beta_values[["cg09516963"]])
pltB <- ggplot(plt_frame, aes(x = Sex, y = CpG)) +
  geom_violin(col = col_sex) +
  ylab(expression(paste("5mC (", italic("DYRK2"), ")"))) +
  theme_bw(base_size = font_size_scatter) + 
  theme(panel.grid = element_blank())


pltB


# C -----------------------------------------------------------------------

thresh <- 5000
ref <- get_transcript_tib(chr = "chr12", 
                          start = 68039028, stop = 68059927, thresh = thresh)

ref_sub <- ref %>% filter(Transcript_ID == "NM_003583")

plt_locuszoom <- plot_locus_zoom(ref_sub, thresh = thresh, font_size_label = font_size_table)

ewas_sex_sub <- ewas_sex %>% 
  filter(Probe_chr == 12, Probe_position > plt_locuszoom$lim_min, Probe_position < plt_locuszoom$lim_max)

plt_locuszoom_scatter <- ggplot(ewas_sex_sub, aes(x = Probe_position / 1e3, y = Estimate)) +
  geom_point(col = col_age, size = point_size_scatter) + 
  scale_x_continuous(breaks = breaks_pretty(n = 4), limits = c(plt_locuszoom$lim_min / 1e3, plt_locuszoom$lim_max / 1e3)) + 
  theme_classic(base_size = font_size_scatter) +
  theme(axis.text = element_text(color = "black")) +
  xlab("Position (KB)") +
  ylab("Direct effect size")

pltC <- plt_locuszoom$plt + plt_locuszoom_scatter + plot_layout(heights = c(1, 6))
pltC[[1]] <- pltC[[1]] + labs(tag = "C")

pltC

# D -----------------------------------------------------------------------

plt_frame <- geography_enrichments %>% 
  filter(Type == "Direct", Variable == "Sex")

pltD <- plot_odds_ratio_split(plt_frame,  
                              NULL,
                              point_size = point_size_volcano,
                              bar_size = odds_ratio_errorbar_size,
                              colors = c(blues[4], col_sex),
                              base_size = font_size_odds_ratio,
                              legend_position = c(0.78, 0.88)) +
  labs(tag = "D")


pltD



# E -----------------------------------------------------------------------

pltE <- plot_effect_distribution_over_state(ewas_sign_sex,
                                            title = NULL,
                                            states = "Geography_labels",
                                            base_size = font_size_violin,
                                            line_size = line_size_violin,
                                            outlier_point_size = odds_ratio_errorbar_size,
                                            fontsize_count =  font_size_table,
                                            color = col_sex)

pltE[[1]] <- pltE[[1]] +
  labs(tag = "E") + 
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())


pltE



# F -----------------------------------------------------------------------

plt_frame <- geography_enrichments %>% 
  filter(Type == "Mediated", Variable == "Age")

pltF <- plot_odds_ratio_split(plt_frame,  
                              NULL,
                              point_size = point_size_volcano,
                              bar_size = odds_ratio_errorbar_size,
                              colors = c(blues[4], col_sex),
                              base_size = font_size_odds_ratio,
                              legend_position = c(0.85, 0.2)) +
  labs(tag = "F") + 
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

pltF



# G -----------------------------------------------------------------------

pltG <- plot_effect_distribution_over_state(ewas_med_sex_sign,
                                            title = NULL,
                                            states = "Geography_labels",
                                            base_size = font_size_violin,
                                            line_size = line_size_violin,
                                            outlier_point_size = odds_ratio_errorbar_size,
                                            fontsize_count =  font_size_table,
                                            color = col_sex)

pltG[[1]] <- pltG[[1]] +
  labs(tag = "G") + 
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())


pltG



# H ---------------------------------------------------------------------

plt_frame <- state_enrichments %>% 
  filter(Type == "Direct", Variable == "Sex", Cell == "Mononuclear")

pltH <- plot_odds_ratio_split(plt_frame, 
                              NULL,
                              point_size = point_size_volcano,
                              bar_size = odds_ratio_errorbar_size,
                              colors = c(blues[4], col_sex),
                              base_size = font_size_odds_ratio,
                              legend_position = c(0.78, 0.88)) +
  labs(tag = "H")

pltH



# I ---------------------------------------------------------------------

pltI <- plot_effect_distribution_over_state(ewas_sign_sex,
                                            title = NULL,
                                            states = "State_labels",
                                            base_size = font_size_violin,
                                            line_size = line_size_violin,
                                            outlier_point_size = odds_ratio_errorbar_size,
                                            fontsize_count =  font_size_table,
                                            color = col_sex)

pltI[[1]] <- pltI[[1]] +
  labs(tag = "I") + 
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())


pltI


# J ---------------------------------------------------------------------

plt_frame <- state_enrichments %>% 
  filter(Type == "Mediated", Variable == "Sex", Cell == "Mononuclear")

pltJ <- plot_odds_ratio_split(plt_frame, 
                              NULL,
                              point_size = point_size_volcano,
                              bar_size = odds_ratio_errorbar_size,
                              colors = c(blues[4], col_sex),
                              base_size = font_size_odds_ratio,
                              legend_position = c(0.80, 0.10)) +
  labs(tag = "J") + 
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

pltJ



# K ---------------------------------------------------------------------

pltK <- plot_effect_distribution_over_state(ewas_med_sex_sign,
                                            title = NULL,
                                            states = "State_labels",
                                            base_size = font_size_violin,
                                            line_size = line_size_violin,
                                            outlier_point_size = odds_ratio_errorbar_size,
                                            fontsize_count =  font_size_table,
                                            color = col_sex)
pltK[[1]] <- pltK[[1]] + 
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(tag = "K")


pltK




# Assemble plot -----------------------------------------------------------

# design <- "
# AAAABBCCCCCC
# DDEEEEFFGGGG
# HHIIIIJJKKKK
# "

design <- "
ABCCCC
DEEFGG
HIIJKK
"

plt <- pltA + pltB + pltC + pltD + pltE + pltF + pltG + pltH + pltI + pltJ + pltK + plot_layout(design = design)
ggsave("./Plots/Revision_plots/Supp_fig7/Supp_fig7.png", plt, width = 250, height = 200, units = "mm", dpi = 500)
ggsave("./Plots/Revision_plots/Supp_fig7/Supp_fig7.pdf", plt, width = 250, height = 200, units = "mm")
