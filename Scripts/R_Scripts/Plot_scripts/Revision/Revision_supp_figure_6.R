

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(patchwork)
library(scales)
library(ggrepel)
source("./Scripts/R_scripts/Libraries/functions_for_plotting.R")

# Load data ---------------------------------------------------------------

meth_anno_roadmap <- get_anno_meth_roadmap()
meth_anno_geography <- get_anno_meth_geography()
meth_anno_location <- get_anno_meth_location()

ewas_med_age <- readRDS("./Data/Chunk_data/Results/Cell_mediation/m_values_884/Age.rds") %>% 
  filter(Effect == "Mediation") %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  left_join(meth_anno_location) %>%
  left_join(select(meth_anno_roadmap, Probe, State_labels)) %>%
  left_join(meth_anno_geography) %>%
  select(Probe, Estimate, Standard_error, P_value, P_bonferroni, Probe_chr, Probe_position, Probe_gene, State_labels, Geography_labels) %>%
  arrange(P_value)

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

state_enrichments <- readRDS("./Data/RData/Results/EWAS_chromatin_state_enrichments_bonferroni_884.rds")
geography_enrichments <- readRDS("./Data/RData/Results/EWAS_geographic_region_enrichments_bonferroni_884.rds") %>%
  dplyr::rename(Term = Geographic_region_labels)
TFBS_enrichments <- readRDS("./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_dir_EWAS_884_bonferroni.rds")
TFBS_med_enrichments <- readRDS("./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_med_EWAS_884_bonferroni.rds")


# A -----------------------------------------------------------------------

plt_frame <- geography_enrichments %>%
  filter(Type == "Direct", Variable == "Age")

pltA <- plot_odds_ratio_split(plt_frame,
                              NULL,
                              point_size = point_size_volcano,
                              bar_size = odds_ratio_errorbar_size,
                              colors = c(reds[3], col_age),
                              base_size = font_size_odds_ratio,
                              legend_position = c(0.8, 0.88)) +
  labs(tag = "A")


# B -----------------------------------------------------------------------

pltB <- plot_effect_distribution_over_state(ewas_sign_age,
                                            title = NULL,
                                            states = "Geography_labels",
                                            base_size = font_size_violin,
                                            line_size = line_size_violin,
                                            outlier_point_size = odds_ratio_errorbar_size,
                                            fontsize_count =  font_size_table,
                                            color = col_age)

pltB[[1]] <- pltB[[1]] +
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(tag = "B")


# C ---------------------------------------------------------------------

plt_frame <- geography_enrichments %>%
  filter(Type == "Mediated", Variable == "Age")

pltC <- plot_odds_ratio_split(plt_frame,
                              NULL,
                              point_size = point_size_volcano,
                              bar_size = odds_ratio_errorbar_size,
                              colors = c(reds[3], col_age),
                              base_size = font_size_odds_ratio,
                              legend_position = c(0.8, 0.88)) +
  labs(tag = "C")



# D ---------------------------------------------------------------------

pltD <- plot_effect_distribution_over_state(ewas_med_age_sign,
                                            title = NULL,
                                            states = "Geography_labels",
                                            base_size = font_size_violin,
                                            line_size = line_size_violin,
                                            outlier_point_size = odds_ratio_errorbar_size,
                                            fontsize_count =  font_size_table,
                                            color = col_age)

pltD[[1]] <- pltD[[1]] +
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(tag = "D")


# Assemble plot -----------------------------------------------------------

design <- "
  ABB
  CDD
"


plt <- pltA + pltB + pltC + pltD + plot_layout(design = design) &
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave("./Plots/Revision_plots/Supp_fig6/Supp_fig6.png", plt, width = 172, height = 147, units = "mm", dpi = 500)
ggsave("./Plots/Revision_plots/Supp_fig6/Supp_fig6.pdf", plt, width = 172, height = 147, units = "mm")


