

# Load libraries ----------------------------------------------------------


library(tidyverse)
library(patchwork)
library(scales)
library(ggrepel)
library(readxl)
library(robCompositions)
source("./Scripts/R_scripts/Plot_scripts/locuszoom.R")
source("./Scripts/R_scripts/Libraries/functions_for_plotting.R")
source("./Scripts/R_scripts/Libraries/functions_for_annotation.R")
source("./Scripts/R_scripts/Libraries/functions_for_compositions.R")

# Load data ---------------------------------------------------------------

covs <- get_covs_884()
meth <- get_beta_values()
meth <- meth[match(covs$SUBJID, meth$SUBJID),]


meth_anno_roadmap <- get_anno_meth_roadmap_all_cells() %>%
  mutate(State_labels = Chromatin_states_CD4_naive) %>%
  mutate(State_labels = factor(State_labels, names(get_anno_roadmap_translation()))) %>%
  select(!contains("Chromatin_states"))

meth_anno_geography <- get_anno_meth_geography()
meth_anno_location <- get_anno_meth_location()

# Smoking
ewas_med_smoking <- readRDS("./Data/Chunk_data/Results/Cell_mediation/m_values_884/Smoking.rds") %>% 
  filter(Effect == "Mediation") %>% 
  mutate(P_bonferroni = p.adjust(P_value))

ewas_med_smoking_sign <- ewas_med_smoking %>%
  filter(P_bonferroni < 0.05) %>%
  mutate(Is_outlier = is_outlier(Estimate))

ewas_smoking <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Smoking_status.rds") %>%
  filter(Levels != "Ex_Smoker")

ewas_smoking <- left_join(ewas_smoking, meth_anno_location) %>%
  left_join(select(meth_anno_roadmap, Probe, State_labels)) %>%
  left_join(meth_anno_geography) %>%
  select(Probe, Estimate, Standard_error, P_value, P_bonferroni, Probe_chr, Probe_position, Probe_gene, State_labels, Geography_labels) %>%
  arrange(P_value)

ewas_sign_smoking <- filter(ewas_smoking, P_bonferroni < 0.05) %>%
  mutate(Is_outlier = is_outlier(Estimate))

ewas_beta_smoking <- readRDS("./Data/Chunk_data/Results/EWAS/Beta_values/Environment/Random_effect/Correct_for_16_props/Smoking_status.rds") %>%
  filter(Levels != "Ex_Smoker")

ewas_beta_smoking <- left_join(ewas_beta_smoking, meth_anno_location) %>%
  left_join(select(meth_anno_roadmap, Probe, State_labels)) %>%
  left_join(meth_anno_geography) %>%
  select(Probe, Estimate, Standard_error, P_value, P_bonferroni, Probe_chr, Probe_position, Probe_gene, State_labels, Geography_labels) %>%
  arrange(P_value)

ewas_no_corr_smoking <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect/Correct_for_no_cells/Smoking_status.rds") %>%
  filter(Levels != "Ex_Smoker")

ewas_no_corr_smoking <- left_join(ewas_no_corr_smoking, meth_anno_location) %>%
  left_join(select(meth_anno_roadmap, Probe, State_labels)) %>%
  left_join(meth_anno_geography) %>%
  select(Probe, Estimate, Standard_error, P_value, P_bonferroni, Probe_chr, Probe_position, Probe_gene, State_labels, Geography_labels) %>%
  arrange(P_value)

state_enrichments <- readRDS("./Data/RData/Results/EWAS_chromatin_state_enrichments_bonferroni_884.rds")
TFBS_enrichments <- readRDS("./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_dir_EWAS_884_bonferroni.rds")
TFBS_med_enrichments <- readRDS("./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_med_EWAS_884_bonferroni.rds")

# A ----------------------------------------------------------------------

# plt_frame <- ewas_effect_decomposition %>% 
#   filter(Variable == "Smoking", Effect != "Ratio_mediation") %>% 
#   pivot_wider(Probe, names_from = Effect, values_from = Estimate)

# 

plt_frame <- ewas_med_smoking %>%
  select(Probe, Mediated = Estimate, P_med = P_bonferroni) %>%
  left_join(select(ewas_no_corr_smoking, Probe, Probe_gene, Total = Estimate, P_total = P_bonferroni)) %>%
  mutate(Probe_gene = simplify_probe_genes(Probe_gene)) %>%
  filter(P_med < 0.05 | P_total < 0.05) %>% 
  arrange(desc(abs(Total))) 

# 
# cpgs <- c("cg04875128", "cg13108341", "cg16867657", "cg06419846", "cg12480547", "cg27111050", "cg01877366")
# tib_text <- plt_frame %>% 
#   filter(Probe %in% cpgs) %>%
#   mutate(nudge_x = 1.4 * ifelse(Mediated > 0, 1, -1),
#          nudge_y = 0.2 * ifelse(Total > 0, 1, -1))
# 
# font_size = font_size_table * (1 / ggplot2:::.pt)

# pltI <-  plt_frame %>% 
#   ggplot(aes(x = Mediation, y = Total)) +
#   geom_point(color = col_age, size = point_size_scatter) +
#   xlab("Mediation") +
#   ylab("Total") +
#   theme_bw(base_size = font_size_scatter) +
#   theme(panel.grid = element_blank(), 
#         plot.title = element_text(hjust = 0.5)) +
#   labs(tag = "I")


pltA <-  plt_frame %>% 
  ggplot(aes(x = Mediated, y = Total)) +
  geom_point(color = col_smoking, size = point_size_scatter) +
  scale_y_continuous(limits = c(-2.3, 2.3)) + 
  scale_x_continuous(limits = c(-2.3, 2.3)) + 
  xlab("Mediated effect size") +
  ylab("Total effect size") +
  theme_bw(base_size = font_size_scatter) +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  labs(tag = "A")

pltA

# B -----------------------------------------------------------------------

plt_frame <- TFBS_enrichments %>% 
  filter(Variable == "Smoking", Direction == "Neg") %>% 
  top_n(n = 15, conf.low) %>% 
  arrange(desc(conf.low)) %>%
  mutate(Term = factor(TF, TF)) %>% 
  dplyr::rename(Estimate = estimate, Low = conf.low, High = conf.high)



pltB <- plot_odds_ratio(plt_frame, 
                        NULL, 
                        bar_size = odds_ratio_errorbar_size, 
                        color = col_smoking, 
                        base_size = font_size_odds_ratio,
                        limits = c(0, 32),
                        breaks = seq(from = 0, to = 32, by = 4)) +
  labs(tag = "B")

pltB

# C -----------------------------------------------------------------------

cpg <- "cg11335172"
plt_frame <- tibble(CpG = meth[[cpg]],
                    CMV = covs$CMV_serostatus, 
                    NK = 100 * covs$X_NK_of_total.panel4,
                    Smoking = covs$Smoking_status) %>% 
  filter(Smoking != "Ex_Smoker") %>%
  mutate(Smoking = fct_drop(Smoking))

pltC <- ggplot(plt_frame, aes(x = Smoking, y = NK, fill = Smoking)) +
  geom_violin(size = line_size_scatter) +
  scale_fill_manual(values = c(reds[3], col_smoking)) +
  ylab("NK cells [%]") +
  xlab("Smoking status") +
  theme_classic(base_size = font_size_scatter) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  labs(tag = "C")

# D -----------------------------------------------------------------------

pltD <- ggplot(plt_frame, aes(x = CpG, y = NK, color = col_smoking)) +
  geom_point(size = point_size_scatter) +
  geom_smooth(method = "lm", size = line_size_scatter, se = FALSE) +  
  xlab(expression(paste("5mC (", italic("IL18RAP"), ")"))) +
  ylab(NULL) + 
  theme_classic(base_size = font_size_scatter) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  labs(tag = "D")


# E -----------------------------------------------------------------------
probes <- ewas_smoking %>% 
  filter(P_bonferroni <= 0.05) %>% 
  pull(Probe) %>% 
  unique()

years_since_last_smoke <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Years_since_last_smoke.rds")
years_smoking <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Years_smoking.rds") %>% 
  filter(Probe %in% probes)

plt_frame_years_since_last_smoke <- years_since_last_smoke %>% 
  right_join(tibble(Probe = probes))

plt_frame_years_smoking <- years_smoking %>% 
  right_join(tibble(Probe = probes))

plt_frame <- tibble(Probe = probes, 
                    Years_since_last_smoke = plt_frame_years_since_last_smoke$Estimate, 
                    Years_smoking = plt_frame_years_smoking$Estimate)

color = col_smoking
pltE <- plt_frame %>% 
  ggplot(aes(x = Years_since_last_smoke, y = Years_smoking, color = color)) +
  geom_point(size = point_size_scatter) +
  scale_x_continuous(limits = c(-2.6, 2.6)) +
  scale_y_continuous(limits = c(-2.6, 2.6)) +
  scale_color_manual(values = color) +
  xlab("Effect size, 20 years after quitting") +
  ylab("Effect size, 20 years of smoking") +
  theme_classic(base_size = font_size_scatter) +
  theme(axis.text = element_text(color = "black"), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

# F -----------------------------------------------------------------------

years_since_last_smoke_med <- readRDS("./Data/Chunk_data/Results/Cell_mediation/m_values/Years_since_last_smoke.rds") %>% filter(Effect == "Mediation")
years_smoking_med <- readRDS("./Data/Chunk_data/Results/Cell_mediation/m_values/Years_smoking.rds") %>% filter(Effect == "Mediation")

probes <- ewas_med_smoking %>% 
  filter(P_bonferroni <= 0.05) %>% 
  pull(Probe) %>% 
  unique()

plt_frame_years_since_last_smoke_med <- years_since_last_smoke_med %>% 
  right_join(tibble(Probe = probes))

plt_frame_years_smoking_med <- years_smoking_med %>% 
  right_join(tibble(Probe = probes))

plt_frame_med <- tibble(Probe = probes, 
                        Years_since_last_smoke_med = plt_frame_years_since_last_smoke_med$Estimate, 
                        Years_smoking_med = plt_frame_years_smoking_med$Estimate)

lim <- max(max(abs(plt_frame_med$Years_since_last_smoke_med)),
           max(abs(plt_frame_med$Years_smoking_med)))

color = col_smoking
pltF <- plt_frame_med %>% 
  ggplot(aes(x = Years_since_last_smoke_med, y = Years_smoking_med, color = color)) +
  geom_point(size = point_size_scatter) +
  scale_color_manual(values = color) +
  xlim(c(-lim, lim)) +
  ylim(c(-lim, lim)) + 
  xlab("Mediation effect, 20 years after quitting") +
  ylab("Mediation effect, 20 years of smoking") +
  theme_classic(base_size = font_size_scatter) +
  theme(axis.text = element_text(color = "black"), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))



# G -----------------------------------------------------------------------
ewas_CRP <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Log_CRP_levels.rds") %>% 
  left_join(meth_anno_location)

thresh <- 5000
ref <- get_transcript_tib(chr = "chr18", 
                          start = 60775216, 
                          stop = 61069268, 
                          thresh = thresh)

ref_sub <- ref %>% 
  filter(Transcript_ID == "NM_000633")

plt_locuszoom <- plot_locus_zoom(ref_sub, thresh = thresh, font_size_label = 3)

ewas_CRP_sub <- ewas_CRP %>% 
  filter(Probe_chr == "18", 
         Probe_position > plt_locuszoom$lim_min, 
         Probe_position < plt_locuszoom$lim_max)

plt_locuszoom_scatter <- ggplot(ewas_CRP_sub, aes(x = Probe_position / 1e3, y = Estimate)) +
  geom_point(col = col_CRP, size = point_size_scatter) + 
  scale_x_continuous(breaks = breaks_pretty(n = 4), limits = c(plt_locuszoom$lim_min / 1e3, plt_locuszoom$lim_max / 1e3)) +   
  theme_classic(base_size = font_size_scatter) +
  theme(axis.text = element_text(color = "black")) +
  xlab("Position (KB)") +
  ylab("Direct effect size")

pltG <- plt_locuszoom$plt + plt_locuszoom_scatter + plot_layout(heights = c(1, 6))
pltG[[1]] <- pltG[[1]] + labs(tag = "G")

# H -----------------------------------------------------------------------

thresh <- 5000
ref <- get_transcript_tib(chr = "chr21", 
                          start = 43561179, 
                          stop = 43795442, 
                          thresh = thresh)

ref_sub <- ref %>% 
  filter(Transcript_ID == "NM_207628")

plt_locuszoom <- plot_locus_zoom(ref_sub, thresh = thresh, font_size_label = 3)

ewas_CRP_sub <- ewas_CRP %>% 
  filter(Probe_chr == "21", 
         Probe_position > plt_locuszoom$lim_min, 
         Probe_position < plt_locuszoom$lim_max)

plt_locuszoom_scatter <- ggplot(ewas_CRP_sub, aes(x = Probe_position / 1e3, y = Estimate)) +
  geom_point(col = col_CRP, size = point_size_scatter) + 
  scale_x_continuous(breaks = breaks_pretty(n = 4), limits = c(plt_locuszoom$lim_min / 1e3, plt_locuszoom$lim_max / 1e3)) +   
  theme_classic(base_size = font_size_scatter) +
  theme(axis.text = element_text(color = "black")) +
  xlab("Position (KB)") +
  ylab("Direct effect size")

pltH <- plt_locuszoom$plt + plt_locuszoom_scatter + plot_layout(heights = c(1, 6))
pltH[[1]] <- pltH[[1]] + labs(tag = "H")


# Assemble plot -----------------------------------------------------------

design <- "
  AB
  CD
  EF
  GH
"


plt <- pltA + pltB + pltC + pltD + pltE + pltF + pltG + pltH + plot_layout(design = design) &
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave("./Plots/Revision_plots/Supp_fig3/Supp_fig3.png", plt, width = 172, height = 147, units = "mm", dpi = 500)
ggsave("./Plots/Revision_plots/Supp_fig3/Supp_fig3.pdf", plt, width = 172, height = 147, units = "mm")


