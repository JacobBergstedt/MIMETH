# Load packages and functions ---------------------------------------------

library(tidyverse)
library(broom)
library(glue)
library(scales)
library(parallel)
library(missMethyl)
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
source("./Scripts/R_scripts/Libraries/functions_for_odds_ratios.R")
source("./Scripts/R_scripts/Libraries/functions_for_annotation.R")
source("./Scripts/R_scripts/Libraries/functions_for_plotting.R")
source("./Scripts/R_scripts/Libraries/functions_for_enrichments.R")
rename <- dplyr::rename

# Load crucial data -------------------------------------------------------

cell_list <- scan("./Data/prop_controls_16.txt", what = character())
cell_lines <- get_cell_list_roadmap()
cell_lines_sub <- cell_lines
cell_lines_sub_label <- gsub("_", " ", gsub("_cells", "", cell_lines_sub))
cell_lines_sub <- paste0("Chromatin_states_", cell_lines_sub)
states <- names(get_anno_roadmap_translation())

# beta_values <- get_beta_values()
# m_values <- get_m_values()
covs <- get_covs_884()

meth_anno_roadmap <- get_anno_meth_roadmap_all_cells()
# roadmap_desc <- read_tsv("./Data/Roadmap_15_states/roadmap_descriptions.tsv")
# meth_anno <- get_meth_annotation()
meth_anno_geography <- get_anno_meth_geography()
meth_location <- get_anno_meth_location(fix = FALSE)
# ewas_others <- readRDS("./Data/RData/Results/EWAS/M_values/Environment/Random_effect/limited_evidence_sign.rds")

# Load EWAS data ----------------------------------------------------------

ewas_age <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Age.rds") %>% 
  left_join(meth_anno_roadmap) %>% 
  left_join(meth_location) %>%
  mutate(Probe_gene = simplify_probe_genes(Probe_gene))

ewas_sex <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Sex.rds") %>% 
  left_join(meth_anno_roadmap) %>% 
  left_join(meth_location) %>%
  mutate(Probe_gene = simplify_probe_genes(Probe_gene))

ewas_CMV <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/CMV_serostatus.rds") %>% 
  left_join(meth_anno_roadmap) %>% 
  left_join(meth_location) %>%
  mutate(Probe_gene = simplify_probe_genes(Probe_gene))

ewas_smoking <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Smoking_status.rds") %>% 
  left_join(meth_anno_roadmap) %>% 
  left_join(meth_location) %>%
  filter(Levels != "Ex_Smoker") %>% 
  mutate(Probe_gene = simplify_probe_genes(Probe_gene))


# Mediation probes --------------------------------------------------------

ewas_med_age <- readRDS("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/Age.rds") %>% 
  filter(Effect == "Mediation") %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  left_join(meth_anno_roadmap) %>% 
  left_join(meth_location)

ewas_med_sex <- readRDS("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/Sex.rds") %>% 
  filter(Effect == "Mediation") %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  left_join(meth_anno_roadmap) %>% 
  left_join(meth_location)

ewas_med_smoking <- readRDS("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/Smoking.rds") %>% 
  filter(Effect == "Mediation") %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  left_join(meth_anno_roadmap) %>% 
  left_join(meth_location)

ewas_med_CMV <- readRDS("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/CMV.rds") %>% 
  filter(Effect == "Mediation") %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  left_join(meth_anno_roadmap) %>% 
  left_join(meth_location)


ewas_list <- list(Mediated_Smoking = ewas_med_smoking, 
                  Direct_Smoking = ewas_smoking, 
                  Mediated_Age = ewas_med_age, 
                  Direct_Age = ewas_age, 
                  Mediated_CMV = ewas_med_CMV, 
                  Direct_CMV = ewas_CMV, 
                  Mediated_Sex = ewas_med_sex, 
                  Direct_Sex = ewas_sex)


# Cross-tab Bonferroni ----------------------------------------------------

res <- ewas_list %>% 
  mclapply(run_EWAS_enrichment_pipeline_bonferroni, mc.cores = 12) %>% 
  bind_rows(.id = "Type") %>% 
  separate(Type, into = c("Type", "Variable"), sep = "_") %>% 
  mutate(Term = factor(State, states))

saveRDS(res, "./Data/RData/Results/EWAS_chromatin_state_enrichments_bonferroni_884.rds")




# Outside of CGI ----------------------------------------------------------

meth_anno_geo <- get_anno_meth_geography()
meth_anno_roadmap_non_CGI <- meth_anno_roadmap %>% 
  left_join(meth_anno_geo) %>% 
  filter(Geography != "Island")

ewas_age_non_CGI <- left_join(meth_anno_roadmap_non_CGI, ewas_age)

res_direct_age_non_cgi <- run_EWAS_enrichment_pipeline_bonferroni(ewas_age_non_CGI)
saveRDS(res_direct_age_non_cgi, "./Data/RData/Results/EWAS_chromatin_state_enrichments_bonferroni_884_age_non_CGI.rds")


# GO enrichment -----------------------------------------------------------


res_go_EWAS <- mclapply(ewas_list, gometh_enrichment_EWAS, mc.cores = 8) %>% 
  bind_rows(.id = "EWAS") %>%
  filter(ONTOLOGY == "BP") %>% 
  separate(col = Terms, into = c("Direction", "Thresh")) %>% 
  separate(col = EWAS, into = c("Type", "Variable")) %>% 
  group_by(Type, Variable, Direction, Thresh) %>% 
  arrange(FDR, .by_group = TRUE) %>% 
  ungroup()
  
saveRDS(res_go_EWAS, "./Data/RData/Results/EWAS_go_enrichments.rds")


keep <- ewas_list[[1]]$Chromatin_states_Mononuclear_cells %in% c("TSS", "Fl. TSS")
ewas_list_TSS <- map(ewas_list, ~ .[keep, ])

res_go_EWAS_TSS <- mclapply(ewas_list_TSS[-6], gometh_enrichment_EWAS, mc.cores = 6) %>% 
  bind_rows(.id = "EWAS") %>%
  filter(ONTOLOGY == "BP") %>% 
  separate(col = Terms, into = c("Direction", "Thresh")) %>% 
  separate(col = EWAS, into = c("Type", "Variable")) %>% 
  group_by(Type, Variable, Direction, Thresh) %>% 
  arrange(FDR, .by_group = TRUE) %>% 
  ungroup()

saveRDS(res_go_EWAS_TSS, "./Data/RData/Results/EWAS_go_enrichments_TSS.rds")

# Facet plots -------------------------------------------------------------

plot_state_enrichments_all_cells <- function(tib, trans = FALSE, title) {
  plt <- tib %>% 
    ggplot(aes(y = Estimate, 
               ymin = Low, 
               ymax = High, 
               x = fct_rev(Term), 
               color = Direction)) +
    geom_errorbar(size = odds_ratio_errorbar_size, width = odds_ratio_errorbar_size) +
    geom_hline(yintercept = 1, size = line_size_blacks) +
    ylab(NULL) +
    xlab(NULL) +
    geom_point(size = odds_ratio_point_size) +
    ggtitle(title) +
    scale_color_manual(values = blues[c(3, 6)])
  
  if (trans) plt <- plt + scale_y_continuous(trans = "log10")
  
  plt <- plt +
    coord_flip() +
    theme_bw(base_size = font_size_odds_ratio) +
    theme(axis.text = element_text(color = "black"),
          legend.position = c(0.88, 0.94),
          plot.title = element_text(hjust = 0.5),
          panel.grid.minor.x = element_blank(),
          legend.key.size = unit(0.1, "lines"),
          strip.background = element_rect(fill = "white")) +
    facet_wrap(vars(Cell), nrow = 1)
}



# AGE ---------------------------------------------------------------------

plt <- res %>% 
  filter(Type == "Direct", Variable == "Age") %>%
  plot_state_enrichments_all_cells(title = "Direct aging effects")

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Age/Direct_all_cell_lines.pdf", plt, width = 10, height = 5, units = "cm")

plt <- res %>% 
  filter(Type == "Mediated", Variable == "Age") %>% 
  plot_state_enrichments_all_cells(title = "Mediated aging effects")

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Age/Mediated_all_cell_lines.pdf", plt, width = 10, height = 5, units = "cm")

plt <- res %>% 
  filter(Type == "Direct", Variable == "Age") %>% 
  plot_state_enrichments_all_cells(trans = TRUE, title = "Direct aging effects")
ggsave("./Plots/Odds_ratios_bonferroni/Roadmap/Age/Direct_all_cell_lines_transformed.pdf", plt, width = 10, height = 5, units = "cm")

plt <- res %>% 
  filter(Type == "Mediated", Variable == "Age") %>% 
  plot_state_enrichments_all_cells(trans = TRUE, title = "Mediated aging effects")
ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Age/Mediated_all_cell_lines_transformed.pdf", plt, width = 10, height = 5, units = "cm")


# CMV ---------------------------------------------------------------------


plt <- res %>% 
  filter(Type == "Direct", Variable == "CMV") %>%
  plot_state_enrichments_all_cells(title = "CMV, direct effects")

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/CMV/Direct_all_cell_lines.pdf", plt, width = 10, height = 5, units = "cm")

plt <- res %>% 
  filter(Type == "Mediated", Variable == "CMV") %>% 
  plot_state_enrichments_all_cells(title = "CMV, mediated effects")

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/CMV/Mediated_all_cell_lines.pdf", plt, width = 10, height = 5, units = "cm")

plt <- res %>% 
  filter(Type == "Direct", Variable == "CMV") %>% 
  plot_state_enrichments_all_cells(trans = TRUE, title = "CMV, direct effects")
ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/CMV/Direct_all_cell_lines_transformed.pdf", plt, width = 10, height = 5, units = "cm")

plt <- res %>% 
  filter(Type == "Mediated", Variable == "CMV") %>% 
  plot_state_enrichments_all_cells(trans = TRUE, title = "CMV, mediated effects")
ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/CMV/Mediated_all_cell_lines_transformed.pdf", plt, width = 10, height = 5, units = "cm")


# Smoking -----------------------------------------------------------------


plt <- res %>% 
  filter(Type == "Direct", Variable == "Smoking") %>%
  plot_state_enrichments_all_cells(title = "Smoking, direct effects")

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Smoking/Direct_all_cell_lines.pdf", plt, width = 10, height = 5, units = "cm")

plt <- res %>% 
  filter(Type == "Mediated", Variable == "Smoking") %>% 
  plot_state_enrichments_all_cells(title = "Smoking, mediated effects")

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Smoking/Mediated_all_cell_lines.pdf", plt, width = 10, height = 5, units = "cm")

plt <- res %>% 
  filter(Type == "Direct", Variable == "Smoking") %>% 
  plot_state_enrichments_all_cells(trans = TRUE, title = "Smoking, direct effects")
ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Smoking/Direct_all_cell_lines_transformed.pdf", plt, width = 10, height = 5, units = "cm")

plt <- res %>% 
  filter(Type == "Mediated", Variable == "Smoking") %>% 
  plot_state_enrichments_all_cells(trans = TRUE, title = "Smoking, mediated effects")
ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Smoking/Mediated_all_cell_lines_transformed.pdf", plt, width = 10, height = 5, units = "cm")


# Smoking -----------------------------------------------------------------


plt <- res %>% 
  filter(Type == "Direct", Variable == "Sex") %>%
  plot_state_enrichments_all_cells(title = "Sex, direct effects")

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Sex/Direct_all_cell_lines.pdf", plt, width = 10, height = 5, units = "cm")

plt <- res %>% 
  filter(Type == "Mediated", Variable == "Sex") %>% 
  plot_state_enrichments_all_cells(title = "Sex, mediated effects")

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Sex/Mediated_all_cell_lines.pdf", plt, width = 10, height = 5, units = "cm")

plt <- res %>% 
  filter(Type == "Direct", Variable == "Sex") %>% 
  plot_state_enrichments_all_cells(trans = TRUE, title = "Sex, direct effects")
ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Sex/Direct_all_cell_lines_transformed.pdf", plt, width = 10, height = 5, units = "cm")

plt <- res %>% 
  filter(Type == "Mediated", Variable == "Sex") %>% 
  plot_state_enrichments_all_cells(trans = TRUE, title = "Sex, mediated effects")
ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Sex/Mediated_all_cell_lines_transformed.pdf", plt, width = 10, height = 5, units = "cm")



# Plot enrichments --------------------------------------------------------

for (cell in unique(res$Cell)) {
  for (variable in unique(res$Variable)) {
    for (type in unique(res$Type)) {
      plt <- res %>% 
        filter(Variable == variable, Type == type, Cell == cell) %>%
        plot_odds_ratio_split(paste(variable, type, cell, sep = ","),
                              point_size = point_size_volcano,
                              bar_size = odds_ratio_errorbar_size,
                              colors = blues[c(3, 6)],
                              base_size = font_size_odds_ratio,
                              legend_position = c(0.88, 0.94))
      
      ggsave(paste0("./Plots/Odds_ratios_bonferroni_884/Roadmap/", variable, "/", type, "_", cell, ".pdf"), plt, width = 5, height = 5, units = "cm")      
    }
  }
}

