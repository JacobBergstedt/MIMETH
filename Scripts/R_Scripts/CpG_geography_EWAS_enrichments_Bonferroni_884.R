# Load packages and functions ---------------------------------------------

library(tidyverse)
library(broom)
library(glue)
library(scales)
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
source("./Scripts/R_scripts/Libraries/functions_for_odds_ratios.R")
source("./Scripts/R_scripts/Libraries/functions_for_annotation.R")
source("./Scripts/R_scripts/Libraries/functions_for_plotting.R")

build_tables_for_geography <- function(tib, thresh = 0.05) {
  
  
  tib_pos <- tib %>% 
    filter(!(P_bonferroni < thresh & Estimate < 0)) %>% 
    mutate(Case = (P_bonferroni < thresh & Estimate > 0))
  
  tib_neg <- tib %>% 
    filter(!(P_bonferroni < thresh & Estimate > 0)) %>% 
    mutate(Case = (P_bonferroni < thresh & Estimate < 0))
  
  
  res_pos <- vector(mode = "list", length = length(cell_lines_sub))
  res_neg <- vector(mode = "list", length = length(cell_lines_sub))
  names(res_pos) <- cell_lines_sub_label
  names(res_neg) <- cell_lines_sub_label
  
  pos_is_case <- tib_pos$Case
  neg_is_case <- tib_neg$Case
  regions_pos <- tib_pos$Geography
  regions_neg <- tib_neg$Geography

  
  res_neg <- map(set_names(geography_regions, geography_regions), ~ regions_neg == .) %>% 
    map(function(is_region) xtabs(~ is_region + neg_is_case))
  
  res_pos <- map(set_names(geography_regions, geography_regions), ~ regions_pos == .) %>% 
    map(function(is_region) xtabs(~ is_region + pos_is_case))
  
  list(neg = res_neg, pos = res_pos)
  
}

screen <- function(x) {
  
  x <- as.data.frame(x)
  target_hits <- x %>% filter(as.logical(x[[2]])) %>% pull(Freq)
  all(target_hits > 0) & any(target_hits > 1)
  
}

do_fisher_tests <- function(table_list) {
  
  unlist2(table_list) %>%
    keep(screen) %>% 
    map(fisher.test) %>% 
    map_dfr(tidy, .id = "ID") %>%
    separate(col = ID, into = c("Direction", "Geographic_region"), sep = "__") %>% 
    dplyr::rename(Estimate = estimate, Low = conf.low, High = conf.high, P_value = p.value)
  
}

run_geography_enrichment_pipeline <- function(tib) {
  
  build_tables_for_geography(tib) %>% 
    do_fisher_tests()
}



# Load crucial data -------------------------------------------------------

meth_anno_geography <- get_anno_meth_geography()
meth_location <- get_anno_meth_location(fix = FALSE)

geography_regions <- levels(meth_anno_geography$Geography)
names(geography_regions) <- levels(meth_anno_geography$Geography_labels)


# Load EWAS data ----------------------------------------------------------

ewas_age <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Age.rds") %>% 
  left_join(meth_location) %>% 
  left_join(meth_anno_geography) %>% 
  mutate(Probe_gene = simplify_probe_genes(Probe_gene))

ewas_sex <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Sex.rds") %>% 
  left_join(meth_location) %>% 
  left_join(meth_anno_geography) %>% 
  mutate(Probe_gene = simplify_probe_genes(Probe_gene))

ewas_CMV <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/CMV_serostatus.rds") %>% 
  left_join(meth_location) %>% 
  left_join(meth_anno_geography) %>% 
  mutate(Probe_gene = simplify_probe_genes(Probe_gene))

ewas_smoking <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Smoking_status.rds") %>% 
  left_join(meth_location) %>% 
  left_join(meth_anno_geography) %>% 
  filter(Levels != "Ex_Smoker") %>% 
  mutate(Probe_gene = simplify_probe_genes(Probe_gene))

thresh <- 0.1
age_pos_probes <- ewas_age %>% filter(P_bonferroni < thresh, Estimate > 0) %>% pull(Probe)
sex_pos_probes <- ewas_sex %>% filter(P_bonferroni < thresh, Estimate > 0) %>% pull(Probe)
CMV_pos_probes <- ewas_CMV %>% filter(P_bonferroni < thresh, Estimate > 0) %>% pull(Probe)
smoking_pos_probes <- ewas_smoking %>% filter(P_bonferroni < thresh, Estimate > 0) %>% pull(Probe)

age_neg_probes <- ewas_age %>% filter(P_bonferroni < thresh, Estimate < 0) %>% pull(Probe)
sex_neg_probes <- ewas_sex %>% filter(P_bonferroni < thresh, Estimate < 0) %>% pull(Probe)
CMV_neg_probes <- ewas_CMV %>% filter(P_bonferroni < thresh, Estimate < 0) %>% pull(Probe)
smoking_neg_probes <- ewas_smoking %>% filter(P_bonferroni < thresh, Estimate < 0) %>% pull(Probe)


# Mediation probes --------------------------------------------------------

ewas_med_age <- readRDS("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/Age.rds") %>% 
  filter(Effect == "Mediation") %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  left_join(meth_location) %>% 
  left_join(meth_anno_geography)

ewas_med_sex <- readRDS("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/Sex.rds") %>% 
  filter(Effect == "Mediation") %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  left_join(meth_location) %>% 
  left_join(meth_anno_geography)

ewas_med_smoking <- readRDS("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/Smoking.rds") %>% 
  filter(Effect == "Mediation") %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  left_join(meth_location) %>% 
  left_join(meth_anno_geography)

ewas_med_CMV <- readRDS("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/CMV.rds") %>% 
  filter(Effect == "Mediation") %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  left_join(meth_location) %>% 
  left_join(meth_anno_geography)

ewas_med_age_sign <- ewas_med_age %>% filter(P_bonferroni < thresh)
ewas_med_sex_sign <- ewas_med_sex %>% filter(P_bonferroni < thresh)
ewas_med_smoking_sign <- ewas_med_smoking %>% filter(P_bonferroni < thresh)
ewas_med_CMV_sign <- ewas_med_CMV %>% filter(P_bonferroni < thresh)

age_med_pos_probes <- ewas_med_age_sign %>% filter(Estimate > 0) %>% pull(Probe)
age_med_neg_probes <- ewas_med_age_sign %>% filter(Estimate < 0) %>% pull(Probe)
CMV_med_pos_probes <- ewas_med_CMV_sign %>% filter(Estimate > 0) %>% pull(Probe)
CMV_med_neg_probes <- ewas_med_CMV_sign %>% filter(Estimate < 0) %>% pull(Probe)

# Run geography enrichments -----------------------------------------------

res <- list(Mediated_Smoking = ewas_med_smoking, 
            Direct_Smoking = ewas_smoking, 
            Mediated_Age = ewas_med_age, 
            Direct_Age = ewas_age, 
            Mediated_CMV = ewas_med_CMV, 
            Direct_CMV = ewas_CMV, 
            Mediated_Sex = ewas_med_sex, 
            Direct_Sex = ewas_sex) %>% 
  map_dfr(run_geography_enrichment_pipeline, .id = "Type") %>% 
  separate(Type, into = c("Type", "Variable"), sep = "_") %>% 
  mutate(Geographic_region = factor(Geographic_region, geography_regions),
         Geographic_region_labels = fct_recode(Geographic_region, !!!geography_regions))

saveRDS(res, "./Data/RData/Results/EWAS_geographic_region_enrichments_bonferroni_884.rds")


# Plot enrichments --------------------------------------------------------

res_plt <- res %>%
  dplyr::rename(Term = Geographic_region_labels)

for (variable in unique(res_plt$Variable)) {
  for (type in unique(res_plt$Type)) {
    plt <- res_plt %>%
      filter(Variable == variable, Type == type) %>%
      plot_odds_ratio_split(paste(variable, type, sep = ","),
                            point_size = point_size_volcano,
                            bar_size = odds_ratio_errorbar_size,
                            colors = blues[c(3, 6)],
                            base_size = font_size_odds_ratio,
                            legend_position = c(0.88, 0.94))

    ggsave(paste0("./Plots/Odds_ratios_bonferroni_884/Geography/", type, "/", variable, ".pdf"), plt, width = 5, height = 5, units = "cm")
  }
}

