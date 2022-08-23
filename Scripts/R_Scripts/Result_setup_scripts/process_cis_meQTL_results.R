library(tidyverse)
library(parallel)
source("./Scripts/R_scripts/Libraries/functions_for_collections.R")
source("./Scripts/R_scripts/Libraries/functions_for_meQTL_mapping2.R")
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")

# BH adjustment
cis <- readRDS("./Data/RData/Results/MeQTL/cis_m_values_884.rds")
cis_adjusted <- two_stage_adjust(cis, 
                                 top_level = "Probe", 
                                 level = 0.05, 
                                 bottom_level = "SNP", 
                                 Bonferroni_bottom_adjustment = TRUE)

cis_adjusted_probes <- distinct(cis_adjusted, Probe, .keep_all = TRUE)
saveRDS(cis_adjusted_probes, "./Data/RData/Results/MeQTL/cis_adjusted_probes_884.rds")



# Bonferroni
cis <- readRDS("./Data/RData/Results/MeQTL/cis_m_values_884.rds")
cis$P_bonferroni <- p.adjust(cis$Pvalue)

# cis_adjusted <- two_stage_adjust(cis, 
#                                  top_level = "Probe", 
#                                  level = 0.05, 
#                                  bottom_level = "SNP", 
#                                  Bonferroni_bottom_adjustment = TRUE, 
#                                  manual_top_adjustment_factor = 643546)

cis_adjusted_probes <- distinct(cis_adjusted, Probe, .keep_all = TRUE)
saveRDS(cis_adjusted_probes, "./Data/RData/Results/MeQTL/cis_adjusted_probes_884_bonferroni.rds")