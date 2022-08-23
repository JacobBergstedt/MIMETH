library(tidyverse)
library(glue)
library(vroom)
library(GenomicRanges)
source("./Scripts/R_scripts/Libraries/functions_for_annotation.R")

map_labex <- readRDS("./Data/RData/Genotypes/LabExMI_imputation_1000x5699237_annotated_map_with_ancestral.rds") %>% 
  mutate(SNP_chr = paste0("chr", SNP_chr))

map_labex_roadmap <- map_labex %>%
  add_15_state_annotation(chr_col = "SNP_chr", position_col = "SNP_position", genomic_feature = "SNP_ID")


map_labex <- inner_join(map_labex, map_labex_roadmap)
saveRDS(map_labex, "./Data/RData/Genotypes/LabExMI_imputation_1000x5699237_annotated_map.rds")

