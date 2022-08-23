library(lumi)
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")

beta_values <- readRDS("./Data/RData/Methylation/MIMETH.minfi.final.betaMatrix.no_outliers.rds")
beta_values <- beta_values[beta_values$SUBJID != 883, ]
loc <- get_meth_location_annotation() %>% 
  filter(!Chromosome %in% c("chrX", "chrY"))

ancestry_IDs <- read_tsv("./Data/LabExMI_989x661149_norelat.mds") %>% 
  mutate(SUBJID = as.character(IID)) %>% select(SUBJID)

meth_ID <- beta_values["SUBJID"] %>% inner_join(ancestry_IDs) %>% pull(SUBJID)
beta_autosomes <- beta_values[beta_values$SUBJID %in% meth_ID, c("SUBJID", loc$Probe)]

saveRDS(beta_autosomes, "./Data/RData/Methylation/MIMETH.minfi.final.betaMatrix.autosomes.no_outliers.rds")
saveRDS(beta_autosomes[1:100], "./Data/RData/Methylation/meth_beta_values_small.rds")
