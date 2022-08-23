library(tidyverse)
library(parallel)
source("./Scripts/R_scripts/Libraries/functions_for_annotation.R")
meth <- readRDS("./Data/RData/Methylation/MIMETH.minfi.final.MMatrix.no_outliers.rds")

# One person wants to be removed from study
meth <- meth[meth$SUBJID != 883, ]

loc <- get_anno_meth_location() %>% 
  filter(!Probe_chr %in% c("chrX", "chrY"))

ancestry_ID <- read_tsv("./Data/LabExMI_989x661149_norelat.mds") %>% 
  mutate(SUBJID = as.character(IID)) %>% select(SUBJID)

meth_ID <- inner_join(meth["SUBJID"], ancestry_ID) %>% pull(SUBJID)
meth_autosomes <- meth[meth$SUBJID %in% meth_ID, c("SUBJID", loc$Probe)]

saveRDS(meth_autosomes, "./Data/RData/Methylation/MIMETH.minfi.final.MMatrix.autosomes.no_outliers.rds")
saveRDS(meth_autosomes$SUBJID, "./Data/RData/Methylation/subject_ids.rds")
saveRDS(meth_autosomes[, 1:101], "./Data/RData/Methylation/meth_1e2.rds")
saveRDS(meth_autosomes[, 1:1001], "./Data/RData/Methylation/meth_1e3.rds")
saveRDS(meth_autosomes[, 1:10001], "./Data/RData/Methylation/meth_1e4.rds")
