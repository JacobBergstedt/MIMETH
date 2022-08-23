library(tidyverse)
library(FlowSorted.Blood.EPIC)
library(FlowSorted.BloodExtended.EPIC)
library(FlowSorted.Blood.450k)
data(IDOLOptimizedCpGs)
data(FlowSorted.BloodExtended.EPIC)
flowsorted_epic <- readRDS("./Data/RData/Methylation/FlowSorted.Blood.EPIC")
rgset <- readRDS("Data/RData/Methylation/MIMETH.minfi.RGset_969.rds")
ss <- readRDS("./Data/RData/Methylation/Annotation/MIMETH.969_sample_sheet.rds")

# IDOL --------------------------------------------------------------------

cells_IDOL <- estimateCellCounts2(rgset, 
                             compositeCellType = "Blood", 
                             processMethod = "preprocessNoob", 
                             probeSelect = "IDOL", 
                             cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
                             referencePlatform = "IlluminaHumanMethylationEPIC", 
                             referenceset = "flowsorted_epic", 
                             IDOLOptimizedCpGs = IDOLOptimizedCpGs,
                             returnAll = FALSE)

cells_IDOL <- as.data.frame(cells_IDOL$counts) %>% 
  rownames_to_column(var = "SentrixID") %>% 
  left_join(ss[c("SentrixID", "SUBJID")]) %>% 
  dplyr::select(-SentrixID)
saveRDS(cells_IDOL, "./Data/RData/Cells/idol_proportions.rds")

# Houseman ----------------------------------------------------------------

cells_Houseman <- estimateCellCounts(rgset, 
                                     compositeCellType = "Blood", 
                                     processMethod = "preprocessNoob", 
                                     probeSelect = "IDOL", 
                                     cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
                                     referencePlatform = "IlluminaHumanMethylationEPIC",
                                     returnAll = FALSE)


cells_Houseman <- as.data.frame(cells_Houseman) %>% 
  rownames_to_column(var = "SentrixID") %>% 
  left_join(ss[c("SentrixID", "SUBJID")]) %>% 
  dplyr::select(-SentrixID) %>% 
  dplyr::select(SUBJID, everything())
saveRDS(cells_Houseman, "./Data/RData/Cells/Houseman_proportions.rds")



# IDOL extended -----------------------------------------------------------


cells_IDOL_extended <- estimateCellCounts2(rgset, 
                                           compositeCellType = "Blood", 
                                           cellTypes = c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK"),
                                           processMethod = "preprocessNoob",
                                           referencePlatform = "IlluminaHumanMethylationEPIC", 
                                           referenceset = "FlowSorted.BloodExtended.EPIC", 
                                           IDOLOptimizedCpGs = IDOLOptimizedCpGsBloodExtended, 
                                           returnAll = FALSE)

cells_IDOL_extended <- as.data.frame(cells_IDOL_extended$counts) %>% 
  rownames_to_column(var = "SentrixID") %>% 
  left_join(ss[c("SentrixID", "SUBJID")]) %>% 
  dplyr::select(-SentrixID)
saveRDS(cells_IDOL_extended, "./Data/RData/Cells/IDOL_extended_proportions.rds")
