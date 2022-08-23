library(tidyverse)
library(magrittr)
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")


logit <- function(x) log2(x / (1 - x))

transform_to_open_interval <- function(x){
  n <- length(x)
  (x * (n - 1) + 0.5) / n  
}

get_20_PCs <- function(x) {
  id <- x$SUBJID
  x <- x %>% 
    select(-SUBJID) %>% 
    as.matrix()
  x <- prcomp(x, center = TRUE, scale. = TRUE)$x[, 1:20]
  bind_cols(SUBJID = id, as_tibble(x))
}

# Load methylation data ---------------------------------------------------
load("Data/RData/globals.RData")
meth_ID <- get_m_values_sample() %>% select(SUBJID)

# Ancestry principal components -------------------------------------------
ancestry_pcs <- read_tsv("./Data/LabExMI_989x661149_norelat.mds") %>%
  mutate(SUBJID = as.character(IID)) %>%
  select(-FID, -IID) %>%
  rename_at(vars(-SUBJID), ~ paste0("Ancestry_", .))

# Organize proportion data ------------------------------------------------
props <- readRDS("./Data/RData/Cells/proportions_for_methylation_adjustment.rds")
sum_props <- rowSums(props[-1])
props[-1] <- props[-1] / sum_props

# props_PCs <- get_20_PCs(props)
# colnames(props_PCs) <- c("SUBJID", paste0("props_PC", 1:20))


# Counts ------------------------------------------------------------------
# counts <- readRDS("./Data/RData/Cells/counts_imputed.rds")
# counts_PCs <- get_20_PCs(counts)
# names(counts_PCs) <- c("SUBJID", paste0("counts_PC", 1:20))

# IDOL cells --------------------------------------------------------------
idol_cells <- readRDS("./Data/RData/Cells/idol_proportions.rds") %>% 
  rename_at(vars(-SUBJID), ~ paste0("IDOL_", .))


# IDOL extended --------------------------------------------------------------
idol_cells_extended <- readRDS("./Data/RData/Cells/IDOL_extended_proportions.rds") %>% 
  rename_at(vars(-SUBJID), ~ paste0("IDOL_extended_", .))


# Houseman cells --------------------------------------------------------------
houseman_cells <- readRDS("./Data/RData/Cells/Houseman_proportions.rds") %>% 
  rename_at(vars(-SUBJID), ~ paste0("Houseman_", .))

# eCRF --------------------------------------------------------------------
ecrf <- readRDS("./Data/RData/Environment/ecrf_curated_2019.rds") %>% 
  select(-Ancestry_PC1, -Ancestry_PC2)

ecrf$MCHC[is.na(ecrf$MCHC)] <- mean(ecrf$MCHC, na.rm = TRUE)

# Methylation batch variables ---------------------------------------------
meth_batch <- read_tsv("./Data/mimeth.DNA.extraction.txt") %>%
  mutate(SUBJID = as.character(SUBJID)) %>% 
  select(-Sample_ID_1.5mltube) %>% 
  dplyr::rename(Batch_sample_well = Sample_Well, 
         Batch_sample_plate = Sample_Plate, 
         Batch_sentrix_ID = Sentrix_ID,
         Batch_sentrix_position = Sentrix_Position, 
         Batch_batch = Batch, 
         Batch_box = Box, 
         Batch_original_DNA_concentation = Original_Conc) %>% 
  mutate(Batch_sentrix_position_linear = as.numeric(factor(Batch_sentrix_position)) - 1)

# Merge -------------------------------------------------------------------
covariates <- props %>% 
  inner_join(ecrf) %>%  
  inner_join(idol_cells) %>% 
  inner_join(houseman_cells) %>% 
  inner_join(ancestry_pcs) %>%
  inner_join(meth_batch) %>% 
  inner_join(idol_cells_extended)

covariates_all_samples <- covariates
covariates <- inner_join(meth_ID, covariates)

# Write to file -----------------------------------------------------------
write(names(select(ecrf, -SUBJID, -Day_of_sampling, -Date_of_sampling, -Height)),
      "./Data/covariate_list.txt")

saveRDS(covariates, "./Data/RData/Environment/covariates_884.rds")

