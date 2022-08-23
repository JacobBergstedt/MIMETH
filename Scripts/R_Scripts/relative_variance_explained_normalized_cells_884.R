
library(tidyverse)
library(glue)
library(parallel)
library(MatrixEQTL)
library(tidyverse)
library(snpStats)
library(splines)
library(dtplyr)
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
source("./Scripts/R_scripts/Libraries/functions_for_variance_explained2.R")
source("./Scripts/R_scripts/Libraries/functions_for_meQTL_mapping2.R")
source("./Scripts/MiscFunctions.R")
options(stringsAsFactors = FALSE)
options(width = 200)

n_cores <- 5
rep <- 1:5

meth_locs <- get_anno_meth_location()
cell_list_15_cells <- get_15_props()
fm_cells <- " ~ " %>% str_add_str(paste0(cell_list_15_cells, collapse = " + "))

# This analysis is per chromosome. The following sets up lists containing needed information for each chromosome -----------------------

read_meth_chr <- function(chr) readRDS(glue::glue("./Data/Chunk_data/Methylation/Beta/Per_chromosome/meth_chromosome_{chr}.rds"))
read_snp_mat_chr <- function(chr) readRDS(glue::glue("./Data/Chunk_data/Genotype/Chromosome_{chr}/LabExMI_imputation_958x5699237_chr{chr}_snp_matrix.rds"))

mf <- get_covs_884() %>% dplyr::select(SUBJID, contains("IDOL"), contains("Houseman"), all_of(cell_list_15_cells), Sex, Smoking_status, Log_CRP_levels, Age, CMV_serostatus) %>% 
  as.data.frame()

mf_panel5 <- get_panel5_cells() %>% 
  rename_with(~ paste0("Panel5_", str_sub(., end = -8)), .cols = !SUBJID)

row_sums <- select(mf_panel5, !SUBJID) %>% 
  rowSums()

mf_panel5 <- mf_panel5 %>% mutate(across(!SUBJID, ~ . / row_sums))
mf <- mf %>% inner_join(mf_panel5)
rownames(mf) <- mf$SUBJID

cl <- makeCluster(n_cores)
clusterExport(cl = cl, "meth_locs")

meth_list <- parLapply(cl, 1:22, read_meth_chr)
snp_list <- parLapply(cl, 1:22, read_snp_mat_chr)
meth_loc_list <- parLapply(cl, 1:22, setup_meth_loc_chr)
snp_loc_list <- parLapply(cl, 1:22, setup_snp_loc_chr)

stopCluster(cl)

meth_list <- map(meth_list, ~ .[match(rownames(mf), .$SUBJID), ])
snp_list <- map(snp_list, ~ .[rownames(mf),])

# Start computations ------------------------------------------------------

print(n_cores)
print("Start loop!")

gc()

for (i in rep) {
  
  res <- compute_relative_importance_cv(meth_list,
                                        snp_list,
                                        meth_loc_list = meth_loc_list,
                                        snp_loc_list = snp_loc_list,
                                        mf = mf,
                                        fm_cells = fm_cells,
                                        n_cores = n_cores)
  
  saveRDS(res, paste0("./Data/Chunk_data/Results/Variance_explained/Normalized_cells_884/Beta/res_rep_", i, ".rds"))
  
}


