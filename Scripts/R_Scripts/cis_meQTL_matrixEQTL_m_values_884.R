
library(MatrixEQTL)
library(tidyverse)
library(snpStats)
library(parallel)
library(splines)
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
source("./Scripts/R_scripts/Libraries/functions_for_meQTL_mapping2.R")

setup_meth <- function(chr, type, covs) {
  
  meth_chr <- readRDS(glue::glue("./Data/Chunk_data/Methylation/{type}/Per_chromosome/meth_chromosome_{chr}.rds"))
  rownames(meth_chr) <- meth_chr$SUBJID
  meth_chr$SUBJID <- NULL
  meth_chr <- as.matrix(meth_chr[covs$SUBJID, ])
  setup_methylation_for_meQTL_mapping(meth_chr)
}

setup_genotypes <- function(chr, covs) {
  
  snps <- readRDS(glue::glue("./Data/Chunk_data/Genotype/Chromosome_{chr}/LabExMI_imputation_958x5699237_chr{chr}_snp_matrix.rds"))
  snps <- snps[covs$SUBJID, ]
  setup_genotypes_for_meQTL_mapping(snps)
  
}

# This analysis is per chromosome. The following sets up lists containing needed information for each chromosome -----------------------

cell_list_15_cells <- get_15_props()
meth_locs <- get_anno_meth_location()

cl <- makeCluster(11)
clusterExport(cl = cl, "meth_locs")
clusterEvalQ(cl = cl, {
  source("./Scripts/R_scripts/Libraries/functions_for_meQTL_mapping2.R")
  library(MatrixEQTL)
})

covs <- get_covs_884()
meth_list <- parLapply(cl, 1:22, setup_meth, type = "M_values", covs = covs)
snp_list <- parLapply(cl, 1:22, setup_genotypes, covs = covs)
meth_loc_list <- parLapply(cl, 1:22, setup_meth_loc_chr)
snp_loc_list <- parLapply(cl, 1:22, setup_snp_loc_chr)

stopCluster(cl)

# Sets up covariate matrix for matrixEQTL ---------------------------------

fm_control <- get_control_fm_15_cells()
mf <- as.data.frame(covs)
rownames(mf) <- mf$SUBJID
mf <- model.matrix(fm_control, mf)[, -1]
meqtl_covariates <- setup_covariates_for_meQTL_mapping(mf)

# Run loop for each chromosome --------------------------------------------
res <- mapply(run_local_meqtl_for_chr, meth_list, snp_list, meth_loc_list, snp_loc_list, MoreArgs = list(mm = meqtl_covariates), SIMPLIFY = FALSE)
res <- bind_rows(res)

saveRDS(res, "./Data/RData/Results/MeQTL/cis_m_values_884.rds")


