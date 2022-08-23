library(tidyverse)
library(parallel)
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
source("./Scripts/R_scripts/Libraries/functions_for_data_splitting.R")
meth_ID <- get_meth_ids()
snps <- get_snp_matrix(NULL)[meth_ID, ]
snp_locations <- get_snp_location()
meth_locations <- get_meth_location()
invisible(map(1:22, save_genotypes_per_chr, meth_locations = meth_locations, snp_locations = snp_locations))