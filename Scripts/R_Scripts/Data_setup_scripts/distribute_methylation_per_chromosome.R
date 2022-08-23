library(tidyverse)
library(parallel)
source("./Scripts/R_scripts/Libraries/functions_for_data_splitting.R")

m_values <- readRDS("./Data/RData/Methylation/MIMETH.minfi.final.MMatrix.autosomes.no_outliers.rds")
beta_values <- readRDS("./Data/RData/Methylation/MIMETH.minfi.final.betaMatrix.autosomes.no_outliers.rds")
meth_locations <- get_anno_meth_location()

invisible(mclapply(1:22, save_m_values_per_chr, m_values = m_values, meth_locations = meth_locations, mc.cores = 12))
invisible(mclapply(1:22, save_beta_values_per_chr, beta_values = beta_values, meth_locations = meth_locations, mc.cores = 12))



