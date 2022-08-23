library(tidyverse)
library(matrixStats)
library(irlba)
meth <- readRDS("./Data/RData/Methylation/MIMETH.minfi.final.MMatrix.rds")
meth_id <- meth$SUBJID
meth$SUBJID <- NULL
meth <- as.matrix(meth)
meth_PCA <- prcomp_irlba(meth, n = 50, center = TRUE, scale. = TRUE)$x %>% 
  as_tibble()
meth_PCA$SUBJID <- meth_id
saveRDS(meth_PCA, "./Data/RData/Methylation/meth_PCA_50_components.rds")


