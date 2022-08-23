
# Initialize --------------------------------------------------------------
library(tidyverse)
library(stabs)
library(glmnet)
library(parallel)
source("./Scripts/R_scripts/Libraries/functions_for_CMV_prediction.R")

# Load data ---------------------------------------------------------------

meth <- readRDS("./Data/RData/Methylation/MIMETH.minfi.MMatrix.noob_969.ComBat2.rds")
ss <- readRDS("./Data/RData/Methylation/Annotation/MIMETH.969_sample_sheet.rds")
meth <- meth[, ss$SentrixID]
colnames(meth) <- ss$SUBJID
meth <- t(meth)
covs <- data.frame(SUBJID = rownames(meth)) %>% 
  left_join(readRDS("./Data/RData/Environment/covariates_all_samples.rds") )

y <- covs$CMV_serostatus

# Define globals ----------------------------------------------------------

q <- 50
tol <- 2
alpha <- 0.95
n_rep <- 4
fold_list <- replicate(n = n_rep, sample.int(10, length(y), replace = TRUE), simplify = FALSE)

# Run cross validation ----------------------------------------------------

selection_runs <- map_dfr(fold_list, cv_stability_selection, meth, y, alpha = alpha, q = q, tol = tol)

saveRDS(selection_runs, "./Data/RData/CMV_estimation_accuracy_stabsel.rds")

