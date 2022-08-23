

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
n_rep <- 5
fold_list <- replicate(n = n_rep, sample.int(10, length(y), replace = TRUE), simplify = FALSE)


# Run cross validation ----------------------------------------------------

alpha_seq <- c(0.5, 0.8, 0.95, 0.99)

res <- cv_alpha_logreg(meth, 
                       y, 
                       alpha_seq = alpha_seq, 
                       fold_list = fold_list,
                       n_lambda = 1e2,
                       n_cores = 4)

res_tib <- map_dfr(flatten(flatten(res)), 
                   ~ tibble(Alpha = .$alpha, Lambda_nr = .$lambda_nr, Lambda = .$lambda, CR = .$CR, Df = .$pred_nonzero), 
                   .id = "Resample") %>% 
  group_by(Lambda_nr, Alpha) %>% 
  summarize(Df = Df[1], Lambda = Lambda[1], CR_mean = mean(CR), CR_sd = sd(CR)) %>% 
  ungroup() %>% 
  arrange(desc(CR_mean))

param <- res_tib %>% filter(CR_mean == max(CR_mean))
lambda_max <- param$Lambda
alpha_max <- param$Alpha
m_final <- glmnet(x = meth, y = y, family = "binomial", alpha = alpha_max)
coefs <- predict(m_final, s = lambda_max, type = "coefficients")

saveRDS(res, "./Data/RData/CMV_estimation_list_elastic_net.rds")
saveRDS(res_tib, "./Data/RData/CMV_estimation_tib_elastic_net.rds")
saveRDS(coefs, "./Data/RData/CMV_model_coefs.rds")
