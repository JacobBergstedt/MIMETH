#!/local/gensoft2/exe/R/3.6.0/scripts/Rscript

library(tidyverse)
library(parallel)
library(gamlss)
library(splines)
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")

fit_age_dispersion <- function(y, snp_mat, fm, mf) {
  mf$y <- y
  fm <- fm_add_snps(fm, snp_mat)
  if (!is.null(snp_mat)) mf <- cbind(mf, snp_mat)
  mf <- na.omit(mf)
  null <- gamlss(fm, 
                 data = mf,
                 control = gamlss.control(trace = FALSE))
  alt <- update(null, formula. = ~ Age, what = "sigma")
  LRT <- deviance(null) - deviance(alt)
  p_lrt <- pchisq(LRT, df = 1, lower.tail = FALSE)
  est <- coef(alt, what = "sigma")["Age"]
  tibble(Term = "Age",
         Estimate = est,
         P_value = p_lrt)
}


meth <- get_m_values()
covs <- get_covs_884() %>% 
  mutate(Age = Age / 50)

snp_mat_list <- readRDS("./Data/RData/Genotypes/snp_mat_per_probe.rds") %>% compact()
snp_mat_list <- mclapply(snp_mat_list, function(x) x[covs$SUBJID, , drop = FALSE], mc.cores = 12)



meth <- meth[match(covs$SUBJID, meth$SUBJID), ]
meth <- meth[names(snp_mat_list)]

mf <- covs %>% 
  get_model_frame_15_props() %>% 
  mutate(Age = Age / 50)

fm <- y ~ Sex + Smoking_status + CMV_serostatus + ns(Age, df = 3) + Ancestry_PC1 + Ancestry_PC2


res <- mcmapply(fit_age_dispersion,
                meth[1:1000], 
                snp_mat_list[1:100],
                MoreArgs = list(fm = fm, mf = mf),
                SIMPLIFY = FALSE,
                mc.cores = 12) %>% 
  bind_rows(.id = "Probe") %>% 
  left_join(get_anno_meth_location())

saveRDS(res, "./Data/RData/Results/Dispersion/res_dispersion_gamlss_no_cells_884.rds")