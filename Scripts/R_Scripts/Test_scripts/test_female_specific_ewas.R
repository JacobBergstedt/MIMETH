#!/opt/gensoft/exe/R/3.6.2/scripts/Rscript
library(splines)
library(optparse)
library(dplyr)
library(parallel)
library(purrr)
library(sandwich)
library(lmtest)
library(lme4)
library(pbkrtest)
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_annotation.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")


p <- Sys.time()
cov <- "Age_of_menarche"
fm <- get_control_fm_15_cells() %>%
  update(paste0("~ . - Sex + ", cov))

# Load data
meth <- get_m_values()
covariates <- get_covs_884() %>% 
  mutate(Age = Age / 50)


fem_db <- readRDS("./Data/RData/Environment/women_specific_db.rds") %>%
  mutate(Progesterone_per_day_mg = log(Progesterone_per_day_mg),
         Estrogen_per_day_mg = log(Estrogen_per_day_mg))

covariates <- covariates %>% inner_join(fem_db, by = "SUBJID")

meth <- meth[meth$SUBJID %in% covariates$SUBJID, ]
meth <- meth[match(covariates$SUBJID, meth$SUBJID), ]




snp_mat_list <- readRDS("./Data/RData/Genotypes/snp_mat_per_probe.rds")
snp_mat_list <- map(snp_mat_list, ~ .[covariates$SUBJID, , drop = FALSE])

covariates <- covariates[c("Day_of_sampling", all.vars(fm))]
meth$SUBJID <- NULL

print("Data loaded")
print(Sys.time() - p)

are_cats <- keep(covariates, is.factor) %>%
  keep(~ length(levels(.)) > 2) %>%
  names()

cov_levels <- get_levels(covariates[[cov]])
fm <- update(fm, y ~ .)

terms <- paste0(cov, na.omit(cov_levels))
terms_tib <- tibble(Exposure = cov, Levels = cov_levels, Terms = terms)

snp_mat_4_debug <- snp_mat_list[["cg05738196"]]
fit_ewas_debug_random_effect_snp(fm_add_snps_and_random_effect(fm, snp_mat_4_debug),
                                 na.omit(cbind(y = meth[["cg05738196"]], covariates, snp_mat_4_debug)),
                                 terms_tib)

gc()

res <- mcmapply(infer_ewas_random_effect_snps,
                meth[1:100],
                snp_mat_list[1:100],
                MoreArgs = list(fm = fm, covariates = covariates, terms_tib = terms_tib),
                SIMPLIFY = FALSE,
                mc.cores = 1) %>%
  bind_rows(.id = "Probe")


res <- compute_adjusted_p_values(res)
saveRDS(res, paste0(chunk_path, cov, ".rds"))

