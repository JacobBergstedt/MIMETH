#!/opt/gensoft/exe/R/3.6.2/scripts/Rscript
#SBATCH --mem=120G
#SBATCH --cpus-per-task=6
#SBATCH --partition=geh
#SBATCH --qos=geh

library(optparse)
library(tidyverse)
library(parallel)
library(sandwich)
library(lmtest)
library(splines)
library(lme4)
library(pbkrtest)
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")

option_list = list(
  make_option("--interacting_var", action = "store", type = "character"),
  make_option("--n_cores", action = "store", type = "integer")
)

opt <- parse_args(OptionParser(option_list = option_list))
interacting_var <- opt$interacting_var
n_cores <- opt$n_cores

# Parse options
if (interacting_var != "Smoking_binary") {
  
  fm <- get_control_fm_15_cells() %>%
    update(paste0("~ . + Age * ", interacting_var)) %>% 
    rm_duplicate_age_terms()
    
} else {
  
  fm <- get_control_fm_15_cells() %>%
    update(~ Smoking_binary + . -Smoking_status) %>% 
    update(paste0("~ . + Age * ", interacting_var)) %>% 
    rm_duplicate_age_terms()
  
}


# Load input data
meth <- get_m_values()
snp_mat_list <- get_snp_mat_per_probe()

covariates <- get_covs_884() %>% 
  mutate(Smoking_binary = fct_collapse(Smoking_status, Non_Smoker = c("Non_Smoker","Ex_Smoker"))) %>% 
  mutate(Age = Age / 50)

snp_mat_list <- mclapply(snp_mat_list, function(x) x[covariates$SUBJID, , drop = FALSE], mc.cores = n_cores)
meth <- meth[match(covariates$SUBJID, meth$SUBJID),]
meth$SUBJID <- NULL
covariates <- covariates[c("Day_of_sampling", all.vars(fm))]

fm <- update(fm, y ~ .)

print(fm)

cov_levels <- get_levels(covariates[[interacting_var]])
cov_terms <- paste0(interacting_var, na.omit(cov_levels))
interaction_terms <- paste0(cov_terms, ":Age")

terms_tib <- tibble(Exposure = paste0("Age:", interacting_var), 
                    Terms = interaction_terms,
                    Interacting_term1 = "Age",
                    Interacting_term2 = interacting_var)


snp_mat_4_debug <- snp_mat_list[["cg05738196"]]
fit_ewas_debug_random_effect_snp(fm_add_snps_and_random_effect(fm, snp_mat_4_debug),
                                 na.omit(cbind(y = meth[["cg05738196"]], covariates, snp_mat_4_debug)),
                                 terms_tib)

res <- mcmapply(infer_ewas_random_effect_snps,
                meth,
                snp_mat_list,
                MoreArgs = list(fm = fm, covariates = covariates, terms_tib = terms_tib),
                SIMPLIFY = FALSE,
                mc.cores = n_cores) %>%
  bind_rows(.id = "Probe")

res <- compute_adjusted_p_values(res)

saveRDS(res, 
        paste0("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Age_interactions/age_", 
               tolower(interacting_var), ".rds"))


