## Packages
library(MASS)
library(lme4)
library(lmerTest)
library(tidyverse)
library(parallel)
library(splines)
library(pbkrtest)
library(sandwich)
library(broom)
library(irlba)
library(readxl)
source("./Scripts/R_scripts/Libraries/functions_for_mediation2.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_compositions.R")
library(robCompositions)
select <- dplyr::select
fit_compositional_model <- function(y, covariates = NULL, balances = NULL, day_of_sampling) {
  
  db <- data.frame(y = y, Day_of_sampling = day_of_sampling)
  if (!is.null(covariates)) db <- cbind(db, covariates)
  if (!is.null(balances)) db <- cbind(db, balances)
  
  fm <- as.formula(paste0("y ~ ", paste0(names(select(db, -y, -Day_of_sampling)), collapse = " + "), " + (1|Day_of_sampling)"))
  m <- lmer(fm, db)
  res <- coef(summary(m)) 
  res[-1, ] %>% 
    as.data.frame() %>% 
    rownames_to_column("Variable") %>% 
    as_tibble() %>% 
    select(Variable, Estimate, Standard_error = "Std. Error", P_value = "Pr(>|t|)")
  
}

get_50_PCs <- function(x, n_cores = 10) {
  
  impute_with_mean <- function(probe) {
    probe[is.na(probe)] <- mean(probe, na.rm = TRUE)
    probe
  }
  
  id <- x$SUBJID
  x$SUBJID <- NULL
  x <- mclapply(x, impute_with_mean, mc.cores = n_cores)
  x <- do.call("cbind", x)
  PCA <- prcomp_irlba(x, n = 50, center = TRUE, scale. = TRUE)
  list(Prop_Var = summary(PCA)$importance["Proportion of Variance", ], 
       PCA = bind_cols(SUBJID = id, as_tibble(PCA$x)))
  
}


# Run analysis ------------------------------------------------------------

## Data
cell_list <- get_cell_list()
meth <- get_m_values()

covariates <- get_covs_884() %>%
  select(SUBJID, all_of(cell_list), Age, Sex, Ancestry_PC1, Ancestry_PC2, CMV_serostatus, Day_of_sampling, Smoking_status, starts_with("IDOL")) %>% 
  filter(!if_any(all_of(cell_list), .fns = ~ . <= 0))

cells <- covariates %>% 
  select(all_of(cell_list))

IDOL <- covariates %>% 
  select(starts_with("IDOL")) %>% select(!contains("extended"))


meth <- meth[match(covariates$SUBJID, meth$SUBJID),]
meth_PCs <- get_50_PCs(meth)

covariates <- covariates %>% select(-all_of(cell_list), -SUBJID, -starts_with("IDOL"))
PCs <- meth_PCs$PCA
PC_prop_var <- meth_PCs$Prop_Var %>% enframe(name = "PC", value = "Prop_var")  

cell_balances <- balance_preds(cells, sbp = read_xlsx("./Data/Cell_subset_partitions_3_with_DC.xlsx"))
IDOL_balances <- balance_preds(IDOL, sbp = read_xlsx("./Data/IDOL_partitions.xlsx"))

res_only_cells <- map_dfr(PCs[-1], fit_compositional_model, balances = cell_balances, day_of_sampling = covariates$Day_of_sampling, .id = "PC") %>% 
  left_join(PC_prop_var)

res_cells_covs <- map_dfr(PCs[-1], fit_compositional_model, covariates = covariates, balances = cell_balances, day_of_sampling = covariates$Day_of_sampling, .id = "PC") %>% 
  left_join(PC_prop_var)

res_only_covs <- map_dfr(PCs[-1], fit_compositional_model, covariates = covariates, day_of_sampling = covariates$Day_of_sampling, .id = "PC") %>% 
  left_join(PC_prop_var)

res_only_IDOL <- map_dfr(PCs[-1], fit_compositional_model, covariates = NULL, balances = IDOL_balances, day_of_sampling = covariates$Day_of_sampling, .id = "PC") %>% 
  left_join(PC_prop_var)

res_IDOL_covs <- map_dfr(PCs[-1], fit_compositional_model, covariates = covariates, balances = IDOL_balances, day_of_sampling = covariates$Day_of_sampling, .id = "PC") %>% 
  left_join(PC_prop_var)

res <- bind_rows(add_column(res_only_cells, Adjustment = "Only 16 cells"),
                 add_column(res_cells_covs, Adjustment = "Covariates and 16 cells"),
                 add_column(res_only_covs, Adjustment = "Only covariates"),
                 add_column(res_only_IDOL, Adjustment = "IDOL"),
                 add_column(res_IDOL_covs, Adjustment = "Covariates and IDOL"))


saveRDS(res, "./Data/RData/Results/meth_PCA_regression_results_with_random_effect_with_DCs.rds")

