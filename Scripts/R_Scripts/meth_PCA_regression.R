## Packages
library(MASS)
library(lme4)
library(tidyverse)
library(parallel)
library(lmerTest)
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
fit_model <- function(y, covariates = NULL, cells = NULL, day_of_sampling) {
  
  
  db <- data.frame(y = y, Day_of_sampling = day_of_sampling)
  if (!is.null(covariates)) db <- cbind(db, covariates)
  if (!is.null(cells)) db <- cbind(db, cells)
  
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
cell_list_16_cells <- get_cell_list()
cell_list <- cell_list_16_cells[!cell_list_16_cells %in% c("X_VIABLE_NEUTROPHILS_OF_TOTAL.panel7")]

meth <- get_m_values()

covariates <- get_covs_884() %>%
  select(SUBJID, all_of(cell_list), Age, Sex, Ancestry_PC1, Ancestry_PC2, CMV_serostatus, Smoking_status, Day_of_sampling, starts_with("IDOL")) %>% 
  mutate(Age = Age / 50)

cells <- covariates %>% 
  select(all_of(cell_list))

IDOL <- covariates %>% 
  select(starts_with("IDOL"), -contains("extended")) %>% 
  select(-IDOL_Neu)

IDOL_ext <- covariates %>% 
  select(starts_with("IDOL_extended")) %>% 
  select(-IDOL_extended_Neu)


meth <- meth[match(covariates$SUBJID, meth$SUBJID),]
meth_PCs <- get_50_PCs(meth)

covariates <- covariates %>% select(-all_of(cell_list), -SUBJID, -starts_with("IDOL"))
PCs <- meth_PCs$PCA
PC_prop_var <- meth_PCs$Prop_Var %>% enframe(name = "PC", value = "Prop_var")  


res_cells_covs <- map_dfr(PCs[-1], fit_model, covariates = covariates, cells = cells, day_of_sampling = covariates$Day_of_sampling, .id = "PC") %>% 
  left_join(PC_prop_var)

res_only_covs <- map_dfr(PCs[-1], fit_model, covariates = covariates, day_of_sampling = covariates$Day_of_sampling, .id = "PC") %>% 
  left_join(PC_prop_var)

res_IDOL_covs <- map_dfr(PCs[-1], fit_model, covariates = covariates, cells = IDOL, day_of_sampling = covariates$Day_of_sampling, .id = "PC") %>% 
  left_join(PC_prop_var)

res_IDOL_ext_covs <- map_dfr(PCs[-1], fit_model, covariates = covariates, cells = IDOL_ext, day_of_sampling = covariates$Day_of_sampling, .id = "PC") %>% 
  left_join(PC_prop_var)

res <- bind_rows(add_column(res_cells_covs, Adjustment = "Covariates and 16 cells"),
                 add_column(res_only_covs, Adjustment = "Only covariates"),
                 add_column(res_IDOL_covs, Adjustment = "Covariates and IDOL"),
                 add_column(res_IDOL_ext_covs, Adjustment = "Covariates and IDOL_ext"))


saveRDS(res, "./Data/RData/Results/meth_PCA_regression_non_compositional_with_DC.rds")

