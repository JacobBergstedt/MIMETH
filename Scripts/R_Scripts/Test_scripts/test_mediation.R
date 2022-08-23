
## Packages

library(MASS)
library(tidyverse)
library(parallel)

# diff <- left_join(dplyr::select(age_ewas_no_corr, Probe, Total = Estimate), dplyr::select(age_ewas, Probe, Direct = Estimate)) %>% left_join(dplyr::select(res_mediation_age, Probe, Mediation), by = "Probe")
# diff <- diff %>% mutate(Diff = Total - Direct) %>% arrange(desc(Diff))


## Functions

estimate_effect_on_cells <- function(db, treatment) {
  
  m <- set_names(fm_cells, cell_list_16_cells) %>% 
    map(~ lm(., db))
  
  beta <-  map_dbl(m, ~ coefficients(.)[treatment])
  sims <- vapply(m, function(x) mvrnorm(1e3, coefficients(x), vcov(x))[, treatment], double(1000))
  list(sims = sims, beta = beta[cell_list_16_cells])
  
}

estimate_total_effect <- function(db, treatment) {
  m <- lm(fm_total, db)
  total <- coefficients(m)[treatment]
  sims <- mvrnorm(1000, coefficients(m), vcov(m))[, treatment]
  list(sims = sims, beta = total)
}
  
estimate_cell_effect <- function(db, treatment) {
  
  m <- lm(fm_direct, db)
  beta <- coefficients(m)
  sims <- mvrnorm(1000, beta, vcov(m))[, cell_list_16_cells]
  list(sims = sims, beta = beta[cell_list_16_cells])
  
}

estimate_mediation <- function(db, treatment) {
  
  cell_effect <- estimate_cell_effect(db, treatment)
  te <- estimate_total_effect(db, treatment)
  effect_on_cells <- estimate_effect_on_cells(db, treatment)
  
  med <- sum(cell_effect$beta * effect_on_cells$beta)
  sim_med <- rowSums(cell_effect$sims * effect_on_cells$sims)
  
  if (sign(med) == sign(te$beta)) {
    
    ratio_med <- if_else(abs(med) > abs(te$beta), 1, med / te$beta)
    sim_ratio <- sim_med / te$sims  
    p_ratio <- 2 * pnorm(abs(ratio_med / sd(sim_ratio)), lower.tail = FALSE)
    se_ratio <- sd(sim_ratio)
    
  } else {
    
    p_ratio <- NA
    se_ratio <- NA
    ratio_med <- NA
    
  }
  
  tibble(Effect = c("Mediation", "Total", "Ratio_mediation"),
         Estimate = c(med, te$beta, ratio_med),
         Standard_error = c(sd(sim_med), sd(te$sims), se_ratio),
         P_value = c(2 * pnorm(abs(med / sd(sim_med)), lower.tail = FALSE),
                     2 * pnorm(abs(te$beta / sd(te$sims)), lower.tail = FALSE),
                     p_ratio))
  
}

run_mediation <- function(meth_site, mf, treatment) {
  
  db <- na.omit(cbind(y = meth_site, mf))
  estimate_mediation(db, treatment)
  
}

## Globals

cell_list_16_cells <- scan("./Data/prop_controls_16.txt", what = character())
fm_cells <- map(cell_list_16_cells, ~ paste0(., " ~ Age + Sex + CMV_serostatus + Smoking_status + Ancestry_PC1 + Ancestry_PC2"))
fm_total <- y ~ Age + Sex + CMV_serostatus + Smoking_status + Ancestry_PC1 + Ancestry_PC2
fm_direct <- get_control_fm() %>% update(y ~ Age + .) %>% rm_duplicate_age_terms()

## Data

mf <- get_covs()[c("Sex", "Smoking_status", "Day_of_sampling", "Ancestry_PC1", "Ancestry_PC2", "CMV_serostatus", "Age", cell_list_16_cells)] %>% 
  mutate(Age = Age / 50)
meth <- get_m_values()
meth["SUBJID"] <- NULL


## Main

res_age <- mclapply(meth, run_mediation, mf = mf, treatment = "Age", mc.cores = 12) %>% 
  bind_rows(.id = "Probe") %>% 
  add_column(Variable = "Age")

# res_CMV <- mclapply(meth, run_mediation, mf = mf, treatment = "CMV_serostatusPositive", mc.cores = 12) %>% 
#   bind_rows(.id = "Probe") %>% 
#   add_column(Variable = "CMV")
# 
# res_sex <- mclapply(meth, run_mediation, mf = mf, treatment = "SexFemale", mc.cores = 12) %>% 
#   bind_rows(.id = "Probe") %>% 
#   add_column(Variable = "Sex")
# 
# res_smoking <- mclapply(meth, run_mediation, mf = mf, treatment = "Smoking_statusSmoker", mc.cores = 12) %>% 
#   bind_rows(.id = "Probe") %>% 
#   add_column(Variable = "Smoking")
# 
# 
# res <- bind_rows(res_age, res_CMV, res_sex, res_smoking)
