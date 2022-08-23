estimate_effect_on_cells <- function(fm, db, treatment, cell_list, n_sim) {
  
  m <- set_names(fm, cell_list) %>% 
    map(~ lm(., db))
  
  beta <-  map_dbl(m, ~ coefficients(.)[treatment])
  sims <- vapply(m, simulate_beta, double(n_sim), n_sim = n_sim, vars = treatment)
  list(sims = sims, beta = beta[cell_list])
  
}

estimate_total_effect <- function(fm, db, treatment, n_sim) {
  
  m <- lm(fm, db)
  total <- coefficients(m)[treatment]
  sims <- simulate_beta(m, n_sim, treatment)
  list(sims = sims, beta = total)
  
}

estimate_cell_effect <- function(fm, db, cell_list, n_sim) {
  
  m <- lm(fm, db)
  beta <- coefficients(m)
  sims <- simulate_beta(m, n_sim, vars = cell_list)
  list(sims = sims, beta = beta[cell_list])
  
}

estimate_mediation <- function(fm_list, db, treatment, effect_on_cells, cell_list, n_sim) {
  
  cell_effect <- estimate_cell_effect(fm_list$fm_direct, db, cell_list = cell_list, n_sim = n_sim)
  te <- estimate_total_effect(fm_list$fm_total, db, treatment, n_sim = n_sim)
  
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

run_mediation_for_site <- function(meth_site, fm_list, mf, treatment, effect_on_cells, cell_list = cell_list, n_sim = n_sim) {
  
  db <- na.omit(cbind(y = meth_site, mf))
  estimate_mediation(fm_list, db, treatment, effect_on_cells = effect_on_cells, cell_list = cell_list, n_sim = n_sim)
  
}

run_mediation <- function(meth, fm_list, mf, treatment, n_sim, cell_list, n_cores = 12) {

  
  effect_on_cells <- estimate_effect_on_cells(fm_list$fm_cell_responses,  
                                              db = mf, 
                                              treatment = treatment, 
                                              cell_list = cell_list, 
                                              n_sim = n_sim)
  mclapply(meth, run_mediation_for_site, 
           fm_list = fm_list, mf = mf, treatment = treatment, effect_on_cells = effect_on_cells, cell_list = cell_list, n_sim = n_sim, mc.cores = n_cores)
  
}

