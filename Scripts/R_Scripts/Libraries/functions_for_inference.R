
estimate_cell_mediation <- function(meth_site, model_frame, fm, fm_cells, adjustments, cell_adjustments) {
 
  db <- na.omit(cbind(y = meth_site, model_frame))
  
  m <- lm(fm_cells, db)
  m_total <- lm(fm, db)
  
  relimp_direct <- calc.relimp(m, type = "lmg", always = cell_adjustments)$lmg
  relimp_total <- calc.relimp(m_total, type = "lmg", always = adjustments)$lmg
  list(Direct = relimp_direct, Total = relimp_total)
  
}

get_model_frame <- function(db) {
  fm <- get_control_fm()
  get_all_vars(fm, db)
}

get_model_frame_15_props <- function(db) {
  fm <- get_control_fm_15_cells()
  get_all_vars(fm, db)
}

get_control_fm <- function(to_move_last = NULL) {
  
  move_last <- function(vec, to_move_last) {
    if (!to_move_last %in% vec) stop("Interacting variable not in the canonical proportion controls")
    c(vec[vec != to_move_last], to_move_last)
  }
  
  cell_props <- scan("./Data/prop_controls_16.txt", what = character())
  
  if (!is_null(to_move_last)) cell_props <- move_last(cell_props, to_move_last)
  
  " ~ Sex + Smoking_status + CMV_serostatus + ns(Age, df = 3) + Ancestry_PC1 + Ancestry_PC2 + " %>%
    paste0(paste0(cell_props, collapse = " + ")) %>%
    as.formula()
}

get_control_fm_15_cells <- function(to_move_last = NULL) {
  
  move_last <- function(vec, to_move_last) {
    if (!to_move_last %in% vec) stop("Interacting variable not in the canonical proportion controls")
    c(vec[vec != to_move_last], to_move_last)
  }
  
  cell_props <- get_15_props()
  
  if (!is_null(to_move_last)) cell_props <- move_last(cell_props, to_move_last)
  
  " ~ Sex + Smoking_status + CMV_serostatus + ns(Age, df = 3) + Ancestry_PC1 + Ancestry_PC2 + " %>%
    paste0(paste0(cell_props, collapse = " + ")) %>%
    as.formula()
}

get_total_effect_fm <- function(CpG = NULL) {
  if(is_null(CpG)) ~ Sex + Smoking_status + CMV_serostatus + ns(Age, df = 3) + Ancestry_PC1 + Ancestry_PC2 
  else as.formula(paste0(CpG, "~ Sex + Smoking_status + CMV_serostatus + ns(Age, df = 3) + Ancestry_PC1 + Ancestry_PC2"))
  
}

get_levels <- function(x) {
  if (!is.factor(x)) NA
  else levels(x)[-1]
}

compute_adjusted_p_values <- function(res) {
  res_distinct <- distinct(res, Probe, .keep_all = TRUE) %>% 
    mutate(P_FDR = p.adjust(P_value, "fdr"), P_bonferroni = p.adjust(P_value)) %>% 
    select(P_FDR, P_bonferroni, Probe, Exposure)
  full_join(res, res_distinct)
}

fit_example_model <- function(CpG, meth, covs, linear_age = FALSE) {
  db <- tibble(meth[c("SUBJID", CpG)]) %>% inner_join(covs)
  fm <- get_control_fm() %>% update(paste0(CpG, " ~ . + (1|Day_of_sampling)"))
  if (linear_age) {
    fm <- fm %>% update(. ~ . + Age) %>% rm_duplicate_age_terms()
  }
  lmer(fm, db)
}

fit_example_model_15_cells <- function(CpG, meth, covs, snp_mat = NULL, linear_age = FALSE) {
  
  
  fm <- get_control_fm_15_cells() %>% 
    update(as.formula(paste0(CpG, " ~ . "))) %>% 
    fm_add_snps_and_random_effect(snp_mat)
  
  db <- tibble(meth[c("SUBJID", CpG)]) %>% inner_join(covs)
  
  if (!is_null(snp_mat)) {
    
    db <- inner_join(db, add_column(SUBJID = rownames(snp_mat), as.data.frame(snp_mat)))
  }
  
  if (linear_age) {
    fm <- fm %>% update(. ~ . + Age) %>% rm_duplicate_age_terms()
  }
  
  lmer(fm, db, na.action = na.exclude)
  
}


fit_ewas_random_effect <- function(fm, db, terms_tib) {
  
  m <- lmer(fm, db, 
            control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = formals(isSingular)$tol)))
  
  m_null <- update(m, paste0(". ~ . - ", terms_tib$Exposure[1]))
  p_val <- KRmodcomp(m, m_null)$stats$p.value
  res <- as_tibble(summary(m)$coefficients[unique(terms_tib$Terms), c("Estimate", "Std. Error"), drop = FALSE])
  names(res) <- c("Estimate", "Standard_error")
  res$P_value <- p_val
  bind_cols(res, terms_tib)
}

fit_ewas_debug <- function(fm, db, var_fun, terms_tib) {
  m <- lm(fm, db)
  res <- coeftest(m, vcov. = var_fun)
  print(summary(m))
  print(res)
}

fit_ewas_debug_random_effect <- function(fm, db, terms_tib) {
  m <- lmer(fm, db)
  m_null <- update(m, paste0(". ~ . - ", terms_tib$Exposure[1]))
  p_val <- KRmodcomp(m, m_null)$stats$p.value
  print(summary(m))
  print(KRmodcomp(m, m_null))
}

fit_ewas_debug_random_effect_snp <- function(fm, db, terms_tib) {
  m <- lmer(fm, db)
  m_null <- update(m, paste0(". ~ . - ", terms_tib$Exposure[1]))
  p_val <- KRmodcomp(m, m_null)$stats$p.value
  print(summary(m))
  print(KRmodcomp(m, m_null))
}

fit_EWAS_full <- function(cov, null_fm, covariates, meth, snp_mat_list, n_cores) {
  
  fm <- null_fm %>%
    update(paste0("~ ", cov, " + .")) %>%
    rm_duplicate_age_terms()
  
  meth$SUBJID <- NULL
  covariates <- covariates[c("Day_of_sampling", all.vars(fm))]
  
  # Logic for special cases
  if (cov == "Years_smoking") {
    
    fm <- update(fm, y ~ . - Smoking_status)
    ids <- !is.na(covariates$Smoking_status) & covariates$Smoking_status == "Smoker"
    covariates <- covariates[ids, ]
    covariates$Years_smoking <- covariates$Years_smoking / 20
    meth <- meth[ids, ]
    snp_mat_list <- map(snp_mat_list, ~ .[ids, , drop = FALSE])
    
  } else if (cov == "Years_since_last_smoke") {
    
    fm <- update(fm, y ~ . - Smoking_status)
    ids <- !is.na(covariates$Smoking_status) & covariates$Smoking_status == "Ex_Smoker"
    covariates <- covariates[ids, ]
    covariates$Years_since_last_smoke <- covariates$Years_since_last_smoke / 20
    meth <- meth[ids, ]
    snp_mat_list <- map(snp_mat_list, ~ .[ids, , drop = FALSE])
    
  } else if (cov == "CMV_serology") {
    
    fm <- update(fm, y ~ . - CMV_serostatus)
    ids <- !is.na(covariates$CMV_serostatus) & covariates$CMV_serostatus == "Positive"
    covariates <- covariates[ids, ]
    covariates$CMV_serology <- covariates$CMV_serology / 500
    meth <- meth[ids, ]
    snp_mat_list <- map(snp_mat_list, ~ .[ids, , drop = FALSE])
      
  } else {
    
    fm <- update(fm, y ~ .)
    
  }
  
  terms_tib <- get_terms_tib(cov, covariates)
  
  
  
  res <- mcmapply(infer_ewas_random_effect_snps,
                  meth,
                  snp_mat_list,
                  MoreArgs = list(fm = fm, covariates = covariates, terms_tib = terms_tib),
                  SIMPLIFY = FALSE,
                  mc.cores = n_cores)
  
  res <- bind_rows(res, .id = "Probe")
  
  compute_adjusted_p_values(res)
  
  
}


fm_add_snps <- function(fm, snp_mat) {
  if (!is.null(snp_mat)) {
    update_string <- as.formula(paste0(". ~ . + ", paste0(colnames(snp_mat), collapse = " + ")))
    fm <- update(fm, update_string)
  }
  fm
}


fm_add_snps_and_random_effect <- function(fm, snp_mat) {
  update_string <- paste0(". ~ . + ", paste0(colnames(snp_mat), collapse = " + "), " + (1 | Day_of_sampling)")
  update(fm, update_string)
}

get_terms_tib <- function(cov, covariates) {
  cov_levels <- get_levels(covariates[[cov]])
  terms <- paste0(cov, na.omit(cov_levels))
  tibble(Exposure = cov, Levels = cov_levels, Terms = terms)
}

infer_ewas_random_effect_snps <- function(probe, snps, fm, covariates, terms_tib) {
  fm <- fm_add_snps_and_random_effect(fm, snps)
  if (!is.null(snps)) db <- na.omit(cbind(y = probe, covariates, snps)) 
  else db <- na.omit(cbind(y = probe, covariates)) 
  fit_ewas_random_effect(fm, db, terms_tib)
}

infer_ewas_random_effect <- function(probe, fm, covariates, terms_tib) {
  db <- na.omit(cbind(y = probe, covariates)) 
  fit_ewas_random_effect(fm, db, terms_tib)
}

infer4chunk_control4selection <- function(probe, selected_proportions, vcov.) {
  fm_base <- as.formula(paste0("y ~ ",
                               paste0(selected_proportions, collapse = " + "),
                               " + Age + CMVPositiveSerology + MCHC + Sex + Smoking"))
  db <- cbind(y = probe, controls, proportions[selected_proportions])
  res <- lapply(cov_names, infer4cov_control4selection, fm = fm_base, db = db, vcov.)
  bind_rows(res)
}

infer4cov_control4selection <- function(cov, fm_base, db, vcov.){
  
  if (cov == "control_model") {
    
    res <- fit_ewas_controls(fm_base, db, vcov.)
    
  } else if (cov %in% are_cats) {
    
    db <- cbind(exposure = covs[[cov]], db)
    res <- fit_ewas_cats(fm_base, db, vcov.)
    res$Exposure <- cov
    res
    
  } else {
    
    fm <- update(fm_base, . ~ . + exposure)
    db <- cbind(exposure = covs[[cov]], db)
    res <- fit_ewas(fm, db, vcov.)
    res$Exposure <- cov
    res
    
  }
}

rm_duplicate_age_terms <- function(fm) {
  vars <- formula.tools::rhs.vars(fm)
  if (sum(grepl("Age", vars)) > 1) fm <- update(fm, paste0("~ . - ", vars[grepl("ns", vars)]))
  fm
}

compare_AIC <- function(probe, fm, mf) {
  mf$y <- probe
  m_normal <- gamlss(fm, data = na.omit(mf), family = NO())
  m_logit_normal <- gamlss(fm, data = na.omit(mf), family = LOGITNO())
  tibble(AIC_normal = AIC(m_normal), AIC_logit_normal = AIC(m_logit_normal))
}


two_stage_adjust <- function(tib, 
                             bottom_level = "SNP", 
                             top_level = "Probe", 
                             Bonferroni_bottom_adjustment, 
                             level = 0.01, 
                             manual_top_adjustment_factor = NULL, 
                             manual_local_adjustment_factor = NULL) {
  
  # Nr snps: 5699237
  # Nr probes: 643,546
  
  if ("Pvalue" %in% names(tib)) tib <- rename(tib, P_value = Pvalue)
  tib <- group_by(tib, .data[[top_level]])
  
  if (is_null(manual_local_adjustment_factor)) tib <- mutate(tib, P_local = p.adjust(P_value))
  else tib <- mutate(tib, P_local = pmin(1, P_value * .env[["manual_local_adjustment_factor"]]))
  
  # Summary tibble 
  tib_summary <- summarize(tib, P_intersection = min(P_local))
  
  if (is_null(manual_top_adjustment_factor)) tib_summary <- mutate(tib_summary, P_intersection_adjusted = p.adjust(P_intersection, "BH"))
  else tib_summary <- mutate(tib_summary, P_intersection_adjusted = pmin(P_intersection * .env[["manual_top_adjustment_factor"]], 1))
  
  tib_summary <- mutate(tib_summary, Significant_family = P_intersection_adjusted < level)
  adj_factor <- nrow(tib_summary) / sum(tib_summary$Significant_family)
  
  # Rejoin with original tibble
  tib <- left_join(ungroup(tib), tib_summary, by = top_level)
  
  # Work with significant tibble
  tib_sign <- tib %>% 
    filter(Significant_family) %>%
    group_by(.data[[top_level]])
  
  if (Bonferroni_bottom_adjustment) tib_sign <- mutate(tib_sign, P_adjusted = pmin(1, adj_factor * P_local))
  else tib_sign <- mutate(tib_sign, P_adjusted = pmin(1, p.adjust(P_value, "BH") * adj_factor))
  
  # Get significant local probes and rejoin with original tibble
  tib_sign %>% 
    mutate(Locally_significant = P_adjusted < level) %>% 
    ungroup() %>% 
    select(all_of(c(top_level, bottom_level)), P_adjusted, Locally_significant) %>% 
    right_join(tib, by = c(top_level, bottom_level))
  
}

make_smoking_binary <- function(tib) {
  
  tib %>% 
    mutate(Smoking_status = fct_recode(Smoking_status, No = "Non_Smoker", No = "Ex_Smoker", Yes = "Smoker")) %>% 
    dplyr::rename(Smoker = Smoking_status)
  
}

# two_stage_adjust <- function(tib, bottom_level = "SNP", top_level = "Probe", local_adjustment, level = 0.01, is_trans = FALSE) {
# 
#   if ("Pvalue" %in% names(tib)) tib <- rename(tib, P_value = Pvalue)
#   if (is_trans) {
#     
#     tib <- tib %>% 
#       group_by(.data[[top_level]]) %>% 
#       mutate(P_local = pmin(1, P_value * 5699237))
#     
#   } else {
#     
#     tib <- tib %>% 
#       group_by(.data[[top_level]]) %>% 
#       mutate(P_local = p.adjust(P_value))
#     
#   }
#   
#   tib_summary <- tib %>%
#     summarize(P_intersection = min(P_local)) %>% 
#     mutate(P_intersection_adjusted = p.adjust(P_intersection, "BH"),
#            Significant_family = P_intersection_adjusted < 0.01) %>% 
#     ungroup()
#   
#   adj_factor <- nrow(tib_summary) / sum(tib_summary$Significant_family)
#   tib <- left_join(ungroup(tib), tib_summary, by = top_level)
#   tib_sign <- tib %>%
#     filter(Significant_family) %>% 
#     group_by(.data[[top_level]])
#     
#   if (local_adjustment == "Bonferroni") {
#     
#     tib_sign <- mutate(tib_sign, P_adjusted = pmin(1, adj_factor * P_local), Locally_significant = P_adjusted < level)
#     
#   } else if (local_adjustment == "BH") {
#     
#     tib_sign <- mutate(tib_sign, P_adjusted = pmin(1, p.adjust(P_value, "BH") * adj_factor), Locally_significant = P_adjusted < level)
#     
#   }
#   
#   tib_sign <- tib_sign %>%
#     ungroup() %>%
#     select(top_level, bottom_level, P_adjusted, Locally_significant)
#   
#   left_join(tib, tib_sign, by = c(top_level, bottom_level))
# }


# two_stage_adjust_meQTL_results <- function(tib, type) {
#   
#   tib <- tib %>% lazy_dt()
#   
#   if (type == "trans") {
#     
#     tib <- tib %>% 
#       group_by(Probe) %>% 
#       mutate(P_local = pmin(1, Pvalue * 5699237)) %>% 
#       ungroup()
#     
#   } else {
#     
#     tib <- tib %>% 
#       group_by(Probe) %>% 
#       mutate(P_local = p.adjust(Pvalue)) %>% 
#       ungroup()
#     
#   }
#   
#   tib_summary <- tib %>% 
#     group_by(Probe) %>% 
#     summarize(P_intersection = min(P_local)) %>% 
#     mutate(P_intersection_adjusted = p.adjust(P_intersection, "BH"),
#            Is_significant = P_intersection_adjusted < 0.01) %>% 
#     ungroup()
#   
#   real_tib_summary <- as_tibble(tib_summary)
#   R <- sum(real_tib_summary$Is_significant)
#   adj_factor <- nrow(real_tib_summary) / R
#   
#   tib <- left_join(tib, tib_summary, by = "Probe")
#   
#   tib_sign <- tib %>% 
#     filter(Is_significant) %>% 
#     group_by(Probe) %>% 
#     mutate(P_adjusted = pmin(1, P_local * adj_factor)) %>% 
#     ungroup() %>%
#     select(Probe, SNP, P_adjusted)
#     
#   left_join(tib, tib_sign, by = c("Probe", "SNP")) %>% as_tibble()
# }
#
# two_stage_adjust_ewas_results <- function(tib) {
#   
#   tib <- tib %>% 
#     lazy_dt() %>% 
#     group_by(Exposure) %>% 
#     mutate(P_local = p.adjust(P_value, "BH")) %>% 
#     ungroup()
#   
#   
#   tib_summary <- tib %>% 
#     group_by(Exposure) %>% 
#     summarize(P_intersection = min(P_local)) %>% 
#     mutate(P_intersection_adjusted = p.adjust(P_intersection, "BH"),
#            Is_significant = P_intersection_adjusted < 0.01) %>% 
#     ungroup()
#   
#   real_tib_summary <- as_tibble(tib_summary)
#   R <- sum(real_tib_summary$Is_significant)
#   adj_factor <- nrow(real_tib_summary) / R
#   
#   tib <- left_join(tib, tib_summary, by = "Exposure")
#   
#   tib_sign <- tib %>% 
#     filter(Is_significant) %>% 
#     group_by(Exposure) %>% 
#     mutate(P_adjusted = pmin(1, P_local * adj_factor)) %>% 
#     ungroup() %>%
#     select(Probe, Exposure, P_adjusted)
#   
#   left_join(tib, tib_sign, by = c("Probe", "Exposure")) %>% as_tibble()
# }

# map_meth2cisgenotype <- function(probe, snps) {
#   fit_snp <- function(genotype) {
#     db <- na.omit(cbind(y = probe, covariates, x = genotype))
#     fit_ewas(fm, db, var_fun = var_fun, terms_tib = terms_tib)
#   }
#   genotype4probe <- genotypes[ , snps, drop = FALSE]
#   apply(genotype4probe, 2, fit_snp) %>% bind_rows(.id = "SNP")
# }

# setup_proportion_interaction_inference <- function(prop, cov, fm) {
#   fm <- update(fm, paste0(". ~ . + ", cov, " * ", prop))
#   cov_levels <- get_levels(covariates[[cov]])
#   cov_terms <- paste0(cov, na.omit(cov_levels))
#   interaction_terms <- paste0(cov_terms, ":", prop)
#   terms_tib <- tibble(Terms = interaction_terms,
#                       Exposure = interaction_terms,
#                       Levels = cov_levels)
#   list(fm = fm, terms_tib = terms_tib)
# }

# select_immunophenotypes <- function(probe) {
#   keep <- !is.na(probe)
#   props <- props[keep, ]
#   probe <- probe[keep]
#   sel <- stabsel(props, probe,
#                  intercept = TRUE,
#                  fitfun = glmnet.lasso,
#                  args.fitfun = list(alpha = alpha),
#                  q = q,
#                  B = B,
#                  papply = lapply,
#                  PFER = tol)
#   names(sel$selected)
# }

# compare_ewas_artifacts <- function() {
#   
#   collect_from_dir <- function(path) {
#     
#     read_artifact <- function(path, name) {
#       artifact <- readRDS(path)
#       names(artifact)[2] <- name
#       artifact
#     }
#     
#     artifact_files <- list.files(path)
#     good_files <- tibble(String = artifact_files, 
#            Date = as_date(str_extract(artifact_files, "2019-[0-9]{2}-[0-9]{2}")), 
#            Control = str_extract(artifact_files, "(?<=correct_for_).*(?=_2019)")) %>% 
#       group_by(Control) %>% 
#       summarize(Good_file = String[which.max(Date)])  
#     
#     map2(paste0(path, good_files$Good_file),  good_files$Control, read_artifact) %>% 
#         reduce(left_join)
#   }
#   
#   specify_versions <- function(tab, prefix) {
#     names(tab)[-1] <- paste0(prefix, "_", names(tab)[-1])
#     tab
#   }
#   
#   versions <- c("Constant_variance", "Sandwich")
#   paste0("./Data/RData/Results/EWAS/M_values/", versions, "/Environment/Artifacts/") %>% 
#     map(collect_from_dir) %>% 
#     map2(versions, specify_versions) %>% 
#     reduce(left_join)
# }

# # Functions to estimate cell specific methylation rate --------------------
# 
# build_stan_res_summary <- function(post_samples, fit_sum) {
#   tibble(Cell = major_props, 
#          P_less_than_0.2 = colMeans2(post_samples < 0.2),
#          P_higher_than_0.8 = colMeans2(post_samples > 0.8)) %>% 
#     bind_cols(as_tibble(fit_sum$summary[1:16, c("mean", "sd", "2.5%", "97.5%")]))
# }
# 
# post_comparisons <- function(post_samples) {
#   comparisons_spec <- rbind(unique_combs(major_props[1], major_props[-1]), 
#                             c("X_CD4_naive_of_total.panel1", "X_CD4_CM_of_total.panel1"),
#                             c("X_CD4_naive_of_total.panel1", "X_CD4_EM_of_total.panel1"),
#                             c("X_CD4_naive_of_total.panel1", "X_CD4_EMRA_of_total.panel1"),
#                             c("X_CD8_naive_of_total.panel1", "X_CD8_CM_of_total.panel1"),
#                             c("X_CD8_naive_of_total.panel1", "X_CD8_EM_of_total.panel1"),
#                             c("X_CD8_naive_of_total.panel1", "X_CD8_EMRA_of_total.panel1"))
#   colnames(comparisons_spec) <- c("Var1", "Var2")
#   n_comps <- nrow(comparisons_spec)
#   
#   comparisons <- matrix(ncol = n_comps, nrow = nrow(post_samples))
#   for (i in 1:n_comps) {
#     comparisons[, i] <- post_samples[, comparisons_spec[i, 1]] - post_samples[, comparisons_spec[i, 2]]
#   }
#   CIs <- colQuantiles(comparisons, probs = c(0.025, 0.975))
#   
#   data.frame(Var1 = comparisons_spec[, 1],
#              Var2 = comparisons_spec[, 2],
#              Post_median = colMedians(comparisons),
#              BF = apply(comparisons, 2, function(x) bayesfactor_parameters(x, prior = distribution_normal(1e3, 0, 0.5))$BF),
#              lower = CIs[, 1],
#              upper = CIs[, 2])
#   
# }
# 
# fit_stan <- function(y, model, stan_data) {
#   
#   stan_data$y <- y
#   fit <- sampling(model, data = stan_data, chains = 4, iter = 2000, cores = 4, refresh = 0)
#   post_samples <- as.matrix(fit)[, 1:16]
#   colnames(post_samples) <- colnames(stan_data$x)
#   stan_res <- build_stan_res_summary(post_samples, summary(fit))
#   comparisons_res <- post_comparisons(post_samples)
#   list(stan_res = stan_res, 
#        comparisons_res = comparisons_res)
#   
# }
# two_stage_adjust <- function(tib, bottom_level = "SNP", top_level = "Probe", local_adjustment, level = 0.01, is_trans = FALSE) {
#   
#   if ("Pvalue" %in% names(tib)) tib <- rename(tib, P_value = Pvalue)
#   
#   tib <- tib %>% 
#     rename(top_level = .data[[top_level]], bottom_level = .data[[bottom_level]]) %>% 
#     lazy_dt() %>% 
#     group_by(top_level)
#   
#   if (is_trans) {
#     
#     tib <- tib %>% 
#       mutate(P_local = pmin(1, P_value * 5699237))
#     
#   } else {
#     
#     tib <- tib %>%
#       mutate(P_local = p.adjust(P_value))
#     
#   }
#   
#   tib_summary <- tib %>%
#     summarize(P_intersection = min(P_local)) %>% 
#     mutate(P_intersection_adjusted = p.adjust(P_intersection, "BH"),
#            Significant_family = P_intersection_adjusted < 0.01) %>% 
#     ungroup()
#   
#   real_tib_summary <- as_tibble(tib_summary)
#   R <- sum(real_tib_summary$Significant_family)
#   adj_factor <- nrow(real_tib_summary) / R
#   
#   tib <- left_join(ungroup(tib), tib_summary, by = "top_level")
#   tib_sign <- tib %>%
#     filter(Significant_family) %>% 
#     group_by(top_level)
#   
#   if (local_adjustment == "Bonferroni") {
#     tib_sign <- mutate(tib_sign, P_adjusted = pmin(1, P_local * adj_factor), Locally_significant = P_adjusted < level)
#   } else if (local_adjustment == "BH") {
#     tib_sign <-  mutate(tib_sign, P_adjusted = pmin(1, p.adjust(P_value, "BH") * adj_factor), Locally_significant = P_adjusted < level)
#   }
#   
#   tib_sign <- tib_sign %>% 
#     ungroup() %>%
#     select(top_level, bottom_level, P_adjusted, Locally_significant)
#   
#   left_join(tib, tib_sign, by = c("top_level", "bottom_level")) %>% 
#     as_tibble() %>% 
#     rename("{top_level}" := top_level, "{bottom_level}" := bottom_level)
#   
# }
# infer_interactions <- function(probe, covariates, props_utils) {
#   map_dfr(props_utils,
#           ~ infer_ewas(probe, .$fm, covariates, var_fun = vcovHC, terms_tib = .$terms_tib))
# }
# 
# infer_ewas <- function(probe, fm, covariates, var_fun, terms_tib) {
#   db <- cbind(y = probe, covariates)
#   fit_ewas(fm, db, var_fun, terms_tib)
# }
# 
# infer_ewas_cats <- function(probe, fm, covariates, var_fun, terms_tib) {
#   db <- cbind(y = probe, covariates)
#   db <- na.omit(db)
#   fit_ewas_cats(fm, db, var_fun, terms_tib)
# }
# 
# infer_ewas_random_effect <- function(probe, fm, covariates, terms_tib) {
#   db <- na.omit(cbind(y = probe, covariates))
#   fit_ewas_random_effect(fm, db, terms_tib)
# }
# 
# fit_ewas <- function(fm, db, var_fun, terms_tib) {
#   m <- lm(fm, db)
#   res <- coeftest(m, vcov. = var_fun)[terms_tib$Terms, c("Estimate", "Std. Error", "Pr(>|t|)"), drop = FALSE]
#   res <- as_tibble(res)
#   names(res) <- c("Estimate", "Standard_error", "P_value")
#   res <- bind_cols(res, terms_tib)
#   res
# }
# 
# fit_ewas_cats <- function(fm, db, var_fun, terms_tib) {
#   m <- lm(fm, db)
#   m_null <- update(m, paste0(". ~ . - ", terms_tib$Exposure[1]))
#   res <- coeftest(m, vcov. = var_fun)[terms_tib$Terms, c("Estimate", "Std. Error")]
#   res <- as_tibble(res)
#   names(res) <- c("Estimate", "Standard_error")
#   res <- bind_cols(res, terms_tib[c("Exposure", "Levels")])
#   res$P_value <- waldtest(m_null, m, vcov = var_fun)[["Pr(>F)"]][2]
#   res
# }



