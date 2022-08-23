
generate_cv_indices <- function() {
  ids <- get_meth_ids() 
  time <- gsub("-|:| ", "", Sys.time())
  path <- glue("./Data/RData/Cross_val_indices/cross_validation_indices_{time}.rds")
  saveRDS(sample.int(5, length(ids), replace = TRUE), path)
  path
}

compute_R2 <- function(obs, pred) {
  cor(obs, pred) ^ 2
}

setup_meqtl_covariates_for_cv <- function(train_ids) {
  control_fm <- get_control_fm()
  covariates <- get_covs() %>%
    as.data.frame()
  
  rownames(covariates) <- covariates$SUBJID
  covariates <- model.matrix(control_fm, covariates)[get_meth_ids(), -1][train_ids, ]
  covariates <- t(covariates)
  meqtl_covariates <- SlicedData$new(covariates)
}

run_car_score_estimation <- function(meth_site, snps, model_frame, fm) {
  db <- cbind(y = meth_site, snps, model_frame)
  m <- lm(fm, db)
  calc.relimp(m, type = "car")
}


compute_relative_importance <- function(train_ids, meth_list, snp_list, 
                                        meth_location_list, snp_location_list, model_frame, fm_cells, n_cores) {
  
  cl <- makeCluster(n_cores)
  
  clusterEvalQ(cl = cl, {
    library(splines)
    library(MatrixEQTL)
    library(data.table)
    library(dplyr)
    library(dtplyr)
    source("./Scripts/R_scripts/Libraries/functions_for_meQTL_mapping2.R")
  })
  
  
  clusterExport(cl = cl, c("compute_R2", "run_computations", "str_add_str"))
  
  best_snp_per_probe_per_chr <- clusterMap(cl = cl,
                                           fun = run_local_genotype_mapping_for_chr,
                                           meth_list, snp_list, meth_location_list, snp_location_list,
                                           MoreArgs = list(train_ids = train_ids, meqtl_covariates = setup_meqtl_covariates_for_cv(train_ids)),
                                           RECYCLE = FALSE)
  
  meth <- reduce(meth_list[-1], ~ cbind(.x, .y[-1]), .init = meth_list[[1]][-1])
  
  get_snp_mat <- function(best_snp_per_probe, snps) {
    set_names(as.list(as.data.frame(snps[, best_snp_per_probe$Best_snp, drop = FALSE])), best_snp_per_probe$Probe)
  }
  
  # Get cis SNP for each probe. If there are none introduce a NULL list entry
  snp_list <- unlist(map2(best_snp_per_probe_per_chr, snp_list, get_snp_mat), recursive = FALSE)
  snp_list[names(meth)[!names(meth) %in% names(snp_list)]] <- list(NULL)
  snp_list <- snp_list[names(meth)]
  
  res <- clusterMap(cl = cl,
                    run_computations,
                    meth,
                    snp_list,
                    MoreArgs = list(train_ids = train_ids, model_frame = model_frame, fm_cells = fm_cells))
  
  stopCluster(cl)
  
  bind_rows(res, .id = "Probe") %>% 
    left_join(bind_rows(best_snp_per_probe_per_chr)) %>% 
    select(Probe, Best_snp, everything())
  
}

run_local_genotype_mapping_for_chr <- function(meth, snps, meth_locations, snp_locations, train_ids, meqtl_covariates) {
  
  meth <- meth[train_ids, ]
  snps <- snps[train_ids, ]
  
  id <- meth$SUBJID
  meth$SUBJID <- NULL
  meqtl_meth <- t(meth)
  colnames(meqtl_meth) <- id
  meqtl_meth <-  SlicedData$new(meqtl_meth)
  meqtl_meth$ResliceCombined(sliceSize = 1e2)
  
  meqtl_genotype <- t(snps)
  meqtl_genotype <-  SlicedData$new(meqtl_genotype)
  meqtl_genotype$ResliceCombined(sliceSize = 1e2)
  
  # Call the main analysis function
  Matrix_eQTL_main(snps = meqtl_genotype,
                   gene = meqtl_meth,
                   cvrt = meqtl_covariates,
                   snpspos = as.data.frame(snp_locations[c("SNP_ID", "SNP_chr", "SNP_position")]),
                   genepos = as.data.frame(meth_locations[c("Probe", "Probe_chr", "Probe_position", "Probe_position")]),
                   pvOutputThreshold = 0,
                   pvOutputThreshold.cis = 1,
                   output_file_name = NULL,
                   output_file_name.cis = NULL,
                   cisDist = 5 * 1e4,
                   useModel = modelLINEAR,
                   verbose = FALSE)$cis$eqtls %>%
    lazy_dt() %>%
    rename(Probe = gene) %>% 
    mutate(snps = as.character(snps), Probe = as.character(Probe)) %>% 
    group_by(Probe) %>%
    summarize(Best_snp = snps[which.min(pvalue)]) %>%
    ungroup() %>%
    as_tibble()
  
}

setup_and_run_local_meQTL_mapping <- function(meth_chr, snps_chr, meth_loc_chr, snp_loc_chr, train_ids, mf, fm_cells) {
  
  # Setup meth
  meqtl_meth <- as.data.frame(meth_chr)
  rownames(meqtl_meth) <- meth_chr$SUBJID
  meqtl_meth$SUBJID <- NULL
  meqtl_meth <- meqtl_meth[train_ids, ]
  meqtl_meth <- setup_methylation_for_meQTL_mapping(as.matrix(meqtl_meth))
  
  # Setup genotypes
  meqtl_snps <- setup_genotypes_for_meQTL_mapping(snps_chr[train_ids, ])
  
  # Setup model matrix
  mm <- model.matrix(update(as.formula(fm_cells), ~ Sex + Smoking_status + ns(Age, df = 3) + CMV_serostatus + .), mf)[train_ids, -1]
  meqtl_mm <- setup_covariates_for_meQTL_mapping(mm)
  
  best_snp_per_probe <- run_local_meqtl_for_chr(meqtl_meth, meqtl_snps, meth_loc_chr, snp_loc_chr, mm = meqtl_mm) %>%
    lazy_dt() %>%
    group_by(Probe) %>%
    summarize(Best_snp = SNP[which.min(Pvalue)]) %>%
    ungroup() %>%
    as_tibble()
  
  set_names(best_snp_per_probe$Best_snp, best_snp_per_probe$Probe)
  
}

compute_rel_importance_chr <- function(meth_chr, snps_chr, meth_loc_chr, snp_loc_chr, train_ids, mf, fm_cells) {
  
  best_snp_per_probe <- setup_and_run_local_meQTL_mapping(meth_chr,
                                                          snps_chr = snps_chr,
                                                          meth_loc_chr = meth_loc_chr,
                                                          snp_loc_chr = snp_loc_chr,
                                                          train_ids = train_ids,
                                                          mf = mf,
                                                          fm_cells = fm_cells)
  
  
  
  meth_chr$SUBJID <- NULL
  res <- vector(mode = "list", length = ncol(meth_chr))
  names(res) <- names(meth_chr)
  
 
  
  for (cpg_site in names(meth_chr)) {

    snp_name <- best_snp_per_probe[cpg_site]
    if (!is.na(snp_name)) SNP <- snps_chr[, snp_name]
    else SNP <- NULL
    
    res[[cpg_site]] <- run_computations(meth_chr[[cpg_site]], SNP, train_ids = train_ids, model_frame = mf, fm_cells = fm_cells)
      
  }
  
  bind_rows(res, .id = "Probe") 
  
}


compute_relative_importance_cv <- function(meth_list, snp_list, meth_loc_list, snp_loc_list, mf, fm_cells, n_cores) {
  
  cross_val_ids <- sample.int(5, nrow(mf), replace = TRUE)
  res <- vector("list", 5)
  
  for (fold in 1:5) {
    
    train_ids <- cross_val_ids != fold
    res[[fold]] <- mcmapply(compute_rel_importance_chr, 
                            meth_list, 
                            snp_list, 
                            meth_loc_list, 
                            snp_loc_list,
                            MoreArgs = list(train_ids, mf = mf, fm_cells = fm_cells),
                            SIMPLIFY = FALSE,
                            mc.cores = n_cores) %>% bind_rows()
    
  }
  
  bind_rows(res, .id = "Fold")
  
}


run_computations <- function(meth_site, SNP, train_ids, model_frame, fm_cells) {
  
  fm_IDOL <- y ~ IDOL_CD8T + IDOL_CD4T + IDOL_NK + IDOL_Bcell + IDOL_Mono
  fm_Houseman <- y ~ Houseman_CD8T + Houseman_CD4T + Houseman_NK + Houseman_Bcell + Houseman_Mono
  fm_panel5 <- y ~ Panel5_X_NK_of_total + Panel5_X_mono_of_total + Panel5_X_CD4_of_total + 
    Panel5_X_CD8_of_total + Panel5_X_CD19pos_of_total + Panel5_X_CD4negCD8neg_of_total + Panel5_X_neutrophils_of_total
  
  fm_cells <- update(as.formula(fm_cells), y ~ . )
  
  
  if (!is.null(SNP)) db <- cbind(y = meth_site, SNP = SNP, model_frame)
  else db <- cbind(y = meth_site, model_frame)
  
  db_train <- na.omit(db[train_ids, ])
  db_test <- na.omit(db[!train_ids, ])
  
  
  # Direct effect analysis
  model_list <- list()
  
  # Fit models for direct effect estimation
  model_list$Total_effect_IDOL <-           lm(fm_IDOL, db_train)
  model_list$Total_effect_Houseman <-       lm(fm_Houseman, db_train)
  model_list$Total_effect_cells <-          lm(fm_cells, db_train)
  model_list$Total_effect_Panel5 <-         lm(fm_panel5, db_train)
  model_list$Intrinsic <-                   update(model_list$Total_effect_cells, . ~ . + Sex + ns(Age, df = 3))
  model_list$Exposures <-                   update(model_list$Total_effect_cells, . ~ . + CMV_serostatus + Smoking_status + Log_CRP_levels)
  
  # Fit models for total effect estimation
  model_list$Total_effect_intrinsic <-      lm(y ~ ns(Age, df = 3) + Sex, db_train) 
  model_list$Total_effect_exposures <-      lm(y ~ CMV_serostatus + Smoking_status + Log_CRP_levels, db_train)
  
  
  if (!is.null(SNP)) {
    model_list$Genetic <- update(model_list$Total_effect_cells, . ~ . + SNP)
    model_list$Total_effect_genetic <- lm(y ~ SNP, db_train) 
  }
  
  predictions <- lapply(model_list, function(m) predict(m, db_test))
  obs <- db_test$y
  R2 <- lapply(predictions, function(pred) 100 * compute_R2(obs = obs, pred = pred))
  
  # Predictive accuracy for full model
  fm_full <- get_control_fm_15_cells() %>% 
    update(y ~ . - Ancestry_PC1 - Ancestry_PC2)
  
  if (!is.null(SNP)) fm_full <- update(fm_full, . ~ . + SNP)
  m <- lm(fm_full, db_train)
  pred <- predict(m, newdata = db_test)
  R2_full <- 100 * compute_R2(obs = obs, pred = pred)
  
  # Assemble results
  res <- tibble(Conditional_effect_intrinsic = step(R2$Intrinsic - R2$Total_effect_cells),      
                Conditional_effect_exposures = step(R2$Exposures - R2$Total_effect_cells),
                Total_effect_intrinsic = R2$Total_effect_intrinsic,
                Total_effect_exposures = R2$Total_effect_exposures,
                Total_effect_panel5 = R2$Total_effect_Panel5,
                Total_effect_cells = R2$Total_effect_cells,
                Total_effect_IDOL = R2$Total_effect_IDOL,
                Total_effect_Houseman = R2$Total_effect_Houseman,
                Full_effect = R2_full)
  
  if (!is.null(SNP)) {
    bind_cols(res, tibble(Conditional_effect_genetic = step(R2$Genetic - R2$Total_effect_cells), 
                          Total_effect_genetic = R2$Total_effect_genetic))
  } else {
    bind_cols(res, tibble(Conditional_effect_genetic = 0,
                          Total_effect_genetic = 0))
  }
}