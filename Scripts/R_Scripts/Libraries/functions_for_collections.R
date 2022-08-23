
# Local functions ---------------------------------------------------------

read_and_filter <- function(path, threshold) {
  filter(readRDS(path), Pvalue < threshold)
}

read_and_find <- function(path, snp, probe) {
  readRDS(path) %>% 
    filter(SNP == snp, Probe == probe)
}

collect_meQTL_results <- function(paths, n_cores, threshold = 1, add_annotation) {
  if (threshold > 0.99) res <- bind_rows(mclapply(paths, readRDS, mc.cores = n_cores))
  else res <- bind_rows(mclapply(paths, read_and_filter, threshold = threshold, mc.cores = n_cores))
  
  if (add_annotation) add_snp_meth_locations(res, get_meth_location(), get_snp_location())
  else res
}

# Collection functions ----------------------------------------------------

collect_trans_meQTLs <- function(n_cores, threshold = 1e-4, type = "M_values", add_annotation = FALSE) {
  paths <- paste0("./Data/Chunk_data/Results/MeQTL/Trans_meQTL/", type, "/chromosome_", 1:22, ".rds")
  collect_meQTL_results(paths, n_cores, threshold, add_annotation)
}

collect_cell_decomposition_cis_meQTL_results <- function(n_cores, threshold = 1e-4, add_annotation = FALSE) {
  grid <- expand.grid(1:22, 1:20)
  paths <- paste0("./Data/Chunk_data/Results/MeQTL/Cell_decomposition_cis_meQTL/Chromosome_", grid$Var1, "/res_chunk_", grid$Var2, ".rds")
  collect_meQTL_results(paths, n_cores, threshold, add_annotation)
}

collect_interaction_local_meQTL_results <- function(n_cores, threshold = 1e-4, add_annotation = FALSE) {
  grid <- expand.grid(c("Age", "Sex", "CMV_serostatus", "Smoking_status"), 1:22)
  paths <- paste0("./Data/Chunk_data/Results/MeQTL/Interaction_meQTL/Local/", grid$Var1, "/chromosome_", grid$Var2, ".rds")
  collect_meQTL_results(paths, n_cores, threshold, add_annotation)
}

collect_interaction_trans_meQTL_results <- function(n_cores, threshold = 1e-4, add_annotation = FALSE) {
  grid <- expand.grid(c("Age", "Sex", "CMV_serostatus"), 1:22)
  paths <- paste0("./Data/Chunk_data/Results/MeQTL/Interaction_meQTL/Trans/", grid$Var1, "/chromosome_", grid$Var2, ".rds")
  collect_meQTL_results(paths, n_cores, threshold, add_annotation)
}

collect_ewas_results <- function(chunk_path, n_cores) {
  paths <- paste0(chunk_path, list.files(chunk_path))
  ewas <- bind_rows(mclapply(paths, readRDS, mc.cores = n_cores))
}

collect_ewas_age_interaction_results <- function(n_cores = 12) {
  chunk_path <- "./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect/Age_interactions/"
  paths <- paste0(chunk_path, list.files(chunk_path))
  bind_rows(mclapply(paths, readRDS, mc.cores = n_cores))
}

collect_variance_explained_results <- function(n_cores = 12) {
  paths <- paste0("./Data/Chunk_data/Results/Variance_explained/Non_normalized_cells_958/res_rep_", 1:4, ".rds")
  res <- mclapply(paths, readRDS, mc.cores = n_cores) %>% 
    bind_rows(.id = "Rep")
}

collect_variance_explained_results_884 <- function(n_cores = 12) {
  
  paths <- paste0("./Data/Chunk_data/Results/Variance_explained/Normalized_cells_884/Beta/res_rep_", 1:4, ".rds")
  res <- mclapply(paths, readRDS, mc.cores = n_cores) %>% 
    bind_rows(.id = "Rep")
  
}
