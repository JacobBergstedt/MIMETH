save_beta_values_per_chr <- function(chr, beta_values, meth_locations) {

  probes <- meth_locations %>% 
    filter(Probe_chr == chr) %>% 
    pull(Probe)
  
  saveRDS(cbind(beta_values["SUBJID"], beta_values[probes]), 
          paste0("./Data/Chunk_data/Methylation/Beta/Per_chromosome/meth_chromosome_", chr, ".rds"))
}

save_m_values_per_chr <- function(chr, m_values, meth_locations) {
  
  probes <- meth_locations %>% 
    filter(Probe_chr == chr) %>% 
    pull(Probe)
  
  saveRDS(cbind(m_values["SUBJID"], m_values[probes]), 
          paste0("./Data/Chunk_data/Methylation/M_values/Per_chromosome/meth_chromosome_", chr, ".rds"))
}

save_cis_snps_for_probes <- function(n_cores, meth_locations, snp_locations) {
  save_cis_snps_per_chr <- function(chr, probe_split, meth_locations, snp_locations) {
    probe_split_chr <- probe_split[[chr]]
    meth_locations_chr <- filter(meth_locations, Probe_chr == chr)
    snp_locations_chr <- filter(snp_locations, SNP_chr == chr)
    snps_chr <- get_snp_matrix(chr) 
    for (chunk_nr in seq_len(length(probe_split_chr))) {
      
      probes <- probe_split_chr[[chunk_nr]]
      probe_loc <- right_join(meth_locations_chr, tibble(Probe = probes), by = "Probe")
      probe_snps <- map(set_names(probe_loc$Probe_position, probe_loc$Probe), 
                        ~ snp_locations_chr$SNP_ID[abs(. - snp_locations_chr$SNP_position) < 5e4])
      snps_chunk <- snps_chr[, unique(unlist(probe_snps))]
      saveRDS(snps_chunk, paste0("./Data/Chunk_data/Genotype/Chromosome_", chr, "/Chunks/snps_chunk_", chunk_nr, ".rds"))
      saveRDS(probe_snps, paste0("./Data/Chunk_data/Genotype/Chromosome_", chr, "/Chunks/probe_snps_chunk_", chunk_nr, ".rds"))  
    }
  }
  probe_split <- split_probes(meth_locations)
  mclapply(1:22, save_cis_snps_per_chr, probe_split = probe_split, meth_locations = meth_locations, snp_locations = snp_locations, mc.cores = n_cores)
}

divide_length_into_groups <- function(nr_groups, len) {
  splits <- kronecker(1:nr_groups, rep(1, floor(len / nr_groups)))
  c(splits, rep(max(splits), len - length(splits)))  
}

split_probes <- function(meth_locations, nr_groups = 20) {
  meth_locations %>% 
    group_by(Probe_chr) %>% 
    mutate(Chunk = divide_length_into_groups(nr_groups, n())) %>% 
    group_map(~ group_by(tibble(Chunk = .$Chunk, Probe = .$Probe), Chunk)) %>% 
    map(~ set_names(group_map(., ~ .$Probe), seq_len(nr_groups))) %>% 
    set_names(1:22)
}


save_probe_chunks <- function(meth, meth_locations, meth_type = "Beta", n_cores = 12) {
  save_methylation_chunks_per_chr <- function(probe_chunks_chr, chr, meth, meth_type) {
    for (chunk_nr in seq_along(probe_chunks_chr)) {
      meth_chunk <- meth[c("SUBJID", probe_chunks_chr[[chunk_nr]])]
      dir_path <- paste0("./Data/Chunk_data/Methylation/", meth_type, "/Per_chromosome/Chromosome_", chr, "/Chunks/")
      saveRDS(meth_chunk, paste0(dir_path, "Chunk_", chunk_nr, ".rds"))
    }
  }
  probe_split <- split_probes(meth_locations)
  mclapply(1:22,
      function(chr) save_methylation_chunks_per_chr(probe_chunks_chr = probe_split[[chr]], chr = chr, meth = meth, meth_type = meth_type), mc.cores = 12)
}

save_genotypes_per_chr <- function(chr, meth_locations, snp_locations) {
  meth_locations_chr <- filter(meth_locations, Probe_chr == chr)
  snp_locations_chr <- filter(snp_locations, SNP_chr == chr)
  snps_chr <- snps[, snp_locations_chr$SNP_ID]
  path <- paste0("./Data/Chunk_data/Genotype/Chromosome_", chr, "/")
  saveRDS(snps_chr, paste0(path, "LabExMI_imputation_958x5699237_chr", chr, "_snp_matrix.rds"))
  saveRDS(snp_locations_chr, paste0(path, "LabExMI_imputation_1000x5699237_chr", chr, "_map.rds"))
}