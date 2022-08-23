clean_matrix_eQTL_results <- function(res, meth_locations = get_meth_location(), snp_locations = get_snp_location()) {
  
  res %>%
    as_tibble() %>% 
    mutate(snps = as.character(snps), gene = as.character(gene)) %>% 
    left_join(snp_locations[c("SNP_ID", "SNP_position", "SNP_gene", "SNP_chr")], by = c("snps" = "SNP_ID")) %>% 
    rename(SNP = snps, Probe = gene, Statistic = statistic, Pvalue = pvalue, Mean = beta) %>%
    left_join(meth_locations[c("Probe", "Probe_chr", "Probe_position", "Probe_gene")], by = "Probe")
  
}

setup_methylation_for_meQTL_mapping <- function(meth) {
  meth <- t(meth)
  meth <-  SlicedData$new(meth)
  meth$ResliceCombined(sliceSize = 1e2)
  meth
}

setup_covariates_for_meQTL_mapping <- function(covariates) {
  covariates <- t(covariates)
  SlicedData$new(covariates)
}

setup_genotypes_for_meQTL_mapping <- function(snps) {
  snps <- t(snps)
  meqtl_snps <-  SlicedData$new(snps)
  meqtl_snps$ResliceCombined(sliceSize = 1e2)
  meqtl_snps
}


fit_lm <- function(SNP, db) {
  db$x <- SNP
  db <- na.omit(db)
  m_null <- lm(y ~ . - x, db)
  m_alt <- update(m_null, . ~ . + x)
  p_val <- anova(m_null, m_alt, test = "F")[2, "Pr(>F)"]
  if (is.na(p_val)) p_val <- 1
  p_val
}

test_conditional_association <- function(db, SNPs) {
  apply(SNPs, 2, fit_lm, db = db)
}

filter_conditionally <- function(tib, keys, meth = meth, geno = geno, p_thresh) {
  
  # Screen for cis SNP first 
  
  
  cis_SNP <- tib$Cis_SNP[1]
  cpg_site <- meth[[keys$Probe]]
  
  if (!is.na(cis_SNP)) {  
    
    db <- data.frame(y = cpg_site, Cis_SNP = geno[, cis_SNP])
    p_vals <- test_conditional_association(db, geno[, tib$SNP, drop = FALSE])
    tib <- tib %>% filter(!SNP %in% names(p_vals)[p_vals > p_thresh])  
    
  }
  
  if (nrow(tib) > 1) {
    
    cnt <- 1
    
    
    while (nrow(tib) > cnt) {
      
      if (!is.na(cis_SNP)) {
        
        db <- as.data.frame(geno[, c(cis_SNP, tib$SNP[1:cnt]), drop = FALSE])  
        names(db) <- c("cis_SNP", paste0("SNP_cond", 1:cnt))
        
      } else {
        
        db <- as.data.frame(geno[, tib$SNP[1:cnt], drop = FALSE])  
        names(db) <- paste0("SNP_cond", 1:cnt)
        
      }
      
      db$y <- cpg_site
      p_vals <- test_conditional_association(db,
                                             geno[, tib$SNP[-(1:cnt)], drop = FALSE])
      tib <- tib %>% filter(!SNP %in% names(p_vals)[p_vals > p_thresh])
      cnt <- cnt + 1
      
    }
    
    if (nrow(tib) != 1) {
      
      
      db <- as.data.frame(geno[, tib$SNP])
      names(db) <- tib$SNP
      db$y <- cpg_site
      m <- lm(y ~ ., db)
      
      keep_snps <- tidy(m) %>% 
        filter(term != "(Intercept)", p.value < p_thresh) %>% 
        pull(term)
      
      tib <- tib %>% filter(SNP %in% keep_snps)
      
    }
  
  }
  
  tib
  
}

filter_for_ld <- function(tib, keys, R_thresh = 0.5, geno) {
  
  cnt <- 1
  cor_with_top <- 1
  
  while (nrow(tib) > cnt) {
    
    cor_with_top <- cor(geno[, tib$SNP[cnt]],
                        geno[, tib$SNP[-(1:cnt)], drop = FALSE],
                        use = "pairwise.complete.obs")
    tib <- tib %>% filter(!SNP %in% colnames(cor_with_top)[abs(cor_with_top) > R_thresh])
    cnt <- cnt + 1
    
  }
  
  tib
  
}

filter_snps_conditionally <- function(tib, meth, geno, p_tresh = 1e-6) {
  
  tib %>%
    group_by(Probe, SNP_chr) %>%
    arrange(desc(abs(Statistic)), SNP) %>%
    group_modify( ~ filter_conditionally(.x, .y, meth = meth, geno = geno, p_thresh = p_tresh)) %>% 
    ungroup()
  
}

run_local_meqtl_for_chr <- function(meth, snps, meth_loc, snp_loc, mm) {
  
  res <- Matrix_eQTL_main(snps = snps,
                          gene = meth,
                          cvrt = mm,
                          snpspos = as.data.frame(snp_loc[c("SNP_ID", "SNP_chr", "SNP_position")]),
                          genepos = as.data.frame(meth_loc[c("Probe", "Probe_chr", "Probe_position", "Probe_position")]),
                          pvOutputThreshold = 0,
                          pvOutputThreshold.cis = 1,
                          output_file_name = NULL,
                          output_file_name.cis = NULL,
                          cisDist = 5 * 1e4,
                          useModel = modelLINEAR,
                          verbose = FALSE)$cis$eqtls
  
  clean_matrix_eQTL_results(res, meth_loc, snp_loc)
  
}

# filter_nr_probe_tib <- function(tib, geno, R_thresh = 0.5) {
#   tib %>%
#     group_by(SNP_chr) %>%  -
#     arrange(NR_associated_probes) %>%
#     group_modify( ~ filter_for_ld(.x, .y, R_thresh = R_thresh, geno = geno))
# }

setup_meth_for_local_meQTL_mapping_chr <- function(chr, type = "M_values") {
  meth_chr <- readRDS(glue::glue("./Data/Chunk_data/Methylation/{type}/Per_chromosome/meth_chromosome_{chr}.rds"))
  rownames(meth_chr) <- meth_chr$SUBJID
  meth_chr$SUBJID <- NULL
  setup_methylation_for_meQTL_mapping(as.matrix(meth_chr))
}

setup_genotypes_for_local_meQTL_mapping_chr <- function(chr) {
  snps <- readRDS(glue::glue("./Data/Chunk_data/Genotype/Chromosome_{chr}/LabExMI_imputation_958x5699237_chr{chr}_snp_matrix.rds"))
  setup_genotypes_for_meQTL_mapping(snps)
}

setup_meth_loc_chr <- function(chr) {
  meth_loc_chr <- meth_locs[meth_locs$Probe_chr == chr, ]
}

setup_snp_loc_chr <- function(chr) {
  readRDS(glue::glue("./Data/Chunk_data/Genotype/Chromosome_{chr}/LabExMI_imputation_1000x5699237_chr{chr}_map.rds"))
}
