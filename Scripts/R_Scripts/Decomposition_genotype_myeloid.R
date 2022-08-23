

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(glue)
library(parallel)
library(tidyverse)
library(splines)
library(sandwich)
library(lmtest)
library(broom)
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
options(stringsAsFactors = FALSE)
options(width = 200)

fit_example_decomposition_model <- function(cpg_site, snps_site, mf) {
  
  db <- cbind(cpg_site = cpg_site, snps_site, mf)
  fm <- glue('cpg_site ~ ns(Age, df = 3) + Sex + Smoking_status + CMV_serostatus + ({paste0(colnames(snps_site), collapse = " + ")}) * Myeloid')
  lm(fm, db)
  
}

fit_decomposition_model <- function(cpg_site, snps_site, mf) {
  
  
  db <- cbind(cpg_site = cpg_site, snps_site, mf)
  fm <- glue('cpg_site ~ ns(Age, df = 3) + Sex + Smoking_status + CMV_serostatus + ({paste0(colnames(snps_site), collapse = " + ")}) * Myeloid')
  m <- lm(fm, db)
  terms <- expand.grid(colnames(snps_site), "Myeloid", stringsAsFactors = FALSE)
  terms_int <- paste0(terms$Var1, ":", terms$Var2)
 
  sum_m <- coeftest(m, vcov. = vcovHC)
  keep_ident <- terms_int %in% rownames(sum_m)
  terms_int <- terms_int[keep_ident]
  terms <- terms[keep_ident, ]
  
  tidy(sum_m) %>% 
    filter(term %in% terms_int) %>%
    rename(Estimate = estimate, Standard_error = std.error, P_value = p.value) %>% 
    add_column(SNP = terms$Var1, Variable = terms$Var2, .before = 1)
  
}

# Main --------------------------------------------------------------------

n_cores <- 12
lineage_cells <- get_panel5_cells()

covs <- get_covs_884() %>% 
  mutate(Age = Age / 50) %>%
  dplyr::select(SUBJID, Age, Sex, CMV_serostatus, Smoking_status) %>% 
  inner_join(lineage_cells, by = "SUBJID") %>% 
  mutate(Lymphoid = X_CD19pos_of_total.panel5 + X_NK_of_total.panel5 + X_CD8_of_total.panel5 + X_CD4_of_total.panel5 + X_CD4negCD8neg_of_total.panel5,
         Myeloid = X_PMN_of_total.panel5 + X_mono_of_total.panel5)

snps <- readRDS("./Data/RData/Genotypes/significant_snp_mat_per_probe.rds") %>% compact()
snps <- mclapply(snps, function(x) x[covs$SUBJID, , drop = FALSE], mc.cores = n_cores)
  

# meth <- get_m_values()
meth <- get_beta_values()
meth <- meth[match(covs$SUBJID, meth$SUBJID),]
meth <- meth[names(snps)]


mf <- select(covs, Age, Sex, Smoking_status, CMV_serostatus, Myeloid)

res <- mcmapply(fit_decomposition_model, 
                meth, 
                snps, 
                MoreArgs = list(mf = mf), 
                SIMPLIFY = FALSE, 
                mc.cores = n_cores) %>% 
  bind_rows(.id = "Probe") %>% 
  select(Probe, SNP, Variable, Estimate, Standard_error, P_value) %>% 
  mutate(P_FDR = p.adjust(P_value, "BH"), 
         P_bonferroni = p.adjust(P_value))

saveRDS(res, "./Data/RData/Results/decomposition_myeloid_genotype_sign_SNPs_beta_values_884.rds")
# saveRDS(res, "./Data/RData/Results/decomposition_lymphoid_genotype_sign_SNPs_m_values_884.rds")
