library(glue)
library(tidyverse)
library(lme4)
library(splines)
library(lmtest)
library(lmerTest)
library(parallel)
first <- dplyr::first
select <- dplyr::select
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")

covs <- get_covs_884()
m_values <- get_m_values()
beta_values <- get_beta_values()

beta_values <- beta_values[match(covs$SUBJID, beta_values$SUBJID),]
m_values <- m_values[match(covs$SUBJID, m_values$SUBJID),]



meth_loc <- get_anno_meth_location()
meth_roadmap <- get_anno_meth_roadmap()
chromatin_state_enrichments <- readRDS("./Data/RData/Results/EWAS_chromatin_state_enrichments_bonferroni_884.rds")
TFBS_enrichments_dir <- readRDS("./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_dir_EWAS_884_bonferroni.rds")
TFBS_enrichments_med <- readRDS("./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_med_EWAS_884_bonferroni.rds")


# MeQTL -------------------------------------------------------------------

# TFBS_direct <- readRDS("./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_dir_EWAS.rds")
# TFBS_med <- readRDS("./Data/RData/Results/TFBS_enrichments/TFBS_enrichment_med_EWAS.rds")
# 
# chromatin_state_enrichments <- readRDS("./Tables/chromatin_state_enrichments.rds")
# geographic_region_enrichments <- readRDS("./Tables/geographic_region_enrichments.rds")
# TFBS_enrichments <- readRDS("./Tables/")
#   
# eqtlgen <- readRDS("./Data/RData/Genotypes/eqtlgen.rds") %>% 
#   select(SNP_ID, EQTLgen_GeneSymbol, SNP_position, SNP_chr, EQTLgen_GenePos, EQTLgen_GeneChr, EQTLgen_Pvalue) %>% 
#   mutate(EQTLgen_GeneChr = as.integer(EQTLgen_GeneChr))
# 
# map_labex <- readRDS("./Data/RData/Genotypes/LabExMI_imputation_1000x5699237_annotated_map.rds")
# trans_probes <- readRDS("./Data/RData/Methylation/probes_trans_50K.rds")
# meth_loc <- get_anno_meth_location()
# meth_roadmap <- get_anno_meth_roadmap()
# meth_anno_TFBS <- readRDS("./Data/RData/Methylation/Annotation/meth_anno_TFBS.rds") %>% 
#   filter(Probe %in% trans_probes)
# 
# TF_list <- names(select(meth_anno_TFBS, -Probe))
# 
# meth_anno_TFBS <- bind_cols(meth_anno_TFBS["Probe"], 
#                             keep(select(meth_anno_TFBS, -Probe), ~ sum(. == "Yes") > 20))
# 
# trans_adjust_on_snps <- readRDS("./Data/RData/Results/MeQTL/trans_snp_bottom_layer_sign.rds") %>%
#   mutate(Mean = -2 * Mean) %>%
#   left_join(select(map_labex, SNP_ID, Minor_allele, Major_allele, SNP_distance_to_gene_bp), by = c("SNP" = "SNP_ID"))
# 

# cis_meqtl <- readRDS("./Data/RData/Results/MeQTL/cis_adjusted_probes.rds")
# sum(cis_meqtl$Significant_family)
# mean(cis_meqtl$Significant_family)
# 
# cis_meqtl_beta <- readRDS("./Data/RData/Results/MeQTL/cis_beta_probes.rds") %>% 
#   mutate(Mean = -2 * Mean)
# 
# cis_meqtl_beta %>% arrange(desc(abs(Mean))) %>% slice(1:1000) %>% arrange(abs(Mean))
# cis_meqtl_beta %>% filter(abs(Mean) > 0.3)
# 
# trans_adjust_on_snps_independent <- readRDS("./Data/RData/Results/MeQTL/trans_snp_bottom_layer_sign_independent.rds")
# trans_meqtl_beta <- readRDS("./Data/RData/Results/MeQTL/trans_beta_probes.rds") %>% 
#   mutate(Mean = -2 * Mean)
# 
# trans_meqtl_beta %>% arrange(desc(abs(Mean))) %>% slice(1:100) %>% arrange(abs(Mean))
# trans_meqtl_beta %>% filter(abs(Mean) > 0.15)
# 
# n_distinct(trans_adjust_on_snps_independent$SNP)
# n_distinct(trans_adjust_on_snps_independent$Probe)
# 
# wide_effect_meqtls <- trans_adjust_on_snps %>%
#   group_by(SNP) %>%
#   summarize(SNP_gene = first(SNP_gene),
#             Chromosome = first(SNP_chr),
#             Position = first(SNP_position),
#             NR_associated_probes = n(),
#             NR_chr = n_distinct(Probe_chr)) %>%
#   arrange(desc(NR_associated_probes)) %>%
#   group_by(SNP_gene) %>%
#   summarize(Chromosome = first(Chromosome),
#             Gene_pos = median(Position),
#             SNP_pos = Position[which.max(NR_associated_probes)],
#             SNP = SNP[which.max(NR_associated_probes)],
#             NR_associated_probes = max(NR_associated_probes),
#             NR_chr = NR_chr[which.max(NR_associated_probes)]) %>% 
#   left_join(select(trans_adjust_on_snps, SNP, Probe, Probe_chr, Estimate = Mean, Probe_position))
# 
# also_cis_eqtl <- select(wide_effect_meqtls, SNP_gene, SNP) %>% 
#   left_join(select(eqtlgen, SNP = SNP_ID, EQTLgen_GeneSymbol)) %>% 
#   group_by(SNP_gene) %>% 
#   group_modify(~ tibble(Local_eQTL = .y %in% .x$EQTLgen_GeneSymbol))
# 
# wide_effect_meqtls_eqtl_gen <- left_join(wide_effect_meqtls, also_cis_eqtl)
# 
# wide_effect_meqtls_eqtl_gen %>% distinct(SNP_gene, .keep_all = TRUE) %>% pull(Local_eQTL) %>% mean()
# 
# wide_effect_meqtls %>% filter(NR_associated_probes >= 10) %>% distinct(SNP_gene) %>% filter(!grepl("^MSH5-SAPCD1|^FLJ26245|^DQ580909|^ACTR3B|^AX|^AK|^BC|^ANKR|^LOC|^LINC", SNP_gene)) %>% pull(SNP_gene) %>% length()
# 



# CMV ---------------------------------------------------------------------

ewas_med_CMV <- readRDS("./Data/Chunk_data/Results/Cell_mediation/m_values_884/CMV.rds") %>%
  filter(Effect == "Mediation") %>% 
  left_join(meth_loc) %>% 
  left_join(meth_roadmap) %>% 
  mutate(P_bonferroni = p.adjust(P_value))


ewas_med_CMV_beta <- readRDS("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/CMV.rds") %>%
  filter(Effect == "Mediation") %>% 
  left_join(meth_loc) %>% 
  left_join(meth_roadmap) %>% 
  mutate(P_bonferroni = p.adjust(P_value))

0.0209 - 1.96 * 0.00155
0.0209 + 1.96 * 0.00155


ewas_CMV <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/CMV_serostatus.rds") %>% 
  left_join(meth_loc) %>% 
  left_join(meth_roadmap)


ewas_med_CMV %>% filter(Probe == "cg15843262")  

# cg16477774
fit_example_model_15_cells("cg16477774", beta_values, covs)


plt_frame_EMRA <- select(covs, SUBJID, CMV_serostatus, X_CD8_EMRA_of_total.panel1, X_CD4_EMRA_of_total.panel1, Age, Sex, Smoking_status, Day_of_sampling) %>% 
  mutate(Age = Age / 50)

plt_frame_EMRA <- add_column(beta = beta_values[["cg15843262"]], plt_frame_EMRA)



coef(summary(lmer(X_CD4_EMRA_of_total.panel1 ~ CMV_serostatus + Age + Sex + Smoking_status + (1 | Day_of_sampling), plt_frame_EMRA)))
coef(summary(lmer(X_CD8_EMRA_of_total.panel1 ~ CMV_serostatus + Age + Sex + Smoking_status + (1 | Day_of_sampling), plt_frame_EMRA)))


coef(summary(lmer(beta ~ X_CD4_EMRA_of_total.panel1 +  Age + Sex + CMV_serostatus + Smoking_status + (1 | Day_of_sampling), plt_frame_EMRA)))
coef(summary(lmer(beta ~ X_CD8_EMRA_of_total.panel1 +  Age + Sex + CMV_serostatus + Smoking_status + (1 | Day_of_sampling), plt_frame_EMRA)))



# AGE ---------------------------------------------------------------------

ewas_med_age <- readRDS("./Data/Chunk_data/Results/Cell_mediation/m_values_884/Age.rds") %>%
  filter(Effect == "Mediation") %>% 
  left_join(meth_loc) %>% 
  left_join(meth_roadmap) %>% 
  mutate(P_bonferroni = p.adjust(P_value))

ewas_age <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Age.rds") %>% 
  left_join(meth_loc) %>% 
  left_join(meth_roadmap)


disp <- readRDS("./Data/RData/Results/Dispersion/res_dispersion_gamlss_884.rds") %>%
  mutate(P_bonferroni = p.adjust(P_value)) %>%
  left_join(meth_roadmap) %>%
  mutate(Is_outlier = is_outlier(Estimate))

disp_no_cells <- readRDS("./Data/RData/Results/Dispersion/res_dispersion_gamlss_no_cells_884.rds") %>%
  mutate(P_bonferroni = p.adjust(P_value)) %>%
  left_join(meth_roadmap) %>%
  mutate(Is_outlier = is_outlier(Estimate))


disp_adj_var <- readRDS("./Data/RData/Results/Dispersion/res_dispersion_gamlss_variance_adjusted_for_cells_884.rds") %>%
  mutate(P_bonferroni = p.adjust(P_value)) %>%
  left_join(meth_roadmap) %>%
  mutate(Is_outlier = is_outlier(Estimate))


# Sex ---------------------------------------------------------------------


ewas_med_sex <- readRDS("./Data/Chunk_data/Results/Cell_mediation/m_values_884/Sex.rds") %>%
  filter(Effect == "Mediation") %>% 
  left_join(meth_loc) %>% 
  left_join(meth_roadmap) %>% 
  mutate(P_bonferroni = p.adjust(P_value))

ewas_sex <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Sex.rds") %>% 
  left_join(meth_loc) %>% 
  left_join(meth_roadmap)

frame <- beta_values[c("SUBJID", "cg09516963")] %>% inner_join(select(covs, Sex, SUBJID))


# Smoking ---------------------------------------------------------------------


ewas_med_smoking <- readRDS("./Data/Chunk_data/Results/Cell_mediation/m_values_884/Smoking.rds") %>%
  filter(Effect == "Mediation") %>% 
  left_join(meth_loc) %>% 
  left_join(meth_roadmap) %>% 
  mutate(P_bonferroni = p.adjust(P_value))

ewas_med_smoking_beta <- readRDS("./Data/Chunk_data/Results/Cell_mediation/Beta_values_884/Smoking.rds") %>%
  filter(Effect == "Mediation") %>% 
  left_join(meth_loc) %>% 
  left_join(meth_roadmap) %>% 
  mutate(P_bonferroni = p.adjust(P_value))

ewas_smoking <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Smoking_status.rds") %>% 
  left_join(meth_loc) %>% 
  left_join(meth_roadmap)

frame <- beta_values[c("SUBJID", "cg09516963")] %>% inner_join(select(covs, Smoking, SUBJID))

m <- fit_example_model_15_cells("cg05575921", beta_values, covs)
100 * confint(m)["Smoking_statusSmoker", ]
m <- fit_example_model_15_cells("cg03636183", beta_values, covs)
100 * confint(m)["Smoking_statusSmoker", ]
m <- fit_example_model_15_cells("cg17739917", beta_values, covs)
100 * confint(m)["Smoking_statusSmoker", ]


probes <- ewas_smoking %>% 
  filter(P_bonferroni <= 0.05) %>% 
  pull(Probe) %>% 
  unique()

years_since_last_smoke <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Years_since_last_smoke.rds")
years_smoking <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Years_smoking.rds") %>% 
  filter(Probe %in% probes)

plt_frame_years_since_last_smoke <- years_since_last_smoke %>% 
  right_join(tibble(Probe = probes))

plt_frame_years_smoking <- years_smoking %>% 
  right_join(tibble(Probe = probes))

plt_frame <- tibble(Probe = probes, 
                    Years_since_last_smoke = plt_frame_years_since_last_smoke$Estimate, 
                    Years_smoking = plt_frame_years_smoking$Estimate)
lm(plt_frame$Years_smoking ~ plt_frame$Years_since_last_smoke)

plt_frame_NK <- select(covs, SUBJID, CMV_serostatus, X_NK_of_total.panel4, Age, Sex, Smoking_status, Day_of_sampling) %>% 
  mutate(Age = Age / 50)

plt_frame_NK <- add_column(beta = beta_values[["cg11335172"]], plt_frame_NK)

coef(summary(lmer(X_NK_of_total.panel4 ~ CMV_serostatus + Age + Sex + Smoking_status + (1 | Day_of_sampling), plt_frame_NK)))
coef(summary(lmer(beta ~ X_NK_of_total.panel4 +  Age + Sex + CMV_serostatus + Smoking_status + (1 | Day_of_sampling), plt_frame_NK)))

#

years_since_last_smoke_med <- readRDS("./Data/Chunk_data/Results/Cell_mediation/m_values/Years_since_last_smoke.rds") %>% filter(Effect == "Mediation")
years_smoking_med <- readRDS("./Data/Chunk_data/Results/Cell_mediation/m_values/Years_smoking.rds") %>% filter(Effect == "Mediation")

probes <- ewas_med_smoking %>% 
  filter(P_bonferroni <= 0.05) %>% 
  pull(Probe) %>% 
  unique()

plt_frame_years_since_last_smoke_med <- years_since_last_smoke_med %>% 
  right_join(tibble(Probe = probes))

plt_frame_years_smoking_med <- years_smoking_med %>% 
  right_join(tibble(Probe = probes))

plt_frame_med <- tibble(Probe = probes, 
                        Years_since_last_smoke_med = plt_frame_years_since_last_smoke_med$Estimate, 
                        Years_smoking_med = plt_frame_years_smoking_med$Estimate)



# CRP and other factors

ewas_CRP_total <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_no_cells/Log_CRP_levels.rds") %>% 
  left_join(meth_loc) %>% 
  left_join(meth_roadmap)

ewas_CRP <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Log_CRP_levels.rds") %>% 
  left_join(meth_loc) %>% 
  left_join(meth_roadmap)


ewas_trig <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Log_triglyceride_levels.rds") %>% 
  left_join(meth_loc) %>% 
  left_join(meth_roadmap)


ewas_HDL <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Log_HDL_levels.rds") %>% 
  left_join(meth_loc) %>% 
  left_join(meth_roadmap)

ewas_log_uric_acid <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_no_genotypes/Correct_for_16_props/Log_uric_acid_levels.rds") %>% 
  left_join(meth_loc) %>% 
  left_join(meth_roadmap)

ewas_raw_fruit <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Raw_fruits_consumption.rds") %>% 
  left_join(meth_loc) %>% 
  left_join(meth_roadmap)



f <- function(CpG, meth, covs, snp_mat = NULL, linear_age = FALSE) {
  
  fm <- get_control_fm_15_cells() %>% 
    update(as.formula(paste0(CpG, " ~ . + Log_CRP_levels "))) %>% 
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


m <- f("cg11099984", beta_values, covs)
100 * confint(m)["Log_CRP_levels", ] 


m <- f("cg06500161", beta_values, covs)
100 * confint(m)["Log_CRP_levels", ] 



f <- function(CpG, meth, covs, snp_mat = NULL, linear_age = FALSE) {
  
  fm <- get_control_fm_15_cells() %>% 
    update(as.formula(paste0(CpG, " ~ . + Log_triglyceride_levels"))) %>% 
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


m <- f("cg06500161", beta_values, covs)
100 * confint(m)["Log_triglyceride_levels", ] 


f <- function(CpG, meth, covs, snp_mat = NULL, linear_age = FALSE) {
  
  fm <- get_control_fm_15_cells() %>% 
    update(as.formula(paste0(CpG, " ~ . + Log_HDL_levels"))) %>% 
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


m <- f("cg06500161", beta_values, covs)
100 * confint(m)["Log_HDL_levels", ] 


f <- function(CpG, meth, covs, snp_mat = NULL, linear_age = FALSE) {
  
  fm <- get_control_fm_15_cells() %>% 
    update(as.formula(paste0(CpG, " ~ . + Log_HDL_levels + Log_CRP_levels + Log_triglyceride_levels"))) %>% 
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


m <- f("cg06500161", beta_values, covs)






f <- function(CpG, meth, covs, snp_mat = NULL, linear_age = FALSE) {
  
  fm <- get_control_fm_15_cells() %>% 
    update(as.formula(paste0(CpG, " ~ . + Raw_fruits_consumption"))) %>% 
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


m <- f("cg23726427", beta_values, covs)
100 * confint(m)["Raw_fruits_consumptionTwice_A_Day_Or_More", ]





# Genotype interaction ----------------------------------------------------

fit_example_decomposition_model <- function(cpg_site, snps_site, mf) {
  
  db <- cbind(cpg_site = cpg_site, snps_site, mf)
  fm <- glue('cpg_site ~ ns(Age, df = 3) + Sex + Smoking_status + CMV_serostatus + ({paste0(colnames(snps_site), collapse = " + ")}) * Myeloid')
  lm(fm, db)
  
}

fit_example_interaction_model <- function(meth_site, snps_site, mf) {
  
  db <- cbind(y = meth_site, mf, snps_site)
  
  fm_snps <- paste0(colnames(snps_site), collapse = " + ")
  fm_interaction_terms <- glue('(Age + Sex + CMV_serostatus + Smoking_status) * ({fm_snps})')
  fm_interaction <- glue('y ~ Ancestry_PC1 + Ancestry_PC2 + {fm_cells} + {fm_interaction_terms}')
  
  lm(fm_interaction, db)
}

# Myeloid proportion
n_cores <- 10
lineage_cells <- get_panel5_cells()

covs <- get_covs_884() %>% 
  mutate(Age = Age / 50) %>%
  dplyr::select(SUBJID, Age, Sex, CMV_serostatus, Smoking_status, X_dendritic_cells.panel8) %>% 
  inner_join(lineage_cells, by = "SUBJID") %>% 
  mutate(Lymphoid = X_CD19pos_of_total.panel5 + X_NK_of_total.panel5 + X_CD8_of_total.panel5 + X_CD4_of_total.panel5 + X_CD4negCD8neg_of_total.panel5,
         Myeloid = X_PMN_of_total.panel5 + X_mono_of_total.panel5)

snps <- readRDS("./Data/RData/Genotypes/significant_snp_mat_per_probe.rds") %>% compact()
snps <- mclapply(snps, function(x) x[covs$SUBJID, , drop = FALSE], mc.cores = n_cores)

meth <- get_beta_values()
meth <- meth[match(covs$SUBJID, meth$SUBJID),]
meth <- meth[names(snps)]

mf <- select(covs, Age, Sex, Smoking_status, CMV_serostatus, Myeloid, X_dendritic_cells.panel8)

cis_meqtl <- readRDS("./Data/RData/Results/MeQTL/cis_adjusted_probes_884_bonferroni.rds")
trans_meqtl <- readRDS("./Data/RData/Results/MeQTL/trans_884_bonferroni_sign_independent.rds")
meqtl_chromatin_enrichment <- readRDS("./Data/RData/Results/MeQTL/meQTL_chromatin_state_enrichments_884_bonferroni.rds")

res_interaction <- readRDS("./Tables/interaction_table.rds")

fit_example_decomposition_model(meth[["cg07195891"]], snps[["cg07195891"]], mf)

db <- cbind(cpg_site = meth[["cg07195891"]], snps[["cg07195891"]], mf)
fm <- glue('cpg_site ~ ns(Age, df = 3) + Sex + Smoking_status + CMV_serostatus + ({paste0(colnames(snps[["cg07195891"]]), collapse = " + ")}) * X_dendritic_cells.panel8')
lm(fm, db)


# Interaction
cell_list_15_cells <- get_15_props()
fm_cells <- paste0(cell_list_15_cells, collapse = " + ")

mf <- get_covs_884() %>% 
  dplyr::select(SUBJID, all_of(cell_list_15_cells), Sex, Smoking_status, Log_CRP_levels, Age, CMV_serostatus, Ancestry_PC1, Ancestry_PC2) %>% 
  mutate(Age = Age / 50, Smoking_status = fct_collapse(Smoking_status, Non_Smoker = c("Non_Smoker", "Ex_Smoker"))) %>%
  as.data.frame()

m <- fit_example_interaction_model(meth[["cg21268422"]], snps[["cg21268422"]], mf)
m <- fit_example_interaction_model(meth[["cg09064148"]], snps[["cg09064148"]], mf)


# Proportion of variance explained ----------------------------------------

vars <- c("Total_effect_intrinsic", "Total_effect_exposures", "Total_effect_cells", "Total_effect_genetic")
prop_var <- readRDS("./Data/RData/Results/Proportion_of_variance_explained/prop_var_explained_genetic_interactions_averaged_884.rds")

total_effect <- prop_var %>% 
  select(-contains("Houseman"), -contains("IDOL")) %>% 
  select(Probe | starts_with("Total") | Full_effect) %>%
  pivot_longer(!Probe, names_to = "Predictor", values_to = "R2") %>% 
  filter(Predictor != "Total_effect_panel5") %>% 
  mutate(Predictor = fct_recode(Predictor,
                                `Full model` = "Full_effect",
                                Intrinsic = "Total_effect_intrinsic", 
                                Exposures = "Total_effect_exposures", 
                                `16 Cells` = "Total_effect_cells",
                                Genetic = "Total_effect_genetic")) %>% 
  group_by(Predictor) %>% 
  mutate(Is_outlier = is_outlier(R2)) %>% 
  ungroup()

total_effect_probe <- total_effect %>% 
  filter(!Predictor %in% c("6 Cells", "Full model")) %>% 
  group_by(Probe) %>% 
  filter(R2 == R2[which.max(R2)]) %>% 
  ungroup()


direct_effect <- prop_var %>%
  select(Probe | starts_with("Conditional")) %>%
  pivot_longer(!Probe, names_to = "Predictor", values_to = "R2") %>%
  mutate(Predictor = fct_recode(Predictor, 
                                Intrinsic = "Conditional_effect_intrinsic", 
                                Exposures = "Conditional_effect_exposures", 
                                Genetic = "Conditional_effect_genetic",
                                `Genetic int.` = "Conditional_effect_genetic_interaction")) %>% 
  group_by(Predictor) %>% 
  mutate(Is_outlier = is_outlier(R2)) %>% 
  ungroup()

direct_effect_probe <- direct_effect %>%
  group_by(Probe) %>%
  summarize(Top_pred = Predictor[which.max(R2)],
            R2_sum = sum(R2)) %>% 
  ungroup()


# OLD Proportion of variance explained ----------------------------------------

vars <- c("Total_effect_intrinsic", "Total_effect_exposures", "Total_effect_cells", "Total_effect_genetic")
prop_var <- readRDS("./Data/RData/Results/Proportion_of_variance_explained/prop_var_explained_genetic_interactions_averaged_884.rds")

total_effect <- prop_var %>% 
  select(-contains("Houseman"), -contains("IDOL")) %>% 
  select(Probe | starts_with("Total") | Full_effect) %>%
  pivot_longer(!Probe, names_to = "Predictor", values_to = "R2") %>% 
  filter(Predictor != "Total_effect_panel5") %>% 
  mutate(Predictor = fct_recode(Predictor,
                                `Full model` = "Full_effect",
                                Intrinsic = "Total_effect_intrinsic", 
                                Exposures = "Total_effect_exposures", 
                                `16 Cells` = "Total_effect_cells",
                                Genetic = "Total_effect_genetic")) %>% 
  group_by(Predictor) %>% 
  mutate(Is_outlier = is_outlier(R2)) %>% 
  ungroup()

total_effect_probe <- total_effect %>% 
  filter(!Predictor %in% c("6 Cells", "Full model")) %>% 
  group_by(Probe) %>% 
  filter(R2 == R2[which.max(R2)]) %>% 
  ungroup()


direct_effect <- prop_var %>%
  select(Probe | starts_with("Conditional")) %>%
  pivot_longer(!Probe, names_to = "Predictor", values_to = "R2") %>%
  mutate(Predictor = fct_recode(Predictor, 
                                Intrinsic = "Conditional_effect_intrinsic", 
                                Exposures = "Conditional_effect_exposures", 
                                Genetic = "Conditional_effect_genetic",
                                `Genetic int.` = "Conditional_effect_genetic_interaction")) %>% 
  group_by(Predictor) %>% 
  mutate(Is_outlier = is_outlier(R2)) %>% 
  ungroup()

direct_effect_probe <- direct_effect %>%
  group_by(Probe) %>%
  summarize(Top_pred = Predictor[which.max(R2)],
            R2_sum = sum(R2)) %>% 
  ungroup()

p <- total_effect %>% filter(Predictor == "Full model") 
mean(p$R2 < 5)
sum(p$R2 < 5)

p <- total_effect %>% filter(Predictor == "Full model") 
mean(p$R2 > 25)
sum(p$R2 > 25)

p <- total_effect %>% filter(Predictor == "Full model", R2 > 25) %>% pull(Probe)
p <- total_effect_probe %>% filter(Probe %in% p)
p <- table(p$Predictor)
100 * (p / sum(p))


p <- total_effect %>% filter(Predictor == "16 Cells")
mean(p$R2 > 25)
sum(p$R2 > 25)


p <- total_effect %>% filter(Predictor == "Full model", R2 > 75) %>% pull(Probe)
p <- total_effect_probe %>% filter(Probe %in% p)
p <- table(p$Predictor)
100 * (p / sum(p))


p <- readRDS("./Data/RData/Results/Proportion_of_variance_explained/prop_var_explained_genetic_interactions_P_value_884.rds") %>% 
  mutate(P_bonferroni = p.adjust(P))

p <- p %>% filter(P_bonferroni < 0.05) %>% pull(Probe)
direct_effect %>% filter(Predictor == "Genetic int.", Probe %in% p, R2 > 5)


