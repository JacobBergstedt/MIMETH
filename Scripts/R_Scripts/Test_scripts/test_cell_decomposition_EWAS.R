library(tidyverse)

CD8_rat <- readRDS("./Data/RData/Results/EWAS/M_values/Environment/Sandwich/cell_decomposition_ewas_CD8_naive_ratio.rds") %>% 
  filter(Variable == "Age")

CD4_rat <- readRDS("./Data/RData/Results/EWAS/M_values/Environment/Sandwich/cell_decomposition_ewas_CD4_naive_ratio.rds") %>% 
  filter(Variable == "Age")

lymph_rat <- readRDS("./Data/RData/Results/EWAS/M_values/Environment/Sandwich/cell_decomposition_ewas_ratio_884.rds") %>% 
  filter(Variable == "Age")

geno_rat <- readRDS("./Data/RData/Results/decomposition_lymphoid_genotype_sign_SNPs_m_values_ratio_884.rds") %>% 
  select(Probe, SNP, Rat = Estimate, Rat_P = P_bonferroni)

geno <- readRDS("./Data/RData/Results/decomposition_lymphoid_genotype_sign_SNPs_884.rds") %>% 
  select(Probe, SNP, Prop = Estimate, Prop_P = P_bonferroni)

CD8 <- readRDS("./Data/RData/Results/EWAS/Beta_values/Environment/Sandwich/cell_decomposition_ewas_CD8_naive.rds") %>% 
  filter(Variable == "Age")

CD4 <- readRDS("./Data/RData/Results/EWAS/Beta_values/Environment/Sandwich/cell_decomposition_ewas_CD4_naive.rds") %>% 
  filter(Variable == "Age")

lymph <- readRDS("./Data/RData/Results/EWAS/Beta_values/Environment/Sandwich/cell_decomposition_EWAS_lymphoid.rds") %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  filter(Variable == "Age")

p <- select(CD8_rat, Probe, CD8_rat = Estimate, CD8_rat_P = P_bonferroni) %>% 
  inner_join(select(CD4_rat, Probe, CD4_rat = Estimate, CD4_rat_P = P_bonferroni)) %>% 
  inner_join(select(lymph_rat, Probe, lymph_rat = Estimate, lymph_rat_P = P_bonferroni)) %>% 
  inner_join(select(CD8, Probe, CD8 = Estimate, CD8_P = P_bonferroni)) %>% 
  inner_join(select(CD4, Probe, CD4 = Estimate, CD4_P = P_bonferroni)) %>% 
  inner_join(select(lymph, Probe, lymph = Estimate, lymph_P = P_bonferroni))

geno_frame <- inner_join(geno, geno_rat)
