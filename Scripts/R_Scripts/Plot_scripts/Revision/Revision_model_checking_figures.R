

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(glue)
library(parallel)
library(corrr)
library(tidyverse)
library(sandwich)
library(lmtest)
library(broom)
library(lme4)
library(splines)
library(patchwork)
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
options(stringsAsFactors = FALSE)
options(width = 200)

# Main --------------------------------------------------------------------

n_cores <- 8
cell_list_15_cells <- get_15_props()
fm_cells <- paste0(cell_list_15_cells, collapse = " + ")

mf <- get_covs_884() %>% 
  dplyr::select(SUBJID, all_of(cell_list_15_cells), Sex, Smoking_status, Log_CRP_levels, Age, CMV_serostatus, Ancestry_PC1, Ancestry_PC2, Day_of_sampling) %>% 
  mutate(Age = Age / 50)

snps <- readRDS("./Data/RData/Genotypes/significant_snp_mat_per_probe.rds")
snps <- purrr::discard(snps, is_null)

snps <- mclapply(snps, function(x) x[mf$SUBJID, , drop = FALSE], mc.cores = n_cores)

meth <- get_m_values()
meth <- meth[match(mf$SUBJID, meth$SUBJID),]


meth_anno_location <- get_anno_meth_location()
ewas_CMV <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/CMV_serostatus.rds")
ewas_CMV <- left_join(ewas_CMV, meth_anno_location) %>% arrange(P_value)

ewas_age <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Age.rds")
ewas_age <- left_join(ewas_age, meth_anno_location) %>% arrange(P_value)

ewas_sex <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/Correct_for_16_props/Sex.rds")
ewas_sex <- left_join(ewas_age, meth_anno_location) %>% arrange(P_value)

cpg_sites <- c(ewas_CMV$Probe[1:3],
               ewas_age$Probe[1:3],
               ewas_sex$Probe[1:3])

models <- map(cpg_sites, ~ fit_example_model_15_cells(., meth, mf, snps[[.]]))
res <- lapply(models, function(x) tibble(Residuals = residuals(x), Age = mf$Age, CMV = mf$CMV_serostatus, Sex = mf$Sex, Smoking = mf$Smoking_status)) %>% 
  bind_rows(.id = "Model")


pltA <- res %>%
  ggplot(aes(sample = Residuals)) + 
  stat_qq() + 
  stat_qq_line() +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  facet_wrap(~ Model, scales = 'free', nrow = 5, ncol = 3)



pltB <- res %>%
  ggplot(aes(x = 50 * Age, y = Residuals)) +
  geom_point() +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  facet_wrap(~ Model, scales = 'free', nrow = 5, ncol = 3)



pltC <- res %>%
  ggplot(aes(x = CMV, y = Residuals)) +
  geom_violin() +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  facet_wrap(~ Model, scales = 'free', nrow = 5, ncol = 3)

pltD <- res %>%
  ggplot(aes(x = Sex, y = Residuals)) +
  geom_violin() +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  facet_wrap(~ Model, scales = 'free', nrow = 5, ncol = 3)


# ASsemble plot -----------------------------------------------------------



