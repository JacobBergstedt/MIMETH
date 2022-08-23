

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
ewas_sex <- left_join(ewas_sex, meth_anno_location) %>% arrange(P_value)

cpg_sites <- c(ewas_CMV$Probe[1:3],
               ewas_age$Probe[1:3],
               ewas_sex$Probe[1:3])

models <- set_names(cpg_sites, cpg_sites) %>% map(~ fit_example_model_15_cells(., meth, mf, snps[[.]]))
res <- lapply(models, function(x) tibble(Residuals = residuals(x), Age = mf$Age, CMV = mf$CMV_serostatus, Sex = mf$Sex, Smoking = mf$Smoking_status)) %>% 
  bind_rows(.id = "Model")

# Create plots ------------------------------------------------------------

pltA <- res %>%
  ggplot(aes(sample = Residuals)) + 
  stat_qq(size = point_size_scatter) + 
  stat_qq_line() +
  theme_bw(base_size = font_size_scatter) +
  theme(axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        strip.background = element_blank()) +
  facet_wrap(~ Model, scales = 'free', nrow = 3, ncol = 3) +
  labs(tag = "A")



pltB <- res %>%
  ggplot(aes(x = 50 * Age, y = Residuals)) +
  geom_point(size = point_size_scatter, col = col_age) +
  theme_bw(base_size = font_size_scatter) +
  theme(axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        strip.background = element_blank()) +
  facet_wrap(~ Model, scales = 'free', nrow = 3, ncol = 3) +
  labs(tag = "B")



pltC <- res %>%
  ggplot(aes(x = CMV, y = Residuals)) +
  geom_violin(fill = col_CMV) +
  theme_bw(base_size = font_size_scatter) +
  theme(axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        strip.background = element_blank()) +
  facet_wrap(~ Model, scales = 'free', nrow = 3, ncol = 3) +
  labs(tag = "C")

pltD <- res %>%
  ggplot(aes(x = Sex, y = Residuals)) +
  geom_violin(fill = col_sex) +
  theme_bw(base_size = font_size_scatter) +
  theme(axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        strip.background = element_blank()) +
  facet_wrap(~ Model, scales = 'free', nrow = 3, ncol = 3) +
  labs(tag = "D")


# Assemble plot -----------------------------------------------------------

design <- "
  AB
  CD
"

plt <- pltA + pltB + pltC + pltD + plot_layout(design = design) &
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave("./Plots/Revision_plots/Supp_fig4/Sup_fig_4.png", plt, width = 180, height = 180, units = "mm", dpi = 500)
ggsave("./Plots/Revision_plots/Supp_fig4/Sup_fig_4.pdf", plt, width = 180, height = 180, units = "mm")


