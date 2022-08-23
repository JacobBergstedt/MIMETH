library(tidyverse)
library(glue)
library(broom)
library(fs)
library(vroom)
library(ggplot2)
library(gridExtra)
library(grid)
library(patchwork)
library(GGally)
library(latex2exp)
library(scales)
library(matrixStats)
library(ggrepel)
source("./Scripts/R_scripts/Plot_scripts/locuszoom.R")
source("./Scripts/R_scripts/Libraries/functions_for_plotting.R")
source("./Scripts/R_scripts/Libraries/functions_for_odds_ratios.R")

setup_plt_frame <- function(db, CpG, genotype, genotype_map) {

  plt_frame <- db
  plt_frame$Genotype <- snps_for_plotting[, snp_name]
  snp_map <- genotype_map %>% filter(SNP_ID == snp_name)
  plt_frame$Genotype <- factor(genotype)
  levels(plt_frame$Genotype) <- c(paste0(snp_map$Minor_allele, snp_map$Minor_allele),
                                  paste0(snp_map$Major_allele, snp_map$Minor_allele),
                                  paste0(snp_map$Major_allele, snp_map$Major_allele))
  plt_frame$Genotype <- fct_rev(plt_frame$Genotype)
  plt_frame$CpG <- CpG
  plt_frame

}


# Load data -----------------------------------------------
# meth <- readRDS("./Data/RData/Methylation/MIMETH.minfi.final.betaMatrix.autosomes.no_outliers.rds")
genotype_map <- readRDS("./Data/RData/Genotypes/LabExMI_imputation_1000x5699237_annotated_map_with_ancestral.rds") %>%
  mutate(Which_allele_is_derived = if_else(Which_allele_is_ancestral == "Minor", "Major", "Minor"))

trans_adjust_on_snps_independent <- readRDS("./Data/RData/Results/MeQTL/trans_884_bonferroni_sign_independent.rds") %>%
  left_join(select(genotype_map, SNP_ID, Minor_allele, Major_allele, SNP_distance_to_gene_bp, Which_allele_is_derived), by = c("SNP" = "SNP_ID")) %>%
  mutate(Mean = case_when(
    Which_allele_is_derived == "Minor" ~ - 2 * Mean,
    Which_allele_is_derived == "Major" ~ 2 * Mean,
    is.na(Which_allele_is_derived) ~ -2 * Mean
  ))

# roadmap_desc <- read_tsv("./Data/Roadmap_15_states/roadmap_descriptions.tsv")
# original_roadmap_levels <- paste0(roadmap_desc$STATE_NU, "_", roadmap_desc$MNEMONIC)

meth_anno_roadmap <- get_anno_meth_roadmap()
meth_anno_loc <- get_anno_meth_location()
meth_anno_geo <- get_anno_meth_geography()

cis_adjusted_probes <- readRDS("./Data/RData/Results/MeQTL/cis_adjusted_probes_884_bonferroni.rds") %>%
  filter(Locally_significant) %>%
  left_join(select(genotype_map, SNP_ID, Which_allele_is_derived), by = c("SNP" = "SNP_ID")) %>%
  mutate(Mean = case_when(
    Which_allele_is_derived == "Minor" ~ - 2 * Mean,
    Which_allele_is_derived == "Major" ~ 2 * Mean,
    is.na(Which_allele_is_derived) ~ -2 * Mean
  ))


chromatin_state_enrichment <- readRDS("./Data/RData/Results/MeQTL/meQTL_chromatin_state_enrichments_884_bonferroni.rds")
geography_enrichments <- readRDS("./Data/RData/Results/MeQTL/meQTL_geography_enrichments_884_bonferroni.rds")


# Load data ---------------------------------------------------------------

res_int <- readRDS("./Data/RData/Results/genotype_interactions_sign_SNPs_884.rds")
count_sign <- res_int %>% dplyr::count(Variable, .drop = FALSE)

covs <- get_covs_884()
lineage_cells <- get_panel5_cells()

mf <- covs %>%
  select(-X_mono_of_total.panel5) %>%
  inner_join(lineage_cells, by = "SUBJID") %>%
  mutate(Myeloid = X_mono_of_total.panel5 + X_neutrophils_of_total.panel5,
         Lymphoid = X_CD19pos_of_total.panel5 + X_NK_of_total.panel5 + X_CD8_of_total.panel5 + X_CD4_of_total.panel5 + X_CD4negCD8neg_of_total.panel5) %>%
  select(SUBJID, Age, Sex, Smoking_status, CMV_serostatus, Myeloid, Lymphoid, Log_CRP_levels)

snps_for_plotting <- readRDS("./Data/RData/Genotypes/snps_for_genotype_interaction_plots.rds")[mf$SUBJID, ]
meth <- get_beta_values()
meth <- meth[match(mf$SUBJID, meth$SUBJID),]


# Plot parameters -----------------------------------------------
# state_colors <- c("red", "orangered", "green3", "forestgreen", "darkgreen", "yellowgreen", "goldenrod1",
#                   "turquoise", "lightslateblue", "coral3", "darksalmon", "khaki3", "gray70", "gray85", "white")
# 
cis_col <- "#74C476"
trans_col <- "#006D2C"
one_digit <- function(x) sprintf("%.1f", x)


# 4 A ---------------------------------------------------------------------

plt_frame <- chromatin_state_enrichment %>% 
  filter(Variable == "Local_CpG") %>% 
  rename(Term = State)


pltA <- plot_odds_ratio(plt_frame, 
                        NULL, 
                        bar_size = odds_ratio_errorbar_size, 
                        color = cis_col, 
                        base_size = font_size_odds_ratio,
                        limits = c(0, 4),
                        breaks = seq(0, 4, 1 / 2)) +
  labs(tag = "A")


ggsave("./Plots/Revision_plots/Figure4/Fig4A_local_meQTL_CpG_chromatin_state_enrichments.pdf", pltA, width = 6, height = 4, units = "cm")


# 4 B ---------------------------------------------------------------------

plt_frame <- chromatin_state_enrichment %>% 
  filter(Variable == "Remote_CpG") %>% 
  rename(Term = State)


pltB <- plot_odds_ratio(plt_frame, 
                        NULL, 
                        bar_size = odds_ratio_errorbar_size, 
                        color = trans_col, 
                        base_size = font_size_odds_ratio) +
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())
  labs(tag = "B")


ggsave("./Plots/Revision_plots/Figure4/Fig4B_remote_meQTL_CpG_chromatin_state_enrichments.pdf", pltB, width = 6, height = 4, units = "cm")


# 4 C ---------------------------------------------------------------------


plt_frame <- chromatin_state_enrichment %>% 
  filter(Variable == "Remote_SNP") %>% 
  rename(Term = State)


pltC <- plot_odds_ratio(plt_frame, 
                        NULL, 
                        bar_size = odds_ratio_errorbar_size, 
                        color = trans_col, 
                        base_size = font_size_odds_ratio,
                        limits = c(0, 30),
                        breaks = seq(0, 30, 2)) +
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(tag = "C")


ggsave("./Plots/Revision_plots/Figure4/Fig4C_remote_meQTL_SNP_chromatin_state_enrichments.pdf", pltC, width = 6, height = 4, units = "cm")

# 4 D ---------------------------------------------------------------------

pltDa <- res_int %>%
  ggplot(aes(x = Variable, y = -log10(P_FDR), color = Variable)) +
  geom_jitter(width = 0.2, size = point_size_scatter) +
  scale_color_manual(values = rev(c(col_age, col_sex, col_smoking, col_CMV, col_CRP, col_cell)), drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  ylab("-log10(Adjusted Pvalue)") +
  xlab(NULL) +
  coord_flip() +
  theme_classic(base_size = font_size_scatter) +
  theme(legend.position = "none",
        axis.text = element_text(color = "black")) +
  labs(tag = "D")

pltDb <- count_sign %>% 
  ggplot(aes(x = Variable, y = 0, label = n, fill = Variable)) +
  geom_tile() + 
  geom_text(color = "white", size = font_size_table * (1 / ggplot2:::.pt)) +
  scale_fill_manual(values = rev(c(col_age, col_sex, col_smoking, col_CMV, col_CRP, col_cell))) +
  coord_flip() +
  xlab(NULL) + 
  ylab(NULL) +
  theme_classic(base_size = font_size_scatter) +
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

pltD <- pltDa + pltDb + plot_layout(widths = c(6, 1))
ggsave("./Plots/Revision_plots/Figure4/Fig4D_interaction_p_values.pdf", pltD, width = 6, height = 6, units = "cm")


# 4 E ---------------------------------------------------------------------

cpg_site <- "cg21268422"
snp_name <- "rs2837990"
plt_frame <- setup_plt_frame(mf, CpG = meth[[cpg_site]], genotype = snps_for_plotting[, snp_name], genotype_map = genotype_map)
pltE <- plt_frame %>% 
  ggplot(aes(x = Age, y = CpG, color = Genotype)) +
  geom_point(size = point_size_scatter) +
  geom_smooth(method = "lm", se = FALSE, size = line_size_scatter) +
  ylab(expression(paste("5mC (", italic("BACE2"), ")"))) +  
  scale_color_manual(values = rev(greens[c(3, 5, 7)]), guide = guide_legend(override.aes = list(size = point_size_scatter))) +
  theme_bw(base_size = font_size_scatter) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.title = element_text(size = font_size_scatter), 
        legend.text  = element_text(size = font_size_scatter),
        legend.key.size = unit(0.2, "lines"),
        legend.position = c(0.85, 0.85)) +
  labs(tag = "E")

ggsave("./Plots/Revision_plots/Figure4/Fig4E_BACE2.pdf", pltE, width = 6, height = 6, units = "cm")

# 4 F ---------------------------------------------------------------------

cpg_site <- "cg07195891"
snp_name <- "rs11055602"
plt_frame <- setup_plt_frame(mf, CpG = meth[[cpg_site]], genotype = snps_for_plotting[, snp_name], genotype_map = genotype_map) %>% 
  filter(!is.na(Genotype))

pltF <- plt_frame %>% 
  ggplot(aes(x = Lymphoid, y = CpG, color = Genotype)) +
  geom_point(size = point_size_scatter) +
  geom_smooth(method = "lm", se = FALSE, size = line_size_scatter) +
  ylab(expression(paste("5mC (", italic("CLEC4C"), ")"))) +  
  scale_color_manual(values = rev(greens[c(3, 5, 7)]), guide = guide_legend(override.aes = list(size = point_size_scatter))) +
  theme_bw(base_size = font_size_scatter) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.title = element_text(size = font_size_scatter), 
        legend.text  = element_text(size = font_size_scatter),
        legend.key.size = unit(0.2, "lines"),
        legend.position = c(0.85, 0.85)) +
  labs(tag = "F")

ggsave("./Plots/Revision_plots/Figure4/Fig4F_CLEC4C.pdf", pltF, width = 6, height = 6, units = "cm")



# Assemble ----------------------------------------------------------------


design <- "
  ABC
  DEF
"

plt <- pltA + pltB + pltC + pltD + pltE + pltF + plot_layout(design = design) &
  theme(plot.tag = element_text(size = 6, face = "bold"))

ggsave("./Plots/Revision_plots/Figure4/Figure_4.png", plt, width = 172, height = 172, units = "mm", dpi = 500)
ggsave("./Plots/Revision_plots/Figure4/Figure_4.pdf", plt, width = 172, height = 172, units = "mm")


# Plot all panels ------------------------------------------------
# p <- b + c + d + e + f +
# plot_layout(ncol = 5, widths = c(4,1,4,4,1))

# ggsave("~/Figure_1.pdf", p, width = 18, height = 5, units = "cm")
# (a + b + c + d + e + f) / g

