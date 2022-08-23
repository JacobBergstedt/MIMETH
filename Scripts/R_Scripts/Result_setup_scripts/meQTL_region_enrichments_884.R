library(tidyverse)
library(glue)
library(broom)
library(scales)
library(patchwork)
source("./Scripts/R_scripts/Libraries/functions_for_plotting.R")

run_state_enrichment <- function(tib, is_case) {
  
  cell_line_list <- paste0("Chromatin_states_", cell_lines_sub)
  res <- vector(mode = "list", length = length(cell_line_list))
  names(res) <- cell_line_list
  
  for (i in 1:length(cell_line_list)) {
    
    cell_states <- tib[[cell_line_list[[i]]]]
    res[[i]] <- map(set_names(states, states), ~ cell_states == .) %>% 
      map(function(is_state) xtabs(~ is_state + is_case)) %>% 
      map(fisher.test) %>% 
      map_dfr(tidy, .id = "State") %>%
      dplyr::rename(Estimate = estimate, Low = conf.low, High = conf.high, P_value = p.value)
    
  }
  
  bind_rows(res, .id = "Cell") %>%
    mutate(State = factor(State, states)) %>% 
    mutate(Cell = str_sub(Cell, start = 18)) %>%
    mutate(Cell = gsub("_", " ", Cell), Cell = gsub(" cells", "", Cell)) %>% 
    mutate(Cell = factor(Cell, cell_lines_sub_label)) %>% 
    mutate(State = factor(State, states))
  
}

run_geography_enrichment <- function(tib, is_case) {
  
  map(set_names(geographic_regions, geographic_regions), ~ tib$Geography_labels == .) %>%
    map(function(is_region) xtabs(~ is_region + is_case)) %>% 
    map(fisher.test) %>% 
    map_dfr(tidy, .id = "Geographic_region") %>%
    dplyr::rename(Estimate = estimate, Low = conf.low, High = conf.high, P_value = p.value)
  
}

prop_overlap <- function(tib,  meth_anno_TFBS) {
  
  prop_overlap_for_TF <- function(sites_TF, target, trans_probes) {
    
    mean(target %in% trans_probes[sites_TF == "Yes"])
    
  }
  
  target <- tib$Probe
  trans_probes <- meth_anno_TFBS$Probe
  meth_anno_TFBS["Probe"] <- NULL
  props <- map_dbl(meth_anno_TFBS, prop_overlap_for_TF, target = target, trans_probes = trans_probes)
  tibble(TF = names(props), Prop_overlap = props)
  
}

# -------------------------------------------------------------------------

meth_loc <- get_anno_meth_location(FALSE)
meth_roadmap <- get_anno_meth_roadmap_all_cells()
meth_anno_TFBS <- readRDS("./Data/RData/Methylation/Annotation/meth_anno_TFBS.rds")
meth_geography <- get_anno_meth_geography()

eqtlgen <- readRDS("./Data/RData/Genotypes/eqtlgen.rds") %>% 
  select(-EQTLgen_Gene) %>% 
  mutate(EQTLgen_GeneChr = as.integer(EQTLgen_GeneChr))

map_labex <- readRDS("./Data/RData/Genotypes/LabExMI_imputation_1000x5699237_annotated_map.rds")

trans_adjust_on_snps <- readRDS("./Data/RData/Results/MeQTL/trans_884_bonferroni_sign.rds") %>% 
  left_join(select(map_labex, SNP_ID, Minor_allele, Major_allele, SNP_distance_to_gene_bp), by = c("SNP" = "SNP_ID"))


trans_adjust_on_snps_independent <- readRDS("./Data/RData/Results/MeQTL/trans_884_bonferroni_sign_independent.rds") %>% 
  left_join(select(map_labex, SNP_ID, Minor_allele, Major_allele, SNP_distance_to_gene_bp), by = c("SNP" = "SNP_ID"))

# trans_adjust_on_snps <- readRDS("./Data/RData/Results/MeQTL/trans_snp_bottom_layer_sign.rds") %>% 
#   left_join(select(map_labex, SNP_ID, Minor_allele, Major_allele, SNP_distance_to_gene_bp), by = c("SNP" = "SNP_ID"))
# 
# trans_adjust_on_snps_independent <- readRDS("./Data/RData/Results/MeQTL/trans_snp_bottom_layer_sign_independent.rds") %>% 
#   left_join(select(map_labex, SNP_ID, Minor_allele, Major_allele, SNP_distance_to_gene_bp), by = c("SNP" = "SNP_ID"))

states <- names(get_anno_roadmap_translation())
roadmap_keys <- get_anno_roadmap_translation()
geographic_regions <- levels(meth_geography$Geography_labels)
cis_meqtl <- readRDS("./Data/RData/Results/MeQTL/cis_adjusted_probes_884_bonferroni.rds")

cell_lines <- get_cell_list_roadmap()
cell_lines_sub <- cell_lines[c(1, 2, 3, 4, 6, 8, 9)]
cell_lines_sub_label <- gsub("_", " ", gsub("_cells", "", cell_lines_sub))

# Run enrichments ---------------------------------------------------------

trans_probes <- readRDS("./Data/RData/Methylation/probes_trans_50K.rds")
meth_roadmap_trans_probes <- meth_roadmap %>% filter(Probe %in% trans_probes)
meth_geography_trans_probes <- meth_geography %>% filter(Probe %in% trans_probes)

is_case <- meth_roadmap_trans_probes$Probe %in% trans_adjust_on_snps_independent$Probe
enrichments_trans_probes <- run_state_enrichment(meth_roadmap_trans_probes, is_case)

is_case <- meth_geography_trans_probes$Probe %in% trans_adjust_on_snps_independent$Probe
geography_enrichments_trans_probes <- run_geography_enrichment(meth_geography_trans_probes, is_case = is_case)

is_case <- map_labex$SNP_ID %in% trans_adjust_on_snps_independent$SNP
enrichments_trans_snps <- run_state_enrichment(select(map_labex, starts_with("Chromatin_states")), is_case)

meth_roadmap_cis <- meth_roadmap %>% filter(Probe %in% cis_meqtl$Probe)
meth_geography_cis <- meth_geography %>% filter(Probe %in% cis_meqtl$Probe)

cis_sign <- cis_meqtl %>% filter(Significant_family)
is_case <- meth_roadmap_cis$Probe %in% cis_sign$Probe
enrichments_cis_probes <- run_state_enrichment(meth_roadmap_cis, is_case)

is_case <- meth_geography_cis$Probe %in% cis_sign$Probe
geography_enrichments_cis_probes <- run_geography_enrichment(meth_geography_cis, is_case)

geography_enrichments <- bind_rows(add_column(Variable = "Local", geography_enrichments_cis_probes, .before = 1),
                                   add_column(Variable = "Remote", geography_enrichments_trans_probes, .before = 1))

# chromatin_state_enrichments_trans_CpG <- enrichments_trans_probes %>% filter(Cell == "CD4 naive") %>% mutate(Variable = "Remote_CpG", .before = 1)
# chromatin_state_enrichments_trans_SNP <- enrichments_trans_snps %>% filter(Cell == "CD4 naive") %>% mutate(Variable = "Remote_SNP", .before = 1)
# chromatin_state_enrichments_cis <- enrichments_cis_probes %>% filter(Cell == "CD4 naive") %>% mutate(Variable = "Local_CpG", .before = 1)
# chromatin_state_enrichments <- bind_rows(chromatin_state_enrichments_trans_CpG,
#                                          chromatin_state_enrichments_trans_SNP,
#                                          chromatin_state_enrichments_cis) %>% select(Variable, State, Estimate, P_value, Low, High)

chromatin_state_enrichments_trans_CpG <- enrichments_trans_probes %>% 
  filter(Cell == "Mononuclear") %>% 
  mutate(Variable = "Remote_CpG", .before = 1)

chromatin_state_enrichments_trans_SNP <- enrichments_trans_snps %>% 
  filter(Cell == "Mononuclear") %>% 
  mutate(Variable = "Remote_SNP", .before = 1)

chromatin_state_enrichments_cis <- enrichments_cis_probes %>% 
  filter(Cell == "Mononuclear") %>% 
  mutate(Variable = "Local_CpG", .before = 1)

chromatin_state_enrichments <- bind_rows(chromatin_state_enrichments_trans_CpG,
                                         chromatin_state_enrichments_trans_SNP,
                                         chromatin_state_enrichments_cis) %>% 
  select(Variable, State, Estimate, P_value, Low, High)


saveRDS(chromatin_state_enrichments, "./Data/RData/Results/MeQTL/meQTL_chromatin_state_enrichments_884_bonferroni.rds")
saveRDS(geography_enrichments, "./Data/RData/Results/MeQTL/meQTL_geography_enrichments_884_bonferroni.rds")


# State enrichment for SNPs -----------------------------------------------

plt <- enrichments_trans_snps %>%
  mutate(State = fct_rev(State)) %>%
  ggplot(aes(x = State, y = Estimate, ymin = Low, ymax = High)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = 1) +
  geom_errorbar(width = 0.5) +
  scale_y_continuous(trans = "log10") +
  coord_flip() +
  xlab(NULL) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_rect(fill = "white")) +
  facet_wrap(vars(Cell), nrow = 1)
ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/trans_snps_meqtl.pdf", plt, width = 20, height = 10, units = "cm")

# State enrichment for CpG sites trans ------------------------------------------

plt <- enrichments_trans_probes %>%
  mutate(State = fct_rev(State)) %>%
  ggplot(aes(x = State, y = Estimate, ymin = Low, ymax = High)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = 1) +
  geom_errorbar(width = 0.5) +
  coord_flip() +
  xlab(NULL) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_rect(fill = "white")) +
  facet_wrap(vars(Cell), nrow = 1)
ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/trans_probes_meqtl.pdf", plt, width = 20, height = 10, units = "cm")

# State enrichment for CpG sites cis ------------------------------------------

plt <- enrichments_cis_probes %>%
  mutate(State = fct_rev(State)) %>%
  ggplot(aes(x = State, y = Estimate, ymin = Low, ymax = High)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = 1) +
  geom_errorbar(width = 0.5) +
  scale_y_continuous(trans = "log10") +
  coord_flip() +
  xlab(NULL) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_rect(fill = "white")) +
  facet_wrap(vars(Cell), nrow = 1)
ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/cis_meqtl.pdf", plt, width = 20, height = 10, units = "cm")

# Plot odds ratios --------------------------------------------------------

plt <- enrichments_cis_probes %>%
  filter(Cell == "Mononuclear") %>%
  rename(Term = State) %>%
  plot_odds_ratio2(title = "Local meQTL CpG sites",
                   bar_size = odds_ratio_errorbar_size,
                   color = col_cis,
                   base_size = font_size_odds_ratio)

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Cis_meQTL/Mononuclear.pdf", plt, width = 4, height = 4, units = "cm")


plt <- enrichments_cis_probes %>%
  filter(Cell == "CD8 naive") %>%
  rename(Term = State) %>%
  plot_odds_ratio2(title = "Local meQTL CpG sites",
                   bar_size = odds_ratio_errorbar_size,
                   color = col_cis,
                   base_size = font_size_odds_ratio)

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Cis_meQTL/CD8_naive.pdf", plt, width = 4, height = 4, units = "cm")

plt <- enrichments_cis_probes %>%
  filter(Cell == "CD4 naive") %>%
  rename(Term = State) %>%
  plot_odds_ratio2(title = "Local meQTL CpG sites",
                   bar_size = odds_ratio_errorbar_size,
                   color = col_cis,
                   base_size = font_size_odds_ratio)

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Cis_meQTL/CD4_naive.pdf", plt, width = 4, height = 4, units = "cm")


plt <- enrichments_trans_probes %>%
  filter(Cell == "Mononuclear") %>%
  rename(Term = State) %>%
  plot_odds_ratio2(title = "Long-range meQTL CpGs",
                   bar_size = odds_ratio_errorbar_size,
                   color = col_trans,
                   base_size = font_size_odds_ratio)

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Trans_meQTL/Probes/Mononuclear.pdf", plt, width = 4, height = 4, units = "cm")


plt <- enrichments_trans_probes %>%
  filter(Cell == "CD4 naive") %>%
  rename(Term = State) %>%
  plot_odds_ratio2(title = "Long-range meQTL CpGs",
                   bar_size = odds_ratio_errorbar_size,
                   color = col_trans,
                   base_size = font_size_odds_ratio)

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Trans_meQTL/Probes/CD4_naive.pdf", plt, width = 4, height = 4, units = "cm")

plt <- enrichments_trans_probes %>%
  filter(Cell == "CD8 naive") %>%
  rename(Term = State) %>%
  plot_odds_ratio2(title = "Long-range meQTL CpGs",
                   bar_size = odds_ratio_errorbar_size,
                   color = col_trans,
                   base_size = font_size_odds_ratio)

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Trans_meQTL/Probes/CD8_naive.pdf", plt, width = 4, height = 4, units = "cm")


plt <- enrichments_trans_snps %>%
  filter(Cell == "CD4 naive") %>%
  rename(Term = State) %>%
  plot_odds_ratio2(title = "Long-range meQTL SNPs",
                   bar_size = odds_ratio_errorbar_size,
                   color = col_trans,
                   base_size = font_size_odds_ratio)

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Trans_meQTL/SNPs/CD4_naive.pdf", plt, width = 4, height = 4, units = "cm")


plt <- enrichments_trans_snps %>%
  filter(Cell == "Mononuclear") %>%
  rename(Term = State) %>%
  plot_odds_ratio2(title = "Long-range meQTL SNPs",
                   bar_size = odds_ratio_errorbar_size,
                   color = col_trans,
                   base_size = font_size_odds_ratio)

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Trans_meQTL/SNPs/Mononuclear.pdf", plt, width = 4, height = 4, units = "cm")


plt <- enrichments_trans_snps %>%
  filter(Cell == "CD8 naive") %>%
  rename(Term = State) %>%
  plot_odds_ratio2(title = "Long-range meQTL SNPs",
                   bar_size = odds_ratio_errorbar_size,
                   color = col_trans,
                   base_size = font_size_odds_ratio)

ggsave("./Plots/Odds_ratios_bonferroni_884/Roadmap/Trans_meQTL/SNPs/CD8_naive.pdf", plt, width = 4, height = 4, units = "cm")

# Sankey plot -------------------------------------------------------------

# state_flows <- trans_adjust_on_snps_independent %>%
#   left_join(select(meth_roadmap, Probe, Probe_state = State)) %>%
#   rename(SNP_state = State) %>%
#   select(SNP_state, Probe_state) %>%
#   group_by(SNP_state, Probe_state) %>%
#   tally(name = "Freq") %>%
#   ungroup() %>%
#   filter(Freq > 10) %>%
#   mutate(SNP_state = factor(SNP_state, states), Probe_state = factor(Probe_state, states)) %>%
#   mutate(SNPs = fct_recode(SNP_state, !!!roadmap_keys), Probes = fct_recode(Probe_state, !!!roadmap_keys))
#
# plt <- ggplot(data = state_flows,
#        aes(axis1 = SNPs, axis2 = Probes, y = Freq)) +
#   scale_x_discrete(limits = c("SNPs", "CpG"), expand = c(.05, .05)) +
#   geom_alluvium() +
#   ylab("NR meQTLs") +
#   geom_stratum() +
#   geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2) +
#   theme_minimal()
#
# ggsave("./Plots/trans_meqtl_flow.pdf", plt, width = 15, height = 10, units = "cm")
# ggsave("~/Dropbox/Presentations/2020/Lab_meeting_October/Figures/trans_meqtl_flow.pdf", plt, width = 15, height = 10, units = "cm")

# Number of chromosomes as metric -----------------------------------------



# # Presentation ------------------------------------------------------------
# 
# 
# plt_snps <- enrichments_trans_snps %>%
#   filter(Cell == "CD4 naive") %>%
#   rename(Term = State) %>%
#   plot_odds_ratio2(title = "SNPs",
#                    bar_size = odds_ratio_errorbar_size,
#                    color = "black",
#                    base_size = 6)
# 
# plt_cpgs <- enrichments_trans_probes %>%
#   filter(Cell == "CD4 naive") %>%
#   rename(Term = State) %>%
#   plot_odds_ratio2(title = "CpGs",
#                    bar_size = odds_ratio_errorbar_size,
#                    color = "black",
#                    base_size = 6) +
#   theme(axis.line.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks = element_blank())
# 
# plt <- plt_snps + plt_cpgs
# ggsave("./Plots/Odds_ratios/Roadmap/Trans_meQTL/cpg_snp_joint_CD4.pdf", plt, width = 9, height = 6, units = "cm")

