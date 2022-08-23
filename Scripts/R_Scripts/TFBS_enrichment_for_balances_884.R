
library(GenomicRanges)
library(rtracklayer)
library(fs)
library(tidyverse)
library(missMethyl)
library(vroom)
library(broom)
library(goseq)
select <- dplyr::select
library(scales)
source("./Scripts/R_scripts/Libraries/functions_for_enrichments.R")
source("./Scripts/R_scripts/Libraries/functions_for_plotting.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")

TF_meth_dep <- scan("./Data/meth_binding_TF.txt", what = character())
TF_CD8 <- c("TBX21", "ID3", "EOMES", "BCL6", "FOXO1", "STAT3", "ZEB1", "BACH1")

# Setup data structures ---------------------------------------------------


# path_to_chain = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
# lift_over_chain = import.chain(path_to_chain)
# 
# meth_location <- get_meth_location(fix = FALSE)
# meth_anno_roadmap <- get_meth_roadmap_annotation()
# 
# TFs <- list.files("./Data/TFBS_2019/pwm_tfbs_per_tf")
# TFBS <- map(TF. ~ add_TFBS(., meth_anno = meth_location, lift_over_chain = lift_over_chain)) %>%
#   reduce(inner_join)
# 
# TFBS <- bind_cols(TFBS["Probe"], discard(TFBS[-1], ~ min(table(.)) < 100)) %>% 
#   rename_with(.fn = ~ paste0("TF_", .), .cols = -Probe)
# 

# Get probe lists ---------------------------------------------------------

# dice <- get_dice() %>% 
#   pivot_longer(cols = -Gene, names_to = "Cell", values_to = "Expression")

meth_location <- get_anno_meth_location(fix = FALSE)
meth_anno_roadmap <- get_anno_meth_roadmap_all_cells() %>% 
  select(Probe, State = Chromatin_states_CD4_naive)

meth_anno_loc_roadmap <- meth_location %>% 
  inner_join(meth_anno_roadmap)

meth_anno_TFBS <- readRDS("./Data/RData/Methylation/Annotation/meth_TFBS_annotation.rds") %>% 
  left_join(meth_anno_roadmap)



# Balance probes ----------------------------------------------------------

keep <- c("Mye_Lym", "NK_Lymph", "B_Lymph", "CD4_Tcells", "CD8Naive_CD8", "CD4Naive_CD4")

res_balances <- readRDS("./Data/RData/Results/EWAS/M_values/Environment/EWAS_balance_effects_compositional_analysis.rds") %>% 
  filter(Balance %in% keep) %>% 
  group_by(Balance) %>% 
  mutate(P_bonferroni = p.adjust(P_value)) %>% 
  ungroup()

get_probes <- function(x, label) {
  
  pos <- x %>% filter(P_bonferroni < 0.05, Estimate > 0) %>% pull(Probe)
  neg <- x %>% filter(P_bonferroni < 0.05, Estimate < 0) %>% pull(Probe)
  out <- list(pos, neg)
  names(out) <- paste0(label, "_", c("pos", "neg"))
  out
  
}

f_split <- factor(res_balances$Balance)
res_list <- res_balances %>% 
  split(f_split) %>% 
  map2(levels(f_split), get_probes) %>% 
  flatten()

res_enrichments <- mclapply(res_list, get_TFBS_enrichment, anno = meth_anno_TFBS, mc.cores = 6) %>% 
  bind_rows(.id = "Term")

