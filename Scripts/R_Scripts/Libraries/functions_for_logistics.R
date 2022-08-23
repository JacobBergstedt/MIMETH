

get_dice <- function() {
  readRDS("./Data/Dice_database/dice.rds")
}

get_snp_matrix_path <- function(chr = NULL) {
  if (is.null(chr)) "./Data/RData/Genotypes/LabExMI_imputation_1000x5699237_snp_matrix.rds"
  else paste0("./Data/Chunk_data/Genotype/Chromosome_", chr, "/LabExMI_imputation_958x5699237_chr", chr, "_snp_matrix.rds")
}

get_snp_location_path <- function(chr = NULL) {
  if (!is.null(chr)) {
    paste0("./Data/Chunk_data/Genotype/Chromosome_", chr, "/LabExMI_imputation_1000x5699237_chr", chr, "_map.rds")
  } else {
    "./Data/RData/Genotypes/LabExMI_imputation_1000x5699237_annotated_map.rds"      
  }
}

get_cell_list <- function() {
  scan("./Data/prop_controls_16.txt", what = character())
}

get_cell_list_lineage <- function() {
  
  c("X_NK_of_total.lineage_panel",
    "X_mono_of_total.lineage_panel",
    "X_CD4_of_total.lineage_panel",
    "X_CD8_of_total.lineage_panel",
    "X_CD19pos_of_total.lineage_panel",
    "X_CD4negCD8neg_of_total.lineage_panel",
    "X_CD8bposCD4pos_of_total.lineage_panel",
    "X_neutrophils_of_total.lineage_panel")
}


get_CD8_cell_list <- function() {
  
  c("X_CD8b_naive_of_CD8b.panel1",
    "X_CD8b_CM_of_CD8b.panel1",
    "X_CD8b_EM_of_CD8b.panel1",
    "X_CD8b_EMRA_of_CD8b.panel1")
}


get_CD8_cells <- function() {
  
  readRDS("./Data/RData/Cells/proportions_CD8_subsets.rds")
  
}

get_panel5_cells <- function() {
  
  readRDS("./Data/RData/Cells/proportions_panel5_subsets.rds")
  
}

get_cell_list_roadmap <- function() {
  read_tsv("./Data/Roadmap_15_states/code_cell_keys.tsv")$Celltype
}

get_enrichment_snp_table <- function() {
  read_tsv("./Scripts/Bash_scripts/GWAS_enrichment/Databases/LabExMI_SNPs.txt")
}

get_snp_mat_per_probe <- function() {
  readRDS("./Data/RData/Genotypes/snp_mat_per_probe.rds")
}

get_snp_mat_per_probe_sample <- function() {
  readRDS("./Data/RData/Genotypes/snp_mat_per_probe_1e3.rds")
}

get_m_values_sample <- function() {
  readRDS("./Data/RData/Methylation/meth_1e3.rds")
}

get_m_values_path <- function(chr = NULL) {
  if (!is.null(chr)) paste0("./Data/Chunk_data/Methylation/M_values/Per_chromosome/meth_chromosome_", chr, ".rds") 
  else "./Data/RData/Methylation/MIMETH.minfi.final.MMatrix.autosomes.no_outliers.rds"
}

get_meth_id_path <- function() {
  "./Data/RData/Methylation/subject_ids.rds"
}

get_anno_meth_annotation_path <- function() {
  "./Data/RData/Methylation/Annotation/meth_annotation.rds"
}

get_anno_meth_location_path <- function() {
  "./Data/RData/Methylation/Annotation/meth_location_annotation.rds"
}

get_snp_matrix <- function(chr = NULL) {
  readRDS(get_snp_matrix_path(chr))
}

get_anno_snp_location <- function(chr = NULL) {
  readRDS(get_snp_location_path(chr))
}

get_m_values <- function(chr = NULL) {
  readRDS(get_m_values_path(chr))
}

get_beta_values <- function(chr = NULL) {
  readRDS(get_beta_values_path(chr))
}

get_beta_values_path <- function(chr = NULL) {
  if (!is.null(chr)) paste0("./Data/Chunk_data/Methylation/Beta/Per_chromosome/meth_chromosome_", chr, ".rds") 
  else "./Data/RData/Methylation/MIMETH.minfi.final.betaMatrix.autosomes.no_outliers.rds"
}

get_meth_ids <- function() {
  readRDS(get_meth_id_path())
}

get_anno_meth <- function() {
  readRDS(get_anno_meth_annotation_path())
}

get_anno_meth_location <- function(fix = TRUE) {
  meth_loc <- readRDS(get_anno_meth_location_path())
  if (fix) {
    meth_loc <- meth_loc %>% 
      filter(Probe_chr != "chrX") %>% 
      mutate(Probe_chr = as.integer(gsub("chr", "", Probe_chr)))
  }
  meth_loc
}

translate_roadmap_labels <- function(states) {
  roadmap_desc <- read_tsv("./Data/Roadmap_15_states/roadmap_descriptions.tsv")
  roadmap_label_keys <- set_names(paste0(roadmap_desc$STATE_NU, "_", roadmap_desc$MNEMONIC),
                                  nm = roadmap_desc$SHORT_DESC)
  factor(states, get_anno_roadmap_states()) %>% 
    fct_recode(!!!roadmap_label_keys)
}

# get_anno_meth_roadmap <- function() {
#   roadmap_desc <- read_tsv("./Data/Roadmap_15_states/roadmap_descriptions.tsv")
#   roadmap_label_keys <- set_names(paste0(roadmap_desc$STATE_NU, "_", roadmap_desc$MNEMONIC),
#                                   nm = roadmap_desc$SHORT_DESC)
#   readRDS("./Data/RData/Methylation/Annotation/meth_anno_roadmap_consensus4.rds") %>%
#     mutate(State_labels = fct_recode(factor(State, roadmap_label_keys), !!!roadmap_label_keys))
# }

get_anno_meth_roadmap <- function() {
  
  readRDS("./Data/RData/Methylation/Annotation/meth_anno_roadmap.rds") %>% 
    dplyr::select(Probe, State_labels = Chromatin_states_Mononuclear_cells)
  
}

get_anno_meth_roadmap_all_cells <- function() {
  
  readRDS("./Data/RData/Methylation/Annotation/meth_anno_roadmap.rds")
  
}

get_anno_meth_geography <- function() {
  
  meth_anno <- get_anno_meth()
  
  meth_anno_geography <- meth_anno %>%
    select(Probe, Relation_to_Island) %>%
    mutate(Geography = factor(Relation_to_Island, c("OpenSea", "N_Shelf", "N_Shore", "Island", "S_Shore", "S_Shelf"))) %>%
    select(-Relation_to_Island)
  
  geography_keys <- set_names(c("OpenSea", "N_Shelf", "N_Shore", "Island", "S_Shore", "S_Shelf"),
                              c("Open sea", "N. shelf", "N. shore", "Island", "S. shore", "S. shelf"))  
  
  meth_anno_geography$Geography_labels <- fct_recode(factor(meth_anno_geography$Geography), !!!geography_keys)
  
  meth_anno_geography
  
}

get_anno_roadmap_states <- function() {
  roadmap_desc = read_tsv("./Data/Roadmap_15_states/roadmap_descriptions.tsv")
  states <- paste0(roadmap_desc$STATE_NU, "_", roadmap_desc$MNEMONIC)
  states[order(as.numeric(str_match(unique(states), "[0-9]*")))]
}

get_anno_roadmap_translation <- function() {
  roadmap_desc = read_tsv("./Data/Roadmap_15_states/roadmap_descriptions.tsv")
  states <- paste0(roadmap_desc$STATE_NU, "_", roadmap_desc$MNEMONIC)
  set_names(states, roadmap_desc$SHORT_DESC)
}


get_covs <- function() {
  readRDS("./Data/RData/Environment/covariates.rds")
}

get_covs_884 <- function() {
  readRDS("./Data/RData/Environment/covariates_884.rds")
}

get_trans <- function() {
  readRDS("./Data/RData/Results/MeQTL/trans.rds")
}

get_cis <- function() {
  readRDS("./Data/RData/Results/MeQTL/cis.rds")
}

get_cell_specific_cis <- function() {
  readRDS("./Data/RData/Results/MeQTL/cell_specific_cis.rds")
}

get_cell_specific_trans <- function() {
  readRDS("./Data/RData/Results/MeQTL/cell_specific_trans.rds")
}

get_immune_cell_results <- function() {
  readRDS("./Data/RData/Results/EWAS/M_values/Constant_variance/Immune_cells/correct_for_standard_controls.rds")
}

get_proportion_list <- function() {
  scan("../../../GWAS/LMM/FACS_GWAS/FACS_DataAnalysis/Data/proportion_list.txt", character())
}

get_15_props <- function() {
  
  c("X_VIABLE_BASOPHILS_OF_TOTAL.panel7",
    "X_VIABLE_EOSINOPHILS_OF_TOTAL.panel7",
    "X_mono_of_total.panel5",
    "X_CD4_naive_of_total.panel1",
    "X_CD4_CM_of_total.panel1",
    "X_CD4_EM_of_total.panel1",
    "X_CD4_EMRA_of_total.panel1",
    "X_CD8_naive_of_total.panel1",
    "X_CD8_CM_of_total.panel1",
    "X_CD8_EM_of_total.panel1",
    "X_CD8_EMRA_of_total.panel1",
    "X_CD8bnegCD4neg_of_total.panel1",
    "X_Bcells_ofTotal.panel6",
    "X_NK_of_total.panel4",
    "X_dendritic_cells.panel8")
  
}

get_snps_chunk <- function(chr, chunk) {
  readRDS(paste0("./Data/Chunk_data/Genotype/Chromosome_", chr, "/Chunks/snps_chunk_", chunk, ".rds"))
}

get_meth_chunk <- function(chr, chunk) {
  readRDS(paste0("./Data/Chunk_data/Methylation/Beta/Per_chromosome/Chromosome_", chr, "/Chunks/Chunk_", chunk, ".rds"))  
}

get_probe_snps_chunk <- function(chr, chunk) {
  readRDS(paste0("./Data/Chunk_data/Genotype/Chromosome_", chr, "/Chunks/probe_snps_chunk_", chunk, ".rds"))
}

get_covariate_list <- function() {
  scan("./Data/covariate_list.txt", what = "character")
}