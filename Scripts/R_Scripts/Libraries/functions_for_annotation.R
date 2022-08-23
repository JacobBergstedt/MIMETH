add_15_state_annotation <- function(target, chr_col, position_col, genomic_feature, 
                                    roadmap_desc = read_tsv("./Data/Roadmap_15_states/roadmap_descriptions.tsv")) {
  
  decide_on_duplicates <- function(states) {
    
    if (length(states) >  1) {
      
      has_TSS <- grepl("^TSS", states)
      has_flTSS <- grepl("^Fl. TSS.", states)
      has_enhancer <- grepl("^Enhancers", states)
      has_gen_enhancer <- grepl("^Genic", states)
      has_transcription <- grepl("^Strong transcr", states)
      has_weak_transcription <- grepl("^Weak transcr", states)
      has_3_5_transcription <- grepl("^Transcr.", states)
      has_biv_TSS <- grepl("^Bivalent TSS", states)
      has_biv_fl_TSS <- grepl("^Fl. Bivalent TSS", states)
      has_biv_enhancers <- grepl("^Bivalent Enhancer", states)
      has_ZNF_genes <- grepl("^ZNF genes", states)
      has_repressed_PC <- grepl("^Repressed PC", states)
      has_weak_repressed_PC <- grepl("^Weak Repressed PC", states)
      has_heterochromatin <- grepl("^Heterochromatin", states)
      has_quiescent <- grepl("^Quiescent", states)
      
      if (any(has_TSS)) {
        states[has_TSS] 
      } else if (any(has_flTSS)) {
        states[has_flTSS]
      } else if (any(has_enhancer)) {
        states[has_enhancer]
      } else if (any(has_gen_enhancer)) {
        states[has_gen_enhancer]
      } else if (any(has_transcription)) {
        states[has_transcription]
      } else if (any(has_weak_transcription)) {
        states[has_weak_transcription]
      } else if (any(has_3_5_transcription)) {
        states[has_3_5_transcription]
      } else if (any(has_biv_TSS)) {
        states[has_biv_TSS]
      } else if (any(has_biv_fl_TSS)) {
        states[has_biv_fl_TSS]
      } else if (any(has_biv_enhancers)) {
        states[has_biv_enhancers]
      } else if (any(has_ZNF_genes)) {
        states[has_ZNF_genes]
      } else if (any(has_repressed_PC)) {
        states[has_repressed_PC]
      } else if (any(has_weak_repressed_PC)) {
        states[has_weak_repressed_PC]
      } else if (any(has_heterochromatin)) {
        states[has_heterochromatin]
      } else if (any(has_quiescent)) {
        states[has_quiescent]
      }
      
    } else {
      states
    }
  }
  
  
  overlaps_with_genomic_feature <- function(state_bed, keys, target_range, genomic_feature) {
    
    state_range <- GRanges(seqnames = state_bed$Chromosome, 
                           ranges = IRanges(start = state_bed$Start, end = state_bed$Stop))
    overlaps <- findOverlaps(state_range, target_range)
    add_column(bind_cols(state_bed[from(overlaps), ], target[to(overlaps), ]), Cell = keys$Cell) %>% 
      select(all_of(genomic_feature), Cell, State)
    
  }
  
  keys <- read_tsv("./Data/Roadmap_15_states/code_cell_keys.tsv")
  paths <- glue('./Data/Roadmap_15_states/{keys$Code}_15_coreMarks_mnemonics.bed.gz')
  names(paths) <- keys$Celltype
  state_bed <- vroom(paths,
                     col_names = c("Chromosome", "Start", "Stop", "State"),
                     id = "Cell") %>%
    mutate(Cell = factor(Cell, paths)) %>%
    mutate(Cell = fct_recode(Cell, !!!paths)) %>%
    filter(!Chromosome %in% c("chrX", "chrY"))
  
  target_range <- GRanges(seqnames = target[[chr_col]], 
                          ranges = IRanges(start = target[[position_col]], end = target[[position_col]]))
  
  roadmap_label_keys <- set_names(paste0(roadmap_desc$STATE_NU, "_", roadmap_desc$MNEMONIC),
                                  nm = roadmap_desc$SHORT_DESC)
  
  state_bed %>% 
    group_by(Cell) %>% 
    group_map(overlaps_with_genomic_feature, target_range = target_range, genomic_feature = genomic_feature) %>% 
    bind_rows() %>% 
    mutate(State_labels = as.character(fct_recode(factor(State, roadmap_label_keys), !!!roadmap_label_keys))) %>% 
    select(-State) %>% 
    group_by(.data[[genomic_feature]], Cell) %>% 
    summarize(State_labels = decide_on_duplicates(State_labels)) %>% 
    pivot_wider(names_from = Cell, names_prefix = "Chromatin_states_", values_from = State_labels) %>% 
    ungroup()
  
}

prepare_region_annotation <- function() {
  
  anno <- get_meth_annotation()
  
  genic_regions <- c("Intergenic", "Enhancer", "TFBS","TSS1500", "TSS200",
                     "5'UTR", "1stExon", "Body", "3'UTR")
  
  
  genic_region_annotation <- anno %>% 
    select(meth_anno, Probe, UCSC_RefGene_Group) %>% 
    separate_rows(UCSC_RefGene_Group, sep = ";") %>% 
    group_by(Probe) %>% 
    summarize(Genomic_region = names(which.max(table(UCSC_RefGene_Group))))
  
  geography_annotation <- anno %>% select(Probe, Relation_to_Island) %>%
    mutate(Geography = factor(Relation_to_Island, c("N_Shelf", "N_Shore", "Island", "S_Shore", "S_Shelf", "OpenSea"))) %>% 
    mutate(Geography = fct_recode(Geography, 
                                  `N. shelf` = "N_Shelf",
                                  `S. shelf` = "S_Shelf",
                                  `N. shore` = "N_Shore",
                                  `S. shore` = "S_Shore",
                                  `Open sea` = "OpenSea")) 
  
  list(geography = geography_annotation, genic = genic_region_annotation)
  
}

read_ewas_catalog <- function() {
  col_spec <- cols(Author = col_character(),
                   Consortium = col_character(),
                   PMID = col_integer(),
                   Date = col_date(),
                   Trait = col_character(),
                   EFO = col_character(),
                   Analysis = col_character(),
                   Source = col_character(),
                   Outcome = col_character(),
                   Exposure = col_character(),
                   Covariates = col_character(),
                   Outcome_Units = col_character(),
                   Exposure_Units = col_character(),
                   Methylation_Array = col_character(),
                   Tissue = col_character(),
                   N = col_integer(),
                   N_Cohorts = col_integer(),
                   Categories = col_character(),
                   Age = col_double(),
                   N_Males = col_integer(),
                   N_Females = col_integer(),
                   N_EUR = col_integer(),
                   N_EAS = col_integer(),
                   N_SAS = col_integer(),
                   N_AFR = col_integer(),
                   N_AMR = col_integer(),
                   N_OTH = col_integer(),
                   CpG = col_character(),
                   Location = col_character(),
                   Gene = col_character(),
                   Type = col_character(),
                   Beta = col_double(),
                   SE = col_double(),
                   P = col_double(),
                   Details = col_character())
  anno <- readr::read_tsv("./Data/EWAS_Catalog_23-11-2018.txt",
                  col_types = col_spec)
  names(anno) <- paste0("Catalog_", names(anno))
  anno
}

create_probe_annotation <- function() { 
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(tidyverse)
  
  data(Locations, envir = environment())
  data(Other, envir = environment())
  data(Manifest, envir = environment())
  data(SNPs.Illumina, envir = environment())
  data(Islands.UCSC, envir = environment())

  meth <- get_m_values()
  
  anno <- cbind(Manifest, Locations, Other, SNPs.Illumina, Islands.UCSC) %>%
    as_tibble() %>%
    dplyr::rename(Probe = Name, Probe_chr = chr, Probe_position = pos) %>%
    right_join(tibble(Probe = names(meth)[-1])) %>% 
    dplyr::rename(Probe_gene = UCSC_RefGene_Name)
  
  # Add roadmap enhancer annotation
  anno_roadmap <- readRDS("./Data/RData/Methylation/Annotation/meth_roadmap_annotation.rds") %>% 
    filter(Cell %in% c("Mononuclear_cells", "Neutrophils"), State == "7_Enh") %>% 
    mutate(Enhancer = "Yes") %>% 
    select(Probe, Enhancer) %>% 
    distinct(Probe, .keep_all = TRUE)
  
  anno <- anno %>% 
    left_join(anno_roadmap) %>% 
    mutate(Enhancer = if_else(is.na(Enhancer), "No", Enhancer))
  
  anno_location <- anno %>% 
    dplyr::select(Probe, Probe_chr, Probe_position, Probe_gene)
  
  saveRDS(anno, "./Data/RData/Methylation/Annotation/meth_annotation.rds")
  saveRDS(anno_location, "./Data/RData/Methylation/Annotation/meth_location_annotation.rds")
}

add_probe_annotation <- function(res, probe_annotation, ewas_catalog) {
  probe_annotation <- probe_annotation[c("Probe","Chromosome",  "Position", "SNP_ID", "GencodeBasicV12_NAME")]
  ewas_catalog <- ewas_catalog[c("Catalog_CpG", "Catalog_Trait", "Catalog_Date",  "Catalog_Gene", "Catalog_PMID", "Catalog_P")]
  res <- inner_join(res, ewas_catalog, by = c("Probe" = "Catalog_CpG"))
  res <- inner_join(res, probe_annotation, by = "Probe")
}

add_snp_meth_locations <- function(tib, meth_locations, snp_locations) {
  inner_join(tib, snp_locations, by = c("SNP" = "SNP_ID")) %>% 
    inner_join(meth_locations, by = "Probe") 
}

simplify_probe_genes <- function(x) {
  map_chr(str_split(x, ";"), ~ .[[1]])
}

# Archive -----------------------------------------------------------------
# 
# 
# assemble_roadmap_25_state_annotation <- function() {
#   read_roadmap_tsv <- function(path_file) {
#     head(read_tsv(path_file, col_names = c("Chromosome", "Start", "Stop", "State")), -1)
#   }
#   match_and_clean <- function(roadmap_tib, meth_location) {
#     road_range <- GRanges(seqnames = roadmap_tib$Chromosome,
#                           ranges = IRanges(start = roadmap_tib$Start, end = roadmap_tib$Stop))
#     meth_range <- GRanges(seqnames = paste0("chr", meth_location$Probe_chr),
#                           ranges = IRanges(start = meth_location$Probe_position, stop = meth_location$Probe_position))
#     overlaps <- findOverlaps(road_range, meth_range)
#     bind_cols(roadmap_tib[from(overlaps), ], meth_location[to(overlaps), ])
#   }
#   keys <- read_tsv("./Data/Roadmap_25states/annotation_25_imputed12marks.txt") %>% select(1:3)
#   names(keys) <- c("State", "Mnemonic", "State_description")
#   paths <- c("Bcells",
#              "monocytes",
#              "mononuclear",
#              "neutrophils",
#              "NKcells",
#              "Tcells",
#              "Thelper",
#              "Thelper_memory",
#              "Thelper_naive",
#              "Tkiller_memory",
#              "Tkiller_naive",
#              "Treg")
#   paths_full <- paste0("./Data/Roadmap_25states/", paths, "_25_imputed12marks_stateno.bed.gz")
#   meth_location <- get_meth_location()
#   roadmap_tibs <- map(paths_full, read_roadmap_tsv)
#   names(roadmap_tibs) <- paths
#   map(roadmap_tibs, match_and_clean, meth_location = meth_location) %>% 
#     bind_rows(.id = "Tissue") %>% 
#     left_join(keys)
# }
# add_probe_annotation <- function(res, probe_annotation) {
#   probe_annotation <- probe_annotation[c("Probe","Chromosome",  "Position", "SNP_ID", "GencodeBasicV12_NAME")] %>% setDT()
#   setkey(probe_annotation, Probe)
#   setkey(res, Probe)
#   res[probe_annotation, nomatch = 0]p %>% group_by(Cell) %>% filter(Probe %in% Probe[duplicated(Probe)]) %>% arrange(.by_group = TRUE) %>% print(n = 500)
# }

