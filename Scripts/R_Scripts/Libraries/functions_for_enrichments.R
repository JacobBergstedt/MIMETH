
add_cell_counts <- function(x) {
  
  x <- as.data.frame(x)
  vec2tib(set_names(x$Freq, c("Not_Assoc_Not_TFBS", "Assoc_Not_TFBS", "Not_Assoc_TFBS", "Assoc_TFBS")))
  
}


get_TFBS_xtabs <- function(probe_list, db) {
  
  # active_states <- c("1_TssA", "2_TssAFlnk", "3_TxFlnk", 
  #                    "4_Tx", "5_TxWk", "7_Enh", "6_EnhG")
  
  db <- db %>% 
    mutate(Target = as.factor(Probe %in% probe_list))
  
  map(select(db, -Probe, -Target, -State), function(x) xtabs(~ db$Target + x))
  
}


get_TFBS_enrichment <- function(probes, anno) {
  
  tabs <- get_TFBS_xtabs(probes, anno) %>% 
    keep(screen_TFBS)
  
  cell_counts <- tabs %>% map_dfr(add_cell_counts, .id = "TF")
  
  tabs %>% 
    map(fisher.test) %>%
    map_dfr(tidy, .id = "TF") %>% 
    mutate(P_FDR = p.adjust(p.value, "BH")) %>% 
    inner_join(cell_counts)
  
  
}

screen_TFBS <- function(x) {
  
  target_hits <- as.data.frame(x) %>% filter(x == "Yes") %>% pull(Freq)
  sum(target_hits) > 200
  
}

build_EWAS_enrichment_tables_bonferroni <- function(tib, thresh = 0.05) {
  
  states <- names(get_anno_roadmap_translation())
  cell_lines <- get_cell_list_roadmap()
  cell_lines_label <- gsub("_", " ", gsub("_cells", "", cell_lines))
  cell_lines <- paste0("Chromatin_states_", cell_lines)
  
  
  tib_pos <- tib %>% 
    filter(!(P_bonferroni < thresh & Estimate < 0)) %>% 
    mutate(Case = (P_bonferroni < thresh & Estimate > 0))
  
  tib_neg <- tib %>% 
    filter(!(P_bonferroni < thresh & Estimate > 0)) %>% 
    mutate(Case = (P_bonferroni < thresh & Estimate < 0))
  
  
  res_pos <- vector(mode = "list", length = length(cell_lines))
  res_neg <- vector(mode = "list", length = length(cell_lines))
  names(res_pos) <- cell_lines_label
  names(res_neg) <- cell_lines_label
  
  pos_is_case <- tib_pos$Case
  neg_is_case <- tib_neg$Case
  
  for (i in 1:length(cell_lines)) {
    
    cell_states <- tib_pos[[cell_lines[[i]]]]
    res_pos[[i]] <- map(set_names(states, states), ~ cell_states == .) %>% 
      map(function(is_state) xtabs(~ is_state + pos_is_case))
    
    cell_states <- tib_neg[[cell_lines[[i]]]]
    res_neg[[i]] <- map(set_names(states, states), ~ cell_states == .) %>% 
      map(function(is_state) xtabs(~ is_state + neg_is_case))
    
  }
  
  list(neg = res_neg, pos = res_pos)
  
}

screen_EWAS <- function(x) {
  
  x <- as.data.frame(x)
  target_hits <- x %>% filter(as.logical(x[[2]])) %>% pull(Freq)
  all(target_hits > 0) & any(target_hits > 1)
  
}



do_fisher_tests_for_EWAS <- function(table_list) {
  
  unlist2(unlist2(table_list)) %>%
    keep(screen_EWAS) %>% 
    map(fisher.test) %>% 
    map_dfr(tidy, .id = "ID") %>%
    separate(col = ID, into = c("Direction", "Cell", "State"), sep = "__") %>% 
    rename(Estimate = estimate, Low = conf.low, High = conf.high, P_value = p.value)
  
}


run_EWAS_enrichment_pipeline_bonferroni <- function(tib) {
  
  build_EWAS_enrichment_tables_bonferroni(tib) %>% 
    do_fisher_tests_for_EWAS()
}


go_enrichment_TFBS <- function(hits, TF_list) {
  
  genes <- set_names(rep(0, length(TF_list)), TF_list)
  genes[names(genes) %in% hits] <- 1
  
  pwf <- nullp(genes,"hg19","geneSymbol", plot.fit = FALSE)
  go_enrichment_TF <- as_tibble(goseq(pwf, "hg19", "geneSymbol"))
  go_enrichment_TF %>% 
    filter(!is.na(term)) %>% 
    mutate(P_FDR = p.adjust(over_represented_pvalue, "BH")) %>% 
    arrange(over_represented_pvalue)
  
}



gometh_enrichment_EWAS <- function(tib) {
  
  probes_control <- tib %>%
    filter(P_bonferroni > 0.05) %>%
    pull(Probe)
  
  target_pos_10K <- tib %>%
    filter(Estimate > 0) %>%
    arrange(P_bonferroni) %>%
    dplyr::slice(1:10000) %>%
    pull(Probe)
  
  target_pos_pval <- tib %>%
    filter(Estimate > 0, P_bonferroni < 0.05) %>%
    pull(Probe)
  
  target_neg_10K <- tib %>%
    filter(Estimate < 0) %>%
    arrange(P_bonferroni) %>%
    dplyr::slice(1:10000) %>%
    pull(Probe)
  
  target_neg_pval <- tib %>%
    filter(Estimate < 0, P_bonferroni < 0.05) %>%
    pull(Probe)
  
  list(Pos_10K = target_pos_10K,
       Pos_Pval = target_pos_pval,
       Neg_10K = target_neg_10K,
       Neg_Pval = target_neg_pval) %>%
    map_dfr(gometh, all.cpg = probes_control, array.type = "EPIC", collection = "GO", .id = "Terms")
  
  
}
