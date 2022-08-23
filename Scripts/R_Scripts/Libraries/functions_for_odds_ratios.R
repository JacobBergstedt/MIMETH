# Declare functions --------------------------------------

fit_logistic <- function(df, keys) {
  
  if (df$Freq[1] != 0) {
    m <- glm(Proportion_affected ~ Treatment, data = df, family = "binomial", weights = Total)
    summary_m <- summary(m)
    ci <- exp(confint(m, level = (1 - 0.01))["TreatmentTRUE", ])
    low <- ci[1]
    high = ci[2]
    coefs <- summary_m$coefficients["TreatmentTRUE", c("Estimate", "Pr(>|z|)")]
  } else {
    coefs <- c(-2, 1e-20)
    low = -2
    high = -1.99
  }
  tibble(Estimate = exp(coefs[1]), P_value = coefs[2], Low = low, High = high, 
         Cases_affected = df$Freq[df$Treatment == TRUE], 
         Cases_total = df$Total[df$Treatment == TRUE],
         Non_cases_affected = df$Freq[df$Treatment == FALSE],
         Non_cases_total = df$Total[df$Treatment == FALSE])
}


compare_cpg_geography_proportions <- function(probes, anno) {
  
  in_group <- anno %>% filter(Probe %in% probes)
  out_group <- anno %>% filter(!Probe %in% probes)
  
  tallies_in <- table(in_group$Geography)
  tallies_out <- table(out_group$Geography)
  tallies <- bind_rows(tibble(Term = names(tallies_in), 
                              Freq = tallies_in, 
                              Total = nrow(in_group), 
                              Proportion_affected = Freq / Total, 
                              Treatment = TRUE),
                       tibble(Term = names(tallies_out),
                              Freq = tallies_out,
                              Total = nrow(out_group),
                              Proportion_affected = Freq / Total,
                              Treatment = FALSE))
  
  tallies %>%
    group_by(Geography) %>% 
    group_modify(fit_logistic) %>% 
    ungroup()
}

compare_roadmap_region_proportions <- function(probes, anno) {
  
  anno$State <- factor(anno$State, 
                       unique(anno$State)[order(as.integer(str_match(unique(anno$State), "[0-9]*")))])
  
  in_group <- anno %>% filter(Probe %in% probes)
  out_group <- anno %>% filter(!Probe %in% probes)
  
  build_proportion_tib <- function(tib, keys) {
    
  
    in_table <- tib$In_table[[1]]
    out_table <- tib$Out_table[[1]]
    
    in_tib <- tibble(Region = names(in_table), 
                     Freq = in_table, 
                     Total = tib$In_n, 
                     Proportion_affected = Freq / Total, 
                     Treatment = TRUE)
    
    out_tib <- tibble(Region = names(out_table),
                      Freq = out_table,
                      Total = tib$Out_n,
                      Proportion_affected = Freq / Total,
                      Treatment = FALSE)
    
    bind_rows(in_tib, out_tib)
  }
  
  tallies_in <- in_group %>% 
    group_by(Cell) %>% 
    summarize(In_n = n(), In_table = list(table(State)))
  
  tallies_out <- out_group %>% 
    group_by(Cell) %>% 
    summarize(Out_n = n(), Out_table = list(table(State)))
  
  tallies <- inner_join(tallies_in, tallies_out) %>% 
    group_by(Cell) %>% 
    group_modify(build_proportion_tib) %>% 
    ungroup()
  
  # with_enough_cases <- tallies %>% 
  #   group_by(Cell, Region) %>% 
  #   summarize(min_prop = min(Proportion_affected)) %>% 
  #   filter(min_prop > 1e-2) %>% 
  #   select(-min_prop)
  
  # tallies %>% 
  #   right_join(with_enough_cases) %>% 
  #   group_by(Cell, Region) %>% 
  #   group_modify(fit_logistic) %>% 
  #   ungroup()
  
  tallies %>% 
    group_by(Cell, Region) %>% 
    group_modify(fit_logistic) %>% 
    ungroup()
}

compare_genic_region_proportions <- function(probes, anno) {
  
  in_group <- anno %>% filter(Probe %in% probes)
  out_group <- anno %>% filter(!Probe %in% probes)
  
  tallies_in <- table(in_group$Region)
  tallies_out <- table(out_group$Region)
  tallies <- bind_rows(tibble(Region = names(tallies_in), 
                              Freq = tallies_in, 
                              Total = nrow(in_group), 
                              Proportion_affected = Freq / Total, 
                              Treatment = TRUE),
                       tibble(Region = names(tallies_out),
                              Freq = tallies_out,
                              Total = nrow(out_group),
                              Proportion_affected = Freq / Total,
                              Treatment = FALSE))
  
  tallies %>%
    group_by(Region) %>% 
    group_modify(fit_logistic) %>% 
    ungroup()
  
}


compare_consensus_roadmap <- function(probes, 
                                      anno = readRDS("./Data/RData/Methylation/Annotation/meth_roadmap_consensus_annotation.rds"), 
                                      desc = read_tsv("./Data/Roadmap_15_states/roadmap_descriptions.tsv")) {
  
  
  anno$Case <- factor(anno$Probe %in% probes)
  anno$State <- fct_relevel(factor(anno$State), "15_Quies")  
  
  states_without_obs <- as.data.frame(xtabs(~ State + Case, anno)) %>% 
    filter(Freq == 0) %>% 
    pull(State) %>% 
    as.character()
  
  if (!is_empty(states_without_obs)) anno <- anno %>% filter(!State %in% states_without_obs) %>% mutate(State = droplevels(State))
  
  m <- glm(Case ~ State, family = "binomial", data = anno, control = glm.control(epsilon = 1e-10, maxit = 1e3))
  
  roadmap_label_keys <- set_names(glue('State{desc$STATE_NU}_{desc$MNEMONIC}'), 
                                  desc$SHORT_DESC)
  
  
  CIs <- exp(confint(m, level = 0.99))[-1, ]
  out <- tidy(m) %>%
    dplyr::rename(Term = term, Estimate = estimate, P_value = p.value, Standard_error = std.error) %>%     
    dplyr::filter(Term != "(Intercept)") %>%
    add_row(Term = "State15_Quies", Estimate = log(1), Standard_error = 0, statistic = 0, P_value = 0) %>% 
    mutate(Estimate = exp(Estimate),
           Term = factor(Term, Term[order(as.integer(str_match(unique(Term), "[0-9]+")))]),
           Low = c(CIs[, 1], 1),
           High = c(CIs[, 2], 1))
  
  if (!is_empty(states_without_obs)) out <- out %>% add_row(Term = paste0("State", states_without_obs), Estimate = NA, Standard_error = NA, statistic = NA, P_value = NA, Low = NA, High = NA)
  
  out %>% mutate(Term = fct_recode(Term, !!!roadmap_label_keys))
  
}


# 
# compare_roadmap <- function(probes, anno) {
#   
#   db <- right_join(tibble(Probe = probes, Case = "Yes"), anno_roadmap) %>% 
#     mutate(Case = if_else(is.na(Case), "No", "Yes"), 
#            Case = factor(Case, c("No", "Yes")), 
#            State = relevel(factor(State), "15_Quies"))
#   
#   db %>% 
#     group_by(Cell) %>% 
#     group_modify(~ tidy(glm(Case ~ State, data = ., family = "binomial"))) %>% 
#     ungroup() %>% 
#     mutate(risk = exp(estimate))
#   
# }

compare_geography <- function(probes, anno) {

  db <- left_join(anno, tibble(Probe = probes, Case = TRUE)) %>% 
    mutate(Case = ifelse(is.na(Case), FALSE, TRUE))
  
  m <- glm(Case ~ Geography, family = "binomial", data = db)
  geography_keys <- set_names(paste0("Geography", c("N_Shelf", "N_Shore", "Island", "S_Shore", "S_Shelf")),
                              c("N. shelf", "N. shore", "Island", "S. shore", "S. shelf"))
  
  CIs <- exp(confint(m, level = 0.99))[-1, ]
  tidy(m) %>%
    dplyr::rename(Term = term, Estimate = estimate, P_value = p.value, Standard_error = std.error) %>%     
    dplyr::filter(Term != "(Intercept)") %>%
    mutate(Estimate = exp(Estimate),
           Term = factor(Term, Term),
           Low = CIs[, 1],
           High = CIs[, 2]) %>%
    mutate(Term = fct_recode(Term, !!!geography_keys))
  
  
}



