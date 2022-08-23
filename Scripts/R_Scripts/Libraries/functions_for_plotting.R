

# Global variables --------------------------------------------------------

# Pick colors -------------------------------------------------------------

oranges <- brewer_pal(palette = "Oranges")(9)
reds <- brewer_pal(palette = "Reds")(9)
blues <- brewer_pal(palette = "Blues")(9)
greens <- brewer_pal(palette = "Greens")(9)
col_age_cell <- blues[4]
col_age <- blues[7]
col_sex <- blues[9]


col_CMV <- reds[8]
col_smoking <- reds[6]
col_CRP <- oranges[6]
col_heartrate <- oranges[4]
col_temperature <- oranges[3]
col_cell <- oranges[5]

col_cis <- greens[6]
col_trans <- greens[8]

col_gen <- greens[7]
col_intrinsic <- blues[6]
col_exposures <- reds[7]
State_labels <- oranges[5]


# Point sizes -------------------------------------------------------------

line_size_blacks <- 0.5
line_size_scatter <- 0.5
line_size_violin <- 0.3
point_size_comparison <- 0.5
point_size_scatter <-  0.5
point_size_volcano <- 0.5
odds_ratio_errorbar_size <- 0.3
odds_ratio_point_size <- 0.3
point_size_outlier <- 0.05

# Font sizes --------------------------------------------------------------

label_size_volcano <- 7
font_size_tag <- 7
font_size_volcano <- 7
font_size_odds_ratio <- 7
font_size_scatter <- 7
font_size_correction_comparison <- 7
font_size_violin <- 7
font_size_table <- 7
font_size_label = font_size_table * (1 / ggplot2:::.pt)


# Functions ---------------------------------------------------------------


plot_beta_effects <- function(effect_sizes, p_values, is_cis = FALSE, is_age = FALSE) {
  
  plt_tib <- tibble(Effect_sizes = effect_sizes, 
                    P_values = p_values)
  
  plt <- ggplot(plt_tib, aes(x = Effect_sizes, y = -log10(P_values))) +
    geom_point(size = 2, alpha = 0.3) +
    scale_color_brewer(type = "qual", palette = 3) +
    xlab("Effect size") +
    ylab(TeX("-log_{10} (P_{adjusted})")) +
    geom_vline(xintercept = min(plt_tib$Effect_sizes), linetype="dotted") +
    geom_vline(xintercept = max(plt_tib$Effect_sizes), linetype="dotted") +
    geom_vline(xintercept = min(plt_tib$Effect_sizes) / 2, linetype="dotted") +
    geom_vline(xintercept = max(plt_tib$Effect_sizes) / 2, linetype="dotted") +
    scale_x_continuous(breaks = seq(-11, 11, 2) / 10, limits = c(-11, 11) / 10) +
    scale_y_continuous(limits = c(0, 310)) +
    theme_minimal(base_size = 16) +
    theme(axis.text = element_text(color = "black"),
          panel.grid = element_blank(),
          legend.position = "none")
  
  plt
}

plot_m_effects_without_counts <- function(tib, estimate_col, P_val_col, color, title1, base_size = 12, point_size, label_size = base_size) {
  
  ggplot(tib, aes(x = .data[[estimate_col]], y = -log10(.data[[P_val_col]]), color = color)) +
    geom_point(size = point_size, alpha = 0.8) +
    ggtitle(title1) + 
    xlab("Effect size") +
    ylab("-log10(Padj)") +
    geom_vline(xintercept = 0, color = "black", size = line_size_blacks) +
    scale_color_manual(values = color) +
    scale_x_continuous(breaks = seq(-10, 10, 4), limits = c(-10, 10)) +
    scale_y_continuous(limits = c(0, 320)) +
    theme_classic(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color = "black"),
          panel.grid = element_blank(),
          legend.position = "none",
          panel.border = element_blank())

  
}

plot_m_effects <- function(tib, estimate_col, P_val_col, color, title1, title2, base_size = 12, point_size, label_size = base_size, y_limit = c(0, 200000)) {

  plt <- ggplot(tib, aes(x = .data[[estimate_col]], y = -log10(.data[[P_val_col]]), color = color)) +
    geom_point(size = point_size, alpha = 0.8) +
    ggtitle(title1) +
    xlab("Effect size") +
    ylab("-log10(Padj)") +
    geom_vline(xintercept = 0, color = "black", size = line_size_blacks) +
    scale_color_manual(values = color) +
    scale_x_continuous(breaks = seq(-10, 10, 4), limits = c(-10, 10)) +
    scale_y_continuous(limits = c(0, 320)) +
    theme_classic(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color = "black"),
          panel.grid = element_blank(),
          legend.position = "none",
          panel.border = element_blank())

  less_than_zero <- sum(tib[[estimate_col]] < 0)
  larger_than_zero <- sum(tib[[estimate_col]] > 0)

  tib_counts <- tibble(N = c(less_than_zero, larger_than_zero), Type = c("Neg.", "Pos."))
  plt1 <- ggplot(tib_counts, aes(x = Type, y = N)) +
    geom_col(fill = color) +
    scale_y_continuous(labels = comma, limits = y_limit) +
    ggtitle(title2) +
    xlab("Direction") +
    ylab(NULL) +
    theme_classic(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color = "black"),
          panel.grid = element_blank(),
          legend.position = "none",
          panel.border = element_blank())

  plt + plt1 + plot_layout(widths = c(4, 1))

}

plot_volcano <- function(tib, estimate_col, P_val_col, color, title1, title2, base_size = 12, point_size, label_size = base_size) {
  
  plt <- ggplot(tib, aes(x = .data[[estimate_col]], y = -log10(.data[[P_val_col]]), color = color)) +
    geom_point(size = point_size, alpha = 0.8) +
    ggtitle(title1) + 
    xlab("Effect") +
    ylab("-log10(Padj)") +
    geom_vline(xintercept = 0, color = "black", size = line_size_blacks) +
    scale_color_manual(values = color) +
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color = "black"),
          panel.grid = element_blank(),
          legend.position = "none")
  
  less_than_zero <- sum(tib[[estimate_col]] < 0)
  larger_than_zero <- sum(tib[[estimate_col]] > 0)
  
  tib_counts <- tibble(N = c(less_than_zero, larger_than_zero), Type = c("Neg.", "Pos."))
  plt1 <- ggplot(tib_counts, aes(x = Type, y = N)) +
    geom_col(fill = color) +
    scale_y_continuous(labels = comma) + 
    ggtitle(title2) +
    xlab("Direction") +
    ylab(NULL) + 
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color = "black"),
          panel.grid = element_blank(),
          legend.position = "none")
  
  plt + plt1 + plot_layout(widths = c(4, 1))
  
}


add_quiescent <- function(tib_odds_ratio) {
  add_row(tib_odds_ratio, Term = "Quiescent", Estimate = 0, Standard_error = 0, P_value = 1, Low = 0, High = 0)
}

plot_odds_ratio <- function(tib_odds_ratio, title, bar_size, color, limits = NULL, breaks = waiver(), base_size) {

  tib_odds_ratio %>%
    ggplot(aes(y = Estimate, ymin = Low, ymax = High, x = fct_rev(Term))) +
    geom_errorbar(size = bar_size, width = 0.3, col = color) +
    geom_hline(yintercept = 1, size = line_size_blacks, color = "black") +
    scale_y_continuous(limits = limits, breaks = breaks) +
    xlab(NULL) +
    ylab(NULL) +
    ggtitle(title) +
    geom_point(size = odds_ratio_point_size, col = color) +
    coord_flip() +
    theme_bw(base_size = base_size) +
    theme(axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5),
          panel.grid.minor.x = element_blank())
}

plot_odds_ratio_split <- function(tib_odds_ratio, title, point_size, bar_size, colors, limits = NULL, breaks = waiver(), base_size, legend_position) {
  
  tib_odds_ratio %>% 
    ggplot(aes(y = Estimate, 
               ymin = Low, 
               ymax = High, 
               x = fct_rev(Term), 
               color = Direction)) +
    geom_errorbar(size = bar_size, width = bar_size) +
    geom_hline(yintercept = 1, size = line_size_blacks, color = "black") +
    # scale_y_log10(limits = c(1e-3, 1e2), breaks = c(1e-3, 1e-2, 1e-1, 1, 1e1, 1e2)) +
    ylab(NULL) +
    xlab(NULL) +
    geom_point(size = point_size) +
    ggtitle(title) +
    scale_color_manual(values = colors) +
    guides(color = guide_legend(title = "Dir.")) + 
    coord_flip() +
    theme_bw(base_size = base_size) +
    theme(axis.text = element_text(color = "black"),
          legend.position = legend_position,
          plot.title = element_text(hjust = 0.5),
          panel.grid.minor.x = element_blank(),
          legend.key.size = unit(0.1, "lines"))
}

plot_odds_ratio_linear <- function(tib_odds_ratio, title, point_size, bar_size, colors, limits = NULL, breaks = waiver(), base_size, legend_position) {
  
  tib_odds_ratio %>% 
    ggplot(aes(y = Estimate, 
               ymin = Low, 
               ymax = High, 
               x = Term, 
               color = Direction)) +
    geom_errorbar(size = bar_size, width = 0.2) +
    geom_hline(yintercept = 1, size = line_size_blacks) +
    ylab(NULL) +
    xlab(NULL) +
    geom_point(size = point_size) +
    ggtitle(title) +
    scale_color_manual(values = colors) +
    coord_flip() +
    theme_bw(base_size = base_size) +
    theme(axis.text = element_text(color = "black"),
          legend.position = legend_position,
          plot.title = element_text(hjust = 0.5),
          panel.grid.minor.x = element_blank(),
          legend.key.size = unit(0.1, "lines"))
}



plot_odds_ratio2 <- function(tib_odds_ratio, title, bar_size, color, limits = NULL, breaks = waiver(), base_size) {
  
  tib_odds_ratio %>%
    ggplot(aes(y = Estimate, ymin = Low, ymax = High, x = fct_rev(Term), col = color)) +
    geom_errorbar(size = bar_size, width = 0.3) +
    geom_hline(yintercept = 1, size = line_size_blacks) +
    scale_color_manual(values = color) +
    scale_y_continuous(limits = limits, breaks = breaks) +
    ylab("Odds ratio") +
    xlab(NULL) +
    ggtitle(title) +
    geom_point(size = odds_ratio_point_size) +
    coord_flip() +
    theme_bw(base_size = base_size) +
    theme(axis.text = element_text(color = "black"), 
          panel.grid = element_blank(),  
          legend.position = "none",
          plot.title = element_text(hjust = 0.5))
}

plot_compare_effects <- function(ewas, ewas_no_corr, title, limits, breaks, col, point_size, base_size) {
  
  ewas_no_corr_sign_probes <- ewas_no_corr %>% 
    filter(P_FDR < 0.05) %>% 
    pull(Probe)
  
  ewas_sign_probes <- ewas %>% 
    filter(P_FDR < 0.05) %>% 
    pull(Probe)
  
  probes <- unique(c(ewas_no_corr_sign_probes, ewas_sign_probes))
  
  plt_frame_with_correction <- left_join(tibble(Probe = probes), ewas)
  plt_frame_without_correction <- left_join(tibble(Probe = probes), ewas_no_corr)

  
  plt_frame <- tibble(Probe = plt_frame_with_correction$Probe,
                      With_corr = plt_frame_with_correction$Estimate,
                      Without_corr = plt_frame_without_correction$Estimate)
  
  # plt_frame$Mediation <- (abs(plt_frame$With_corr - plt_frame$Without_corr) > 0.2) & (abs(plt_frame$With_corr) < 0.8)
  
  
  plt_frame %>% 
    ggplot(aes(x = Without_corr, y = With_corr, alpha = 0.8)) +
    geom_point(size = point_size, color = col) +
    scale_x_continuous(breaks = breaks, limits = limits) +
    scale_y_continuous(breaks = breaks, limits = limits) +
    ggtitle(title) +
    geom_abline(intercept = 0, slope = 1, color = "black", size = line_size_blacks) +
    xlab("Total effect") +
    ylab("Direct effect") +
    theme_classic(base_size = base_size) +
    theme(axis.text = element_text(color = "black"), 
          legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  
}



# plot_R2_violin2 <- function(tib, title, base_size, tag, no_y_axis) {
# 
#   plt1 <- tib %>% 
#     ggplot(aes(x = Predictor, total_effect, y = R2, fill = Predictor)) + 
#     geom_violin() +   
#     scale_y_continuous(breaks = c(0, 5, 10), limits = c(0, 10)) +
#     scale_fill_manual(values = c(cell_col, gen_col, intrinsic_col, exposures_col), drop = FALSE) +
#     ylab(NULL) +
#     xlab(NULL) +
#     theme_classic(base_size = base_size) +
#     theme(panel.grid.major.y = element_line(color = "grey80"),
#           plot.tag = element_text(size = 16),
#           plot.title = element_text(hjust = 0.5),
#           text = element_text(color = "black"),
#           axis.text = element_text(color = "black"),
#           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#           panel.grid = element_blank(),
#           legend.position = "none")
#   
#   if (no_y_axis) {
#     plt1 <- plt1 + 
#       theme(axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank())
#   }
#   
#   plt2 <- tib %>% 
#     ggplot(aes(x = Predictor, total_effect, y = R2, color = Predictor)) +
#     geom_jitter(width = 0.04, size = point_size_scatter) +
#     scale_y_continuous(breaks = seq(10, 100, 20), limits = c(10, 100)) +  
#     scale_color_manual(values = c(cell_col, gen_col, intrinsic_col, exposures_col), drop = FALSE) +
#     ggtitle(title) +
#     labs(tag = tag) +
#     ylab(NULL) +
#     xlab(NULL) +
#     theme_classic(base_size = base_size) +
#     theme(panel.grid.major.y = element_line(color = "grey80"),
#           plot.title = element_text(hjust = 0.5),
#           text = element_text(color = "black"),
#           axis.text = element_text(color = "black"),
#           axis.text.x = element_blank(),
#           axis.line.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           panel.grid = element_blank(),
#           legend.position = "none")
#   
# 
#   
#   if (no_y_axis) {
#     plt2 <- plt2 + 
#       theme(axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank())
#   }
#   
#   plt2 / plt1 + plot_layout(heights = c(2, 3))
#   
# }

plot_R2_violin2 <- function(tib, title, base_size, tag, no_y_axis) {
  
  plt1 <- tib %>% 
    ggplot(aes(x = Predictor, y = R2, fill = Predictor)) + 
    geom_violin() +
    # geom_point(aes(x = Predictor, y = R2), ~ filter(., Is_outlier), size = point_size_scatter, alpha = 0.8) +
    # geom_boxplot(outlier.shape = NA, alpha = 0.05) +
    # scale_y_continuous(breaks = c(0, 5, 10, 15, 20), limits = c(0, 20)) +
    scale_y_continuous(trans = "log1p", breaks = seq(10, 100, 10)) +
    scale_fill_manual(values = c(col_cell, col_gen, col_intrinsic, col_exposures), drop = FALSE) +
    ylab(NULL) +
    xlab(NULL) +
    theme_classic(base_size = base_size) +
    theme(panel.grid.major.y = element_line(color = "grey80"),
          plot.tag = element_text(size = 16),
          plot.title = element_text(hjust = 0.5),
          text = element_text(color = "black"),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid = element_blank(),
          legend.position = "none")
  
  plt1
  
  # if (no_y_axis) {
  #   plt1 <- plt1 +
  #     theme(axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank())
  # }
  # 
  # plt2 <- tib %>%
  #   ggplot(aes(x = Predictor, total_effect, y = R2, color = Predictor)) +
  #   geom_jitter(width = 0.05, size = point_size_scatter) +
  #   scale_y_continuous(breaks = seq(20, 100, 20), limits = c(20, 100)) +
  #   scale_color_manual(values = c(cell_col, gen_col, intrinsic_col, exposures_col), drop = FALSE) +
  #   ggtitle(title) +
  #   ylab(NULL) +
  #   xlab(NULL) +
  #   theme_classic(base_size = base_size) +
  #   theme(panel.grid.major.y = element_line(color = "grey80"),
  #         plot.title = element_text(hjust = 0.5),
  #         text = element_text(color = "black"),
  #         axis.text = element_text(color = "black"),
  #         axis.text.x = element_blank(),
  #         axis.line.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         panel.grid = element_blank(),
  #         legend.position = "none")
  # 
  # 
  # 
  # if (no_y_axis) {
  #   plt2 <- plt2 + 
  #     theme(axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank())
  # }
  # 
  # plt2 / plt1 + plot_layout(heights = c(2, 3))
  
}




plot_compare_p_values2 <- function(ewas, ewas_no_corr, labels = NULL,
                                  label_size, x_lim = c(-2.5, 20), y_lim = c(-2.5, 50),
                                  breaks_x, breaks_y, col, prune = TRUE, line_size = 1) {
  
  form <- function(x) format(round(100 * x, digits = 2), big.mark = ",")
  if (is_null(labels)) labels <- tibble(x = c(-1.2, 8, 18, 4, -1.2), y = c(-1.5, -1.5, 5, 40, 40))
  
  label_size = label_size * (1 / ggplot2:::.pt)
  
  tib <- dplyr::select(ewas, Probe, P_FDR) %>% 
    dplyr::rename(Corr_P = P_FDR)
  
  tib_no_corr <- dplyr::select(ewas_no_corr, Probe, P_FDR) %>% 
    dplyr::rename(No_corr_P = P_FDR)
  
  tib <- inner_join(tib, tib_no_corr)

  
  p_not_sign <- sum(tib$Corr_P > 0.05 & tib$No_corr_P > 0.05) / nrow(ewas)
  p_not_sign_corr <- sum(tib$Corr_P > 0.05 & tib$No_corr_P < 0.05)  / nrow(ewas)
  p_not_sign_no_corr <- sum(tib$Corr_P < 0.05 & tib$No_corr_P > 0.05) / nrow(ewas)
  p_below_line <- sum(tib$Corr_P < 0.05 & tib$No_corr_P < 0.05 & (tib$No_corr_P > tib$Corr_P)) / nrow(ewas)
  p_above_line <- sum(tib$Corr_P < 0.05 & tib$No_corr_P < 0.05 & (tib$No_corr_P < tib$Corr_P)) / nrow(ewas)
  
  labels$label <- c(paste0("% MI = ", form(p_not_sign)),
                    paste0("% MI = ", form(p_not_sign_no_corr)),
                    paste0("% MI = ", form(p_below_line)),
                    paste0("% MI = ", form(p_above_line)),
                    paste0("% MI = ", form(p_not_sign_corr)))
  
  # tib_sign <- tib %>% filter(No_corr_P < 0.05 | Corr_P < 0.05)
  if (prune) tib_plt <- sample_n(tib, 100000, weight = 10 * pmax(-log10(No_corr_P + 1.69e-313), -2 * log10(Corr_P + 1.69e-313)))
  
  ggplot() +
    geom_point(aes(x = -log10(Corr_P), y = -log10(No_corr_P)), data = tib_plt, size = point_size_scatter, color = col) + 
    # geom_point(aes(x = -log10(Corr_P), y = -log10(No_corr_P)), data = tib_in_catalog, size = point_size_scatter, color = front_col) + 
    scale_x_continuous(breaks = breaks_x, limits = x_lim) +
    scale_y_continuous(breaks = breaks_y, limits = y_lim) +
    geom_line(aes(x = x, y = y), 
              data = tibble(y = c(-log10(0.05), min(x_lim[2], y_lim[2])), x = c(-log10(0.05), min(x_lim[2], y_lim[2]))), 
              size = line_size) +
    # geom_label(aes(x = x, y = y, label = label), data = labels, size = label_size, label.padding = unit(0.15, "lines")) +
    geom_hline(yintercept = -log10(0.05), color = "black", size = line_size) + 
    geom_vline(xintercept = -log10(0.05), color = "black", size = line_size) +
    xlab(TeX("Cell adjusted -log_{10} (P_{adj})")) +
    ylab(TeX("No adjustment -log_{10} (P_{adj})"))
}

plot_compare_p_values3 <- function(ewas, ewas_no_corr, x_lim = c(-2.5, 20), y_lim = c(-2.5, 50),
                                   breaks_x, breaks_y, col, prune = TRUE, line_size = 1) {
  
  form <- function(x) format(round(100 * x, digits = 2), big.mark = ",")

  tib <- dplyr::select(ewas, Probe, P_FDR) %>% 
    dplyr::rename(Corr_P = P_FDR)
  
  tib_no_corr <- dplyr::select(ewas_no_corr, Probe, P_FDR) %>% 
    dplyr::rename(No_corr_P = P_FDR)
  
  tib <- inner_join(tib, tib_no_corr)
  
  # tib_sign <- tib %>% filter(No_corr_P < 0.05 | Corr_P < 0.05)
  if (prune) tib_plt <- sample_n(tib, 100000, weight = 10 * pmax(-log10(No_corr_P + 1.69e-313), -2 * log10(Corr_P + 1.69e-313)))
  
  ggplot() +
    geom_point(aes(x = -log10(Corr_P), y = -log10(No_corr_P)), data = tib_plt, size = point_size_scatter, color = col) + 
    # geom_point(aes(x = -log10(Corr_P), y = -log10(No_corr_P)), data = tib_in_catalog, size = point_size_scatter, color = front_col) + 
    scale_x_continuous(breaks = breaks_x, limits = x_lim) +
    scale_y_continuous(breaks = breaks_y, limits = y_lim) +
    geom_line(aes(x = x, y = y), 
              data = tibble(y = c(-log10(0.05), min(x_lim[2], y_lim[2])), x = c(-log10(0.05), min(x_lim[2], y_lim[2]))), 
              size = line_size) +
    # geom_label(aes(x = x, y = y, label = label), data = labels, size = label_size, label.padding = unit(0.15, "lines")) +
    geom_hline(yintercept = -log10(0.05), color = "black", size = line_size) + 
    geom_vline(xintercept = -log10(0.05), color = "black", size = line_size) +
    xlab(TeX("Direct effect, -log_{10} (P_{adj})")) +
    ylab(TeX("Total effect, -log_{10} (P_{adj})"))
}

plot_effect_distribution_over_state <- function(tib, title_main = "Effect size distribution", base_size = 12, line_size, outlier_point_size, fontsize_count = 12, states, color) {
  
  fontsize_count <- fontsize_count / ggplot2:::.pt 
  
  tib_count <- tib %>% 
    group_by(.data[[states]]) %>% 
    summarize(N = n(), Pos = sum(Estimate > 0)) %>%
    mutate(Pos_ratio = round(Pos / N, digits = 2),
           N = format(N, big.mark = ","))
  
  # tib_count <- tib %>% 
  #   group_by(.data[[states]]) %>% 
  #   summarize(Neg = sum(Estimate < 0), Pos = sum(Estimate > 0)) %>%
  #   mutate(Pos_ratio = round(Pos / (Pos + Neg), digits = 2),
  #          Neg = format(Neg, big.mark = ","),
  #          Pos = format(Pos, big.mark = ","))  
  
  # plt_count_neg <- ggplot(tib_count, aes(x = fct_rev(.data[[states]]), y = 0, label = Neg)) + 
  #   geom_text(color = color) +
  #   ggtitle("# Neg.") +
  #   coord_flip() +
  #   xlab(NULL) +
  #   ylab(NULL) +
  #   theme_classic() +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         axis.text = element_text(color = "black"), 
  #         axis.text.x = element_blank(),
  #         axis.line.x = element_blank(),
  #         axis.ticks.x = element_blank())
  
  plt_count <- ggplot(tib_count, aes(x = fct_rev(.data[[states]]), y = 0, label = N)) + 
    geom_text(color = color, size = fontsize_count) +
    ggtitle("#") +
    coord_flip() +
    xlab(NULL) +
    ylab(NULL) +
    theme_classic(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank(),  
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  
  plt_ratio_pos <- ggplot(tib_count, aes(x = fct_rev(.data[[states]]), y = 0, label = Pos_ratio, color = Pos_ratio)) + 
    geom_text(size = fontsize_count) +
    ggtitle("% Pos.") +
    scale_color_gradient(low = blues[3], high = color) +
    coord_flip() +
    xlab(NULL) +
    ylab(NULL) +
    theme_classic(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank(),  
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")
  
  
  plt_effect_distribution_over_regions <- ggplot(tib, aes(x = fct_rev(.data[[states]]), y = Estimate)) +
    geom_violin(color = color, size = line_size) +
    ggtitle(title_main) +
    geom_point(data = ~ filter(., Is_outlier), color = color, size = outlier_point_size) +
    geom_hline(yintercept = 0, size = line_size) +
    coord_flip() +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw(base_size = base_size) + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.x = element_blank(),
          axis.text = element_text(color = "black"))
  
  plt_effect_distribution_over_regions + plt_count + plt_ratio_pos + plot_layout(widths = c(8, 1, 1))
  
}


# plot_compare_p_values <- function(ewas, ewas_no_corr, compare_data, labels = NULL,
#                                   label_size, x_lim = c(-2.5, 20), y_lim = c(-2.5, 50),
#                                   breaks_x, breaks_y, back_col, front_col, prune = TRUE) {
#   
#   form <- function(x) format(round(100 * x, digits = 2), big.mark = ",")
#   if (is_null(labels)) labels <- tibble(x = c(-1.2, 8, 18, 4, -1.2), y = c(-1.5, -1.5, 5, 40, 40))
#   
#   label_size = label_size * (1 / ggplot2:::.pt)
#   
#   tib <- select(ewas, Probe, P_FDR) %>% 
#     rename(Corr_P = P_FDR)
#   
#   tib_no_corr <- select(ewas_no_corr, Probe, P_FDR) %>% 
#     rename(No_corr_P = P_FDR)
#   
#   tib <- inner_join(tib, tib_no_corr)
#   
#   tib_in_catalog <- tib %>% 
#     filter(Probe %in% unique(compare_data$Probe))
#   
#   p_not_sign <- sum(tib$Corr_P > 0.05 & tib$No_corr_P > 0.05) / nrow(ewas)
#   p_not_sign_corr <- sum(tib$Corr_P > 0.05 & tib$No_corr_P < 0.05)  / nrow(ewas)
#   p_not_sign_no_corr <- sum(tib$Corr_P < 0.05 & tib$No_corr_P > 0.05) / nrow(ewas)
#   p_below_line <- sum(tib$Corr_P < 0.05 & tib$No_corr_P < 0.05 & (tib$No_corr_P > tib$Corr_P)) / nrow(ewas)
#   p_above_line <- sum(tib$Corr_P < 0.05 & tib$No_corr_P < 0.05 & (tib$No_corr_P < tib$Corr_P)) / nrow(ewas)
#   
#   p_cat_not_sign <- sum(tib_in_catalog$Corr_P > 0.05 & tib_in_catalog$No_corr_P > 0.05) / nrow(compare_data)
#   p_cat_not_sign_corr <- sum(tib_in_catalog$Corr_P > 0.05 & tib_in_catalog$No_corr_P < 0.05)  / nrow(compare_data)
#   p_cat_not_sign_no_corr <- sum(tib_in_catalog$Corr_P < 0.05 & tib_in_catalog$No_corr_P > 0.05) / nrow(compare_data)
#   p_cat_below_line <- sum(tib_in_catalog$Corr_P < 0.05 & tib_in_catalog$No_corr_P < 0.05 & (tib_in_catalog$No_corr_P > tib_in_catalog$Corr_P)) / nrow(compare_data)
#   p_cat_above_line <- sum(tib_in_catalog$Corr_P < 0.05 & tib_in_catalog$No_corr_P < 0.05 & (tib_in_catalog$No_corr_P < tib_in_catalog$Corr_P)) / nrow(compare_data)
#   
#   labels$label <- c(paste0("% MI = ", form(p_not_sign), "\n % Cat. = ", form(p_cat_not_sign)),
#                     paste0("% MI = ", form(p_not_sign_no_corr), "\n % Cat. = ", form(p_cat_not_sign_no_corr)),
#                     paste0("% MI = ", form(p_below_line), "\n % Cat. = ", form(p_cat_below_line)),
#                     paste0("% MI = ", form(p_above_line), "\n % Cat. = ", form(p_cat_above_line)),
#                     paste0("% MI = ", form(p_not_sign_corr), "\n % Cat. = ", form(p_cat_not_sign_corr)))
#   
#   # tib_sign <- tib %>% filter(No_corr_P < 0.05 | Corr_P < 0.05)
#   if (prune) tib_plt <- sample_n(tib, 100000, weight = 10 * pmax(-log10(No_corr_P + 1.69e-313), -2 * log10(Corr_P + 1.69e-313)))
#   
#   ggplot() +
#     geom_point(aes(x = -log10(Corr_P), y = -log10(No_corr_P)), data = tib_plt, size = point_size_scatter, color = back_col) + 
#     # geom_point(aes(x = -log10(Corr_P), y = -log10(No_corr_P)), data = tib_in_catalog, size = point_size_scatter, color = front_col) + 
#     scale_x_continuous(breaks = breaks_x, limits = x_lim) +
#     scale_y_continuous(breaks = breaks_y, limits = y_lim) +
#     geom_line(aes(x = x, y = y), 
#               data = tibble(y = c(-log10(0.05), min(x_lim[2], y_lim[2])), x = c(-log10(0.05), min(x_lim[2], y_lim[2]))), 
#               size = line_size_scatter) +
#     geom_label(aes(x = x, y = y, label = label), data = labels, size = label_size, label.padding = unit(0.15, "lines")) +
#     geom_hline(yintercept = -log10(0.05), color = "black", size = 1) + 
#     geom_vline(xintercept = -log10(0.05), color = "black", size = 1) +
#     xlab(TeX("Cell adjusted -log_{10} (P_{adj})")) +
#     ylab(TeX("No adjustment -log_{10} (P_{adj})")) +
#     theme_classic()
# }

# plot_R2_violin <- function(tib, title, base_size) {
#   tib %>% 
#     ggplot(aes(x = Predictor, y = R2, fill = Predictor)) + 
#     geom_violin() +   
#     geom_boxplot(fill = "#FFFFFF88", notch = TRUE, width = 0.2, alpha = 0.5, outlier.alpha = 0.5, outlier.size = point_size_scatter) +
#     scale_y_continuous(trans = "log1p", breaks = c(0, 1, 4, 10, 20, 40, 60, 80, 100), limits = c(0, 100)) +
#     scale_fill_manual(values = c(cell_col, gen_col, intrinsic_col, exposures_col)) +
#     ggtitle(title) +
#     ylab(NULL) +
#     xlab(NULL) +
#     theme_classic(base_size = base_size) +
#     theme(plot.title = element_text(hjust = 0.5),
#           text = element_text(color = "black"),
#           axis.text = element_text(color = "black"),
#           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#           panel.grid = element_blank(),
#           legend.position = "none")
# }