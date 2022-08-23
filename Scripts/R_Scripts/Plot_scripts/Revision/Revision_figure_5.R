
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(patchwork)
library(scales)
library(latex2exp)
source("./Scripts/R_scripts/Libraries/functions_for_plotting.R")

# Load data ---------------------------------------------------------------

vars <- c("Total_effect_intrinsic", "Total_effect_exposures", "Total_effect_cells", "Total_effect_genetic")
prop_var <- readRDS("./Data/RData/Results/Proportion_of_variance_explained/prop_var_explained_genetic_interactions_averaged_884.rds")

total_effect <- prop_var %>% 
  select(-contains("Houseman"), -contains("IDOL")) %>% 
  select(Probe | starts_with("Total") | Full_effect) %>%
  pivot_longer(!Probe, names_to = "Predictor", values_to = "R2") %>% 
  filter(Predictor != "Total_effect_panel5") %>% 
  mutate(Predictor = fct_recode(Predictor,
                                `Full model` = "Full_effect",
                                Intrinsic = "Total_effect_intrinsic", 
                                Exposures = "Total_effect_exposures", 
                                `16 Cells` = "Total_effect_cells",
                                Genetic = "Total_effect_genetic")) %>% 
  group_by(Predictor) %>% 
  mutate(Is_outlier = is_outlier(R2)) %>% 
  ungroup()

total_effect_probe <- total_effect %>% 
  filter(!Predictor %in% c("6 Cells", "Full model")) %>% 
  group_by(Probe) %>% 
  filter(R2 == R2[which.max(R2)]) %>% 
  ungroup()

direct_effect <- prop_var %>%
  select(Probe | starts_with("Conditional")) %>%
  pivot_longer(!Probe, names_to = "Predictor", values_to = "R2") %>%
  mutate(Predictor = fct_recode(Predictor, 
                                Intrinsic = "Conditional_effect_intrinsic", 
                                Exposures = "Conditional_effect_exposures", 
                                Genetic = "Conditional_effect_genetic",
                                `Genetic int.` = "Conditional_effect_genetic_interaction")) %>% 
  group_by(Predictor) %>% 
  mutate(Is_outlier = is_outlier(R2)) %>% 
  ungroup()

direct_effect_probe <- direct_effect %>%
  group_by(Probe) %>%
  summarize(Top_pred = Predictor[which.max(R2)],
            R2_sum = sum(R2)) %>% 
  ungroup()

# A -----------------------------------------------------------------------

compute_quantiles <- function(x) {
  p <- ecdf(x)
  p(x)
}


plt_frame <- total_effect %>%
  mutate(Predictor = fct_recode(Predictor, `16 cell comp.` = "16 Cells")) %>% 
  mutate(Predictor = fct_relevel(Predictor, "Full model", "16 cell comp.", "Genetic", "Intrinsic", "Exposures")) %>% 
  group_by(Predictor) %>% 
  mutate(y = compute_quantiles(R2)) %>% 
  ungroup()

pltA <- ggplot(plt_frame, aes(x = R2, y = 100 * (1 - y), color = Predictor)) + 
  geom_line(size = line_size_scatter) + 
  xlab(TeX("Out of sample $R^2$")) + 
  ylab(TeX("% sites with $>R^2$")) +
  scale_x_continuous(limits = c(0, 100), breaks = seq.int(0, 100, 10)) +
  scale_color_manual(values = c("black", col_cell, col_gen, col_gen, col_intrinsic, col_exposures), drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), breaks = seq.int(0, 100, 10)) +
  theme_bw(base_size = font_size_scatter) +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.8, 0.8),
        legend.key.size = unit(0.1, "lines")) +
  labs(tag = "A")

ggsave("./Plots/Revision_plots/Figure5/pltA_total_effects.pdf", pltA, width = 6, height = 6, units = "cm")

# B -----------------------------------------------------------------------

direct_effect$R2[direct_effect$Predictor == "Genetic int." & direct_effect$R2 < 0] <- 0

plt_frame <- direct_effect %>%
  mutate(Predictor = fct_relevel(Predictor, "Genetic", "Intrinsic", "Exposures", "Genetic int.")) %>% 
  group_by(Predictor) %>% 
  mutate(y = compute_quantiles(R2)) %>% 
  ungroup()

pltB <- ggplot(plt_frame, aes(x = R2, y = 100 * (1 - y), color = Predictor)) + 
  geom_line(size = line_size_scatter) + 
  xlab(TeX("Out of sample $R^2$")) + 
  ylab(TeX("% sites with $>R^2$")) +
  scale_x_continuous(limits = c(0, 100), breaks = seq.int(0, 100, 10)) +
  scale_color_manual(values = c(col_gen, col_gen, col_intrinsic, col_exposures), drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), breaks = seq.int(0, 100, 10)) +
  theme_bw(base_size = font_size_scatter) +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.8, 0.8),
        legend.key.size = unit(0.1, "lines")) +
  labs(tag = "B")

       
ggsave("./Plots/Revision_plots/Figure5/pltB_direct_effects.pdf", pltB, width = 6, height = 6, units = "cm")


# C -----------------------------------------------------------------------

nr_sites <- 20000
mid_label <- round(nr_sites / 2)

probe_selection <- total_effect_probe %>% 
  group_by(Predictor) %>% 
  arrange(desc(R2), .by_group = TRUE) %>% 
  slice(1:nr_sites) %>% 
  ungroup()

base <- probe_selection %>%
  summarize(start = nth(Probe, 1), 
            title = nth(Probe, mid_label), 
            end = nth(Probe, nr_sites))

base <- tibble(Predictor = unique(probe_selection$Predictor))

base$start <- c(nth(pull(filter(probe_selection, Predictor == base$Predictor[[1]]), Probe), 1),
                nth(pull(filter(probe_selection, Predictor == base$Predictor[[2]]), Probe), 1),
                nth(pull(filter(probe_selection, Predictor == base$Predictor[[3]]), Probe), 1),
                nth(pull(filter(probe_selection, Predictor == base$Predictor[[4]]), Probe), 1))


base$title <- c(nth(pull(filter(probe_selection, Predictor == base$Predictor[[1]]), Probe), round(0.75 * nr_sites)),
                nth(pull(filter(probe_selection, Predictor == base$Predictor[[2]]), Probe), mid_label),
                nth(pull(filter(probe_selection, Predictor == base$Predictor[[3]]), Probe), mid_label),
                nth(pull(filter(probe_selection, Predictor == base$Predictor[[4]]), Probe), round(0.35 * nr_sites)))


base$end <- c(nth(pull(filter(probe_selection, Predictor == base$Predictor[[1]]), Probe), nr_sites),
              nth(pull(filter(probe_selection, Predictor == base$Predictor[[2]]), Probe), nr_sites),
              nth(pull(filter(probe_selection, Predictor == base$Predictor[[3]]), Probe), nr_sites),
              nth(pull(filter(probe_selection, Predictor == base$Predictor[[4]]), Probe), nr_sites))


plt_frame <-  total_effect %>% right_join(probe_selection)
plt_frame$Probe <- factor(plt_frame$Probe, probe_selection$Probe)

pltC <- ggplot(plt_frame, aes(x = Probe, y = R2, fill = Predictor)) +      
  
  geom_bar(stat = "identity", width = 1.1)  +
  geom_text(data =  base, 
            aes(x = title, y = -10, label = Predictor), 
            hjust = c(1, 1, 0, 0), 
            colour = "black", 
            alpha = 0.8, 
            size = font_size_label, 
            fontface="bold", 
            inherit.aes = FALSE) +
  annotate("text", 
           x = -1, 
           y = c(0, 25, 50, 75, 100), 
           label = c("0", "25", "50", "75", "100") , 
           color = "black", 
           size = font_size_label, 
           angle = 0, 
           fontface = "bold", 
           hjust = 1) +
  annotate("text", 
           x = 40000 - 100, 
           y = c(0, 25, 50, 75, 100), 
           label = c("0", "25", "50", "75", "100") , 
           color = "black", 
           size = font_size_label, 
           angle = 0, 
           fontface = "bold", 
           hjust = 1) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(-100, 100), minor_breaks = NULL) +
  scale_fill_manual(values = c(col_cell, col_exposures, col_gen, col_intrinsic)) +
  theme_minimal(base_size = font_size_scatter) + 
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = line_size_scatter), 
        plot.margin = unit(rep(0, 4), "cm")) +
  coord_polar(theta = "x", clip = "on") +
  labs(tag = "C")



ggsave("./Plots/Revision_plots/Figure5/pltC_circle_plot.pdf", pltC, width = 6, height = 6, units = "cm")


# Assemble plot -----------------------------------------------------------


design <- "
AB
CC
CC
"


plt <- pltA + pltB + pltC + plot_layout(design = design) &
  theme(plot.tag = element_text(size = 6, face = "bold"))

ggsave("./Plots/Revision_plots/Figure5/Figure_5.png", plt, height = 150, width = 120, units = "mm", dpi = 500)
ggsave("./Plots/Revision_plots/Figure5/Figure_5.pdf", plt, height = 150, width = 120, units = "mm")


