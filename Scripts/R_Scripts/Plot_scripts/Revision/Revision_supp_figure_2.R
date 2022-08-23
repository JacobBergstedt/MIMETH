library(tidyverse)
library(ggcorrplot)
library(corrplot)
library(scales)
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
source("./Scripts/R_scripts/Libraries/functions_for_plotting.R")
source("./Scripts/MiscFunctions.R")
covs <- get_covs_884() %>% 
  mutate(X_CD8_mem = X_CD8_CM_of_total.panel1 + X_CD8_EM_of_total.panel1 + X_CD8_EMRA_of_total.panel1,
         X_CD4_mem = X_CD4_CM_of_total.panel1 + X_CD4_EM_of_total.panel1 + X_CD4_EMRA_of_total.panel1)
panel5 <- get_panel5_cells()
str <- names(covs)
idol_extended_labels <- str[grepl("IDOL_extended", str)]
idol_labels <- str[grepl("IDOL", str) & !(grepl("IDOL_extended", str))]
houseman_labels <- str[grepl("Houseman", str)]
panel5_labels <- names(panel5[-1])
FACS_labels <- names(select(covs, starts_with("X_")))

get_plt_frame <- function(x, keys, labels) {
  
  
  p1 <- x[c("SUBJID", keys)] %>% 
    pivot_longer(cols = -SUBJID) %>% 
    mutate(name = fct_recode(factor(name, keys), !!!keys)) %>% 
    rename(Decomp = value)
  
  p2 <- x[c("SUBJID", names(keys))] %>% 
    pivot_longer(cols = -SUBJID) %>% 
    rename(MI = value)
  
  inner_join(p1, p2) %>% 
    left_join(tibble(name = names(keys), Cell = labels)) %>% 
    select(-name)
  
  
}

# IDOL extended -----------------------------------------------------------


keys <- idol_extended_labels[!idol_extended_labels %in% c("IDOL_extended_Bmem", "IDOL_extended_Bnv", "IDOL_extended_Treg")]
names(keys) <- c("X_VIABLE_BASOPHILS_OF_TOTAL.panel7",
                 "X_CD4_mem",
                 "X_CD4_naive_of_total.panel1",
                 "X_CD8_mem",
                 "X_CD8_naive_of_total.panel1",
                 "X_VIABLE_EOSINOPHILS_OF_TOTAL.panel7",
                 "X_mono_of_total.panel5",
                 "X_VIABLE_NEUTROPHILS_OF_TOTAL.panel7",
                 "X_NK_of_total.panel4")
labels <- c("basophils", "CD4 memory", "CD4 naive", "CD8 memory", "CD8 naive", "eosinophils", "monocytes", "neutrophils", "NK")
pltA <- get_plt_frame(covs, keys, labels) %>% 
  ggplot(aes(x = Decomp, y = MI)) +
  ylab("FACS") +
  xlab("IDOL extended") + 
  geom_point(size = point_size_scatter, col = col_cell) +
  facet_wrap(vars(Cell), scales = "free") +
  theme_bw(base_size = font_size_scatter) +
  theme(strip.background = element_rect(fill = "white")) +
  labs(tag = "A")

# IDOL --------------------------------------------------------------------

plt_frame <- panel5 %>% inner_join(covs[c("SUBJID", idol_labels)])
keys <- idol_labels
names(keys) <- c("X_CD8_of_total.panel5",
                 "X_CD4_of_total.panel5",
                 "X_NK_of_total.panel5",
                 "X_CD19pos_of_total.panel5",
                 "X_mono_of_total.panel5",
                 "X_neutrophils_of_total.panel5")
labels <- c("CD8", "CD4", "NK", "B cells", "monocytes", "neutrophils")
pltB <- get_plt_frame(plt_frame, keys, labels) %>% 
  ggplot(aes(x = Decomp, y = MI)) +
  ylab(NULL) +
  xlab("IDOL") + 
  geom_point(size = point_size_scatter, col = col_cell) +
  facet_wrap(vars(Cell), scales = "free") +
  theme_bw(base_size = font_size_scatter) +
  theme(strip.background = element_rect(fill = "white")) +
  labs(tag = "B")



# Houseman ----------------------------------------------------------------


plt_frame <- panel5 %>% inner_join(covs[c("SUBJID", houseman_labels)])
keys <- houseman_labels
names(keys) <- c("X_CD8_of_total.panel5",
                 "X_CD4_of_total.panel5",
                 "X_NK_of_total.panel5",
                 "X_CD19pos_of_total.panel5",
                 "X_mono_of_total.panel5",
                 "X_neutrophils_of_total.panel5")
labels <- c("CD8", "CD4", "NK", "B cells", "monocytes", "neutrophils")
pltC <- get_plt_frame(plt_frame, keys, labels) %>% 
  ggplot(aes(x = Decomp, y = MI)) +
  geom_point(size = point_size_scatter, col = col_cell) +
  ylab(NULL) +
  xlab("Houseman") + 
  facet_wrap(vars(Cell), scales = "free") +
  theme_bw(base_size = font_size_scatter) +
  theme(strip.background = element_rect(fill = "white")) +
  labs(tag = "C")



# Correlation 1 -----------------------------------------------------------

x <- covs %>% inner_join(select(panel5, -X_mono_of_total.panel5, -X_NK_of_total.panel5))
x <- x[c(FACS_labels, panel5_labels[!panel5_labels %in% c("X_neutrophils_of_total.panel5", "X_CD4negCD8neg_of_total.panel5", "X_CD19pos_of_total.panel5", "X_NK_of_total.panel5", "X_mono_of_total.panel5")], idol_extended_labels, idol_labels, houseman_labels)]
names(x) <- c("Neu", "Baso", "Eos", "Mono", "CD4 naive", "CD4 CM", "CD4 EM", "CD4 EMRA", "CD8 naive", "CD8 CM", "CD8 EM", "CD8 EMRA", "CD8 neg CD4 neg", "B", "NK", "DC", "CD8 mem", "CD4 mem", 
              "CD4", "CD8", 
              "IDOL ext bas", "IDOL ext B mem", "IDOL ext B naive", "IDOL ext CD4 mem", "IDOL ext CD4 naive", "IDOL ext CD8 mem", "IDOL ext CD8 naive", "IDOL ext eos", "IDOL ext mono", "IDOL ext neu", "IDOL ext NK", "IDOL ext Treg", 
              "IDOL CD8", "IDOL CD4", "IDOL NK", "IDOL B", "IDOL mono", "IDOL neu", 
              "House CD8", "House CD4", "House NK", "House B", "House Mono", "House neu")

m <- cor(x, method = "spearman")
pltD <- wrap_elements(panel = ~corrplot(m, method = "circle", type = "lower", tl.col="black", tl.srt = 45, tl.cex = 0.5, cl.cex = 0.5), clip = TRUE) +
  labs(tag = "D")



# Correlation 2 -----------------------------------------------------------

pltE <- wrap_elements(panel = ~corrplot(m, method = "circle", order = "hclust", type = "lower", tl.col="black", tl.srt = 45, tl.cex = 0.5, cl.cex = 0.5), clip = TRUE) +
  labs(tag = "E")


# Combine plot ------------------------------------------------------------


design <- "
  AABBCC
  DDDEEE
  DDDEEE
"

plt <- pltA + pltB + pltC + pltD + pltE + plot_layout(design = design) &
  theme(plot.tag = element_text(size = 6, face = "bold"))

ggsave("./Plots/Revision_plots/Supp_fig2/Sup_fig_2.png", plt, width = 180, height = 180, units = "mm", dpi = 500)
ggsave("./Plots/Revision_plots/Supp_fig2/Sup_fig_2.pdf", plt, width = 180, height = 180, units = "mm")




