library(tidyverse)
library(patchwork)

load("./Data/RData/composition_simulation_objects.RData")
plt_frame <- bind_rows(add_column(si_mye_tib, Changed = "Myeloid"),
                       add_column(si_CD4_tib, Changed = "CD4+ T cells"),
                       add_column(si_CD4_diff_tib, Changed = "Diff CD4+ T cells"))
                       
                       
plt2 <- si_mye_tib %>% 
  mutate(term = factor(term, unique(term))) %>% 
  ggplot(aes(x = P_value)) + 
  geom_histogram(fill = "skyblue", colour = "black", bins = 50) + 
  xlab("P value") +
  ylab("Count") +
  facet_wrap(vars(term), ncol = 5) + 
  theme_classic() +
  theme(strip.background = element_rect(fill = "white"))  +
  labs(tag = "B")


plt3 <- si_CD4_tib %>% 
  mutate(term = factor(term, unique(term))) %>% 
  ggplot(aes(x = P_value)) + 
  geom_histogram(fill = "skyblue", colour = "black", bins = 50) + 
  xlab("P value") +
  ylab("Count") +
  facet_wrap(vars(term), ncol = 5) + 
  theme_classic() +
  theme(strip.background = element_rect(fill = "white")) +
  labs(tag = "C")


plt4 <- si_CD4_diff_tib %>% 
  mutate(term = factor(term, unique(term))) %>% 
  ggplot(aes(x = P_value)) + 
  geom_histogram(fill = "skyblue", colour = "black", bins = 50) + 
  xlab("P value") +
  ylab("Count") +
  facet_wrap(vars(term), ncol = 5) + 
  theme_classic() +
  theme(strip.background = element_rect(fill = "white")) +
  labs(tag = "D")

alpha <- sort(exp(rnorm(16, 0, 1)))
alpha[cell_order] <- alpha


X <- as.data.frame(rdirichlet(N, alpha = alpha))
plt1 <- plot_list(X) +
  labs(tag = "A")




plt <- plt1 / (plt2 / plt3 / plt4)
ggsave("./Plots/Revision2_plots/Supp_fig9/sup_fig_9.pdf", plt, width = 200, height = 300, units = "mm")