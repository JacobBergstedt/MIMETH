.libPaths("/pasteur/zeus/projets/p02/IGSR/R_libraries/4.1.0") 
library(CMAverse)
library(tidyverse)
library(parallel)
ewas_med_age <- readRDS("./Data/Chunk_data/Results/Cell_mediation/m_values_884/Age.rds") %>% 
  filter()
cell_list_15_cells <- scan("./Data/prop_controls_15.txt", what = character())
covs <- get_covs_884()
cells <- covs[cell_list_15_cells]
names(cells) <- paste0("M", 1:15)
exposure <- covs[["Age"]] / 50


ewas_med_age <- readRDS("./Data/Chunk_data/Results/Cell_mediation/m_values_884/Age.rds") %>%
  filter(Effect == "Mediation") %>% 
  arrange(P_value)

cpgs <- ewas_med_age$Probe[1:10]
meth <- get_m_values()
meth <- meth[match(covs$SUBJID, meth$SUBJID),]

covs <- select(covs, CMV_serostatus, Sex, Ancestry_PC1, Ancestry_PC2, Smoking_status)
db <- cbind(Exposure = exposure, cells, covs)

fit_CMAverse <- function(y) {
  
  db <- cbind(Y = y, db)
  
  m <- cmest(data = db, 
        model = "rb", 
        outcome = "Y", 
        exposure = "Exposure",
        basec = names(covs),
        mediator = paste0("M", 1:15), 
        EMint = FALSE,
        mreg = as.list(rep("linear", 15)), 
        yreg = "linear",
        astar = 0, 
        a = 1, 
        mval = as.list(rep(0, 15)), 
        estimation = "imputation", 
        inference = "bootstrap",
        nboot = 500)
  
  summary(m)$summarydf
  
}

compare_results <- function(cpg, med, fit_sum) {
  
  med <- med %>% filter(Probe == cpg)
  tibble(CpG = cpg, Our_est =  med$Estimate, Our_sd = med$Standard_error, Their_est = fit_sum["tnie", "Estimate"], Their_sd = fit_sum["tnie", "Std.error"])
  
}

p <- mclapply(meth[cpgs], fit_CMAverse, mc.cores = 10)
o <- map_dfr(cpgs, ~ compare_results(., ewas_med_age, p[[.]]))

saveRDS(o, "./Data/RData/CMAVerse_comparison.rds")