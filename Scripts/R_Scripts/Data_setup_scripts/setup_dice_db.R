library("tidyverse")
library("AnnotationDbi")
library("org.Hs.eg.db")
dice <- read_csv("./Data/Dice_database/mean_tpm_merged.csv")
ensembl_to_symbol <- select(org.Hs.eg.db, 
                            keys = dice$gene,  
                            keytype = "ENSEMBL",
                            columns = c("ENSEMBL", "SYMBOL")) %>% 
  as_tibble() %>% 
  filter(!is.na(SYMBOL)) %>% 
  dplyr::rename(gene = ENSEMBL, Gene = SYMBOL)

dice <- left_join(ensembl_to_symbol, dice) %>% dplyr::select(-gene)
saveRDS(dice, "./Data/Dice_database/dice.rds")