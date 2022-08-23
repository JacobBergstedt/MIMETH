
library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(parallel)

overlap_with_cpg_sites <- function(tib, meth_anno, lift_over_chain) {
  
  tfbs_range <- GRanges(seqnames = tib$Chromosome,
                        ranges = IRanges(start = tib$Start, end = tib$End))
  
  tfbs_range <- GenomicRanges::reduce(unlist(liftOver(tfbs_range, lift_over_chain)))
  
  meth_range <- GRanges(seqnames = meth_anno$Probe_chr,
                        ranges = IRanges(start = meth_anno$Probe_position, end = meth_anno$Probe_position))
  
  overlaps <- findOverlaps(tfbs_range, meth_range)
  is_TFBS <- seq_len(nrow(meth_anno)) %in% to(overlaps)
  ifelse(is_TFBS, "Yes", "No")
  
}


remap <- read_tsv("./Data/TFBS_2019/Remap/remap2020_all_macs2_hg38_v1_0.bed.gz", 
                  col_names = c("Chromosome", "Start", "End", "TF_Cell", "Score", "Strand", "Unknown1", "Unknown2", "Unknown3")) %>% 
  mutate(TF_Cell = replace_bad_periods(TF_Cell)) %>% 
  separate(col = TF_Cell, into = c("Experiment_ID", "TF", "Cell_line"), sep = "\\.") %>% 
  filter(Chromosome %in% paste0("chr", 1:22))

meth_loc <- get_anno_meth_location(fix = FALSE)
path_to_chain = system.file(package = "liftOver", "extdata", "hg38ToHg19.over.chain")
lift_over_chain = import.chain(path_to_chain)

keys <- remap %>% group_by(TF) %>% group_keys()
meth_anno_TFBS <- remap %>% 
  group_by(TF) %>% 
  group_split() %>% 
  set_names(keys$TF) %>% 
  mclapply(overlap_with_cpg_sites, meth_anno = meth_loc, lift_over_chain = lift_over_chain, mc.cores = 8)

meth_anno_TFBS <- add_column(as_tibble(meth_anno_TFBS), meth_loc["Probe"], .before = 1)
saveRDS(meth_anno_TFBS, "./Data/RData/Methylation/Annotation/meth_anno_TFBS.rds")