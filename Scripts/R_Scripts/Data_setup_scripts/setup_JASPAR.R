library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
tib_to_granges <- function(tib) {
  GRanges(seqnames = tib$TF_chr, 
          ranges = IRanges(start = tib$TF_start, end = tib$TF_end), 
          TF = tib$TF, strand = tib$TF_strand)
}
granges_to_tib <- function(granges) {
  ranges <- ranges(granges)
  tibble(TF_chr = as.character(seqnames(granges)),
         TF_start = start(ranges),
         TF_end = end(ranges),
         TF_strand = as.character(strand(granges)),
         TF = granges$TF)
}
hg38_TFBS <- readRDS("./Data/Jaspar_TFBS/JASPAR_hg38.rds")
hg38_TFBS <- tib_to_granges(hg38_TFBS)
hg38_to_hg19_chain <- import.chain("./Data/hg38ToHg19.over.chain")
hg19_lift_TFBS <- unlist(liftOver(hg38_TFBS, hg38_to_hg19_chain))
hg19_1 <- granges_to_tib(hg19_lift_TFBS)
hg19_2 <- readRDS("./Data/Jaspar_TFBS/JASPAR_hg19.rds") %>% 
  dplyr::select(-JASPAR_code)
TFBS <- bind_rows(hg19_1, hg19_2) %>% 
  tib_to_granges()
meth_loc <- get_meth_location()
meth_loc_ranges <- GRanges(seqnames = paste0("chr", meth_loc$Probe_chr), 
                          ranges = IRanges(start = meth_loc$Probe_position, end = meth_loc$Probe_position))
meth_TFBS_overlap <- findOverlaps(meth_loc_ranges, TFBS)
meth_TFBS <- meth_loc[from(meth_TFBS_overlap), ]
meth_TFBS$TF <- TFBS[to(meth_TFBS_overlap), ]$TF
saveRDS(meth_TFBS, "./Data/RData/Methylation/Annotation/meth_TFBS.rds")