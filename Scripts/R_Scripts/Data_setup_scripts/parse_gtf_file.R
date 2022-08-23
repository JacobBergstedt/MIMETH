p <- read_delim("./Data/hg19.refGene.gtf.gz", delim = "\t", 
                col_names = c("Chromosome", "Source", "Feature", "Start", "Stop", "Score", "Strand", "Phase", "Additional"), escape_double = FALSE) %>% 
  filter(Feature == "exon")

o <- map(str_split(gsub("(\")|(;)", "", p$Additional), " "), ~ .[c(2, 4, 6, 8, 10)])
o <- transpose(o)
o <- as_tibble(set_names(map(o, unlist), c("Gene", "Transcript_ID", "Exon_number", "Exon_ID", "Gene_name")))

p <- p %>% 
  select(-Additional) %>% 
  bind_cols(o)

o <- p %>% filter(Gene == "DNMT3A")

# o <- p %>% filter(Chromosome == "chr2", Start > 25428591, Stop < 25592023) %>% arrange(Start)

add_introns <- function(db, keys) {
  db_introns <- db[-1, ]
  db_introns$Stop <- tail(db$Start, -1)
  db_introns$Start <- head(db$Stop, -1)
  db_introns$Feature <- "intron"
  bind_rows(db, db_introns)
}

u <- o %>% group_by(Gene, Transcript_ID) %>% group_modify(add_introns) %>% ungroup() %>% arrange(Start)
u <- u %>% group_by(Gene) %>% mutate(Transcript_nr = as.integer(factor(Transcript_ID, unique(Transcript_ID))))
u$Transcript_nr <- fct_rev(factor(u$Transcript_nr))
u <- u %>% group_by(Gene) %>% mutate(Range_min = min(Start), Range_max = max(Stop)) %>% ungroup()
u <- u %>% arrange(Range_min) %>% print(n = Inf)

nr_transcripts <- u %>% distinct(Transcript_ID, .keep_all = TRUE) %>% dplyr::count(Gene)
nr_transcripts <- set_names(nr_transcripts$n, nr_transcripts$Gene)
groups <- vector(mode = "list", length(nr_transcripts))
names(groups) <- names(nr_transcripts)
gene <- u$Gene[1]
groups[[gene]] <- seq(1, nr_transcripts[gene])

groups_stop <- list()
groups_stop[seq(1, nr_transcripts[gene])] <- u$Range_max[u$Gene == gene][1]
gene_last <- gene

# thresh <- 200000
thresh <- 10000

for (gene in unique(u$Gene)[-1]) {

  start <- u$Range_min[u$Gene == gene][1]
  stop <- u$Range_max[u$Gene == gene][1]

  which_group <- start > (unlist(groups_stop) + thresh)
  rle_which_group <- rle(which_group)
  put_after <- match(TRUE, rle_which_group$values & (rle_which_group$lengths >= nr_transcripts[gene]))

  if (is.na(put_after)) {

    last_true <- tail(rle_which_group$values, 1)

    if ((length(rle_which_group$length) - last_true) != 0) {

      put_in_these_groups <- cumsum(rle_which_group$lengths)[length(rle_which_group$length) - last_true] + seq(1, nr_transcripts[gene])

    } else {

      put_in_these_groups <- seq(1, nr_transcripts[gene])

    }


  } else {
    put_in_these_groups <- c(0, cumsum(rle_which_group$lengths))[put_after] + seq(1, nr_transcripts[gene])
  }

  groups[[gene]] <- put_in_these_groups
  groups_stop[put_in_these_groups] <- stop
  gene_last <- gene

}

u$ID_group <- "1"

for (gene in unique(u$Gene)) {

  transcript_ID <- u$Transcript_ID[u$Gene == gene]
  transcript_ID <- factor(transcript_ID, unique(transcript_ID))
  levels(transcript_ID) <- groups[[gene]]
  transcript_ID <- as.character(transcript_ID)
  u[u$Gene == gene, "ID_group"] <- transcript_ID

}

# u$Name <- paste0(u$Gene, "_", u$Transcript_nr)
# u$Name <- fct_rev(factor(u$Name, unique(u$Name)))
# u$Name <- fct_rev(factor(u$Name, levels = unique(u$Name)))

u$ID_group <- fct_rev(factor(u$ID_group, as.character(seq(1, max(as.integer(u$ID_group))))))

text_df <- u %>%
  group_by(Transcript_ID) %>%
  filter(Start == min(Start)) %>%
  ungroup()

p1 <- ggplot(u) +
  geom_segment(aes(x = Start, xend = Stop, yend = ID_group, y = ID_group), ~ filter(., Feature == "exon"), size = 5) +
  geom_segment(aes(x = Start, xend = Stop, yend = ID_group, y = ID_group), ~ filter(., Feature == "intron"), size = 0.5) +
  geom_text(aes(x = Start, y = ID_group, label = Gene), text_df, size = 3, hjust = "right", nudge_x = -1000) +
  xlim(c(min(u$Start) - thresh, max(u$Stop) + thresh)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

p1

ggsave("~/gene_track.pdf", p1, width = 120, height = 30, units = "cm")

state_bed_sub <- state_bed %>%
  filter(Chromosome == "chr22", Start > (min(u$Start) - thresh), Stop < (max(u$Stop) + thresh))

# p2 <- ggplot(state_bed_sub, aes(x = Start, xend = Stop, y = 1, yend = 1, color = State)) +
#   geom_segment() +
#   xlim(c(min(u$Start) - thresh, max(u$Stop) + thresh)) +
#   theme_classic()

age_ewas_gene <- age_ewas %>% 
  filter(grepl("DNMT3A", Probe_gene)) %>% 
  filter(P_FDR < 0.05)

p2 <- ggplot(age_ewas_gene, aes(x = Probe_position, y = Estimate)) +
  geom_point() + 
  xlim(c(min(u$Start) - thresh, max(u$Stop) + thresh)) +
  theme_classic() +
  ylab(NULL) +
  xlab(NULL)

p1 / p2


cells <- "Mononuclear_cells"
paths <- glue('./Data/Roadmap_15_states/{cells}_15_coreMarks_mnemonics.bed.gz')
names(paths) <- cells
state_bed <- vroom(paths,
                   col_names = c("Chromosome", "Start", "Stop", "State"),
                   id = "Cell") %>%
  mutate(Cell = factor(Cell, paths)) %>%
  mutate(Cell = fct_recode(Cell, !!!paths)) %>%
  filter(!Chromosome %in% c("chrM", "chrX", "chrY"))
