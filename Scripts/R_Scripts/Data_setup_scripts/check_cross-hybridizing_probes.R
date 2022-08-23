## Checking if sex-associated CpG sites are part of cross-hybridizing probes filtered in Pidsley et al. Genome Biol 2016 (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1#Sec22)
## or cross-hybridizing probes filtered in Price et al., Epigenetics Chromatin. 2013 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3740789/ 

# Read list from Price et al., Epigenetics Chromatin. 2013
chp_450 <- read.csv("Data/Quality_filters/Crosshybridizing_450K_Price.csv")
chp_450 <- chp_450 %>%
			select(IlmnID, XY_Hits, Autosomal_Hits)  %>%
			rename(Probe = IlmnID)
chp_450$Probe <- as.character(chp_450$Probe)

# Read list from Pidsley et al. Genome Biol 2016 
chp_EPIC <- read.table("Data/Quality_filters/Crosshybridizing_EPIC_Pidsley.csv", header = FALSE, col.names = c("Probe")) %>%
			add_column(EPIC_Hits = "EPIC_YES")
chp_EPIC$Probe <- as.character(chp_EPIC$Probe)

# Annotate ewas of sex
ewas_sex <- readRDS("./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect/Correct_for_16_props/Sex.rds")
ewas_sex <- ewas_sex %>%
			select(Probe, Estimate, Standard_error, P_FDR) %>%
			left_join(chp_450) %>%
			left_join(chp_EPIC) %>%
			replace_na(list(XY_Hits = "XY_NO", Autosomal_Hits = "A_NO", EPIC_Hits = "EPIC_NO"))


ewas_sex %>% filter(XY_Hits == "XY_YES" | Autosomal_Hits == "A_YES" | EPIC_Hits == "EPIC_YES") %>%
			arrange(P_FDR)


# Export list of probes to exclude
bad_probes <- ewas_sex %>%
				filter(XY_Hits == "XY_YES" | Autosomal_Hits == "A_YES" | EPIC_Hits == "EPIC_YES") %>%
				select(Probe)
saveRDS(bad_probes, file = "Data/Quality_filters/MIMETH.1028_missing_crosshybridizing_probes.rds")