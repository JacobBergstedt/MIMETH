.libPaths("/pasteur/entites/Geh/Shared_Programs/Rlib/3.6.0/")
library(minfi)
library(matrixStats)
library(irlba)
library(sva)
library(lme4)
library(lmerTest)
library(tidyverse)


# Declare functions -------------------------------------------------------

#
# Replace by NA low detection P values and infinite values due to M ratio conversion
# Exclude probes with missingness >= 10/(sample size)
# Replace remaining NAs by probe mean
#
filter_impute_MMatrix <- function(MMatrix, detP) {
	# MMatrix[which(detP > 0.05, arr.ind=TRUE)] <- NA
	MMatrix[detP > 0.05] <- NA
	MMatrix[is.infinite(MMatrix)] <- NA
	probes.miss <- apply(MMatrix, 1, function(x) sum(is.na(x)))
	probes2exclude <- names(probes.miss)[probes.miss >= 10]
	MMatrix <- MMatrix[!rownames(MMatrix) %in% probes2exclude, ]
	NA_pos <- which(is.na(MMatrix), arr.ind = TRUE)
	MMatrix[NA_pos] <- rowMeans(MMatrix, na.rm = T)[NA_pos[, "row"]]
	MMatrix
}

#
# Perform PCA and output 50 PC eigenvalues
#
perform_pca_of_MMatrix <- function(MMatrix) {
	MMatrix <- t(MMatrix)
	PCA <- prcomp_irlba(MMatrix, n = 50, center = TRUE, scale. = TRUE)
	PCA <- as_tibble(PCA$x)
	PCA$Sample_Name <- mimeth.sheet$Sample_Name[mimeth.sheet$SentrixID %in% rownames(MMatrix)]
	PCA
}

#
# Convert a vector of M values to a vector of beta values
#
M2beta <- function(M) {
	beta <- 2^M / (2^M + 1)
	beta
}

#
# Replace by NAs outlier values of vector x
# Outlier values are defined as distant of k * sd from the mean
#
fix_outliers <- function(x, k) {
	M <- mean(x, na.rm = TRUE)
	SD <- sd(x, na.rm = TRUE)
	x[x > M + k * SD] <- NA
	x[x < M - k * SD] <- NA
	x
}


# ### fit_mixed_model with Sentrix_position as random effect
# fit_mixed_model <- function(y, db) {
#   db$y <- y
#   m <- lmer(y ~ X_CD3pos_of_total.panel1 + Age + Sex + CMVPositiveSerology + Smoking + MCHC + (1 | Sentrix_position) +
#               (1 | DateOfSampling), data = db, REML = FALSE)
# 
#   var_comps <- as.data.frame(VarCorr(m))
#   normalize <- sqrt(sum(var_comps$vcov))
# 
#   effect_sentrix <- var_comps %>%
#     dplyr::filter(grp == "Sentrix_position") %>%
#     pull(sdcor) / normalize
# 
#   effect_day <- var_comps %>%
#     dplyr::filter(grp == "DateOfSampling") %>%
#     pull(sdcor) / normalize
# 
#   p_sentrix <- anova(m, update(m, . ~ . - (1 | Sentrix_position)))[["Pr(>Chisq)"]][2]
#   p_day <- anova(m, update(m, . ~ . - (1 | DateOfSampling)))[["Pr(>Chisq)"]][2]
# 
#   summary(m)$coefficients %>%
#     as.data.frame() %>%
#     rownames_to_column(var = "Term") %>%
#     dplyr::filter(Term != "(Intercept)") %>%
#     dplyr::select(Term, Effect = Estimate, P_value = `Pr(>|t|)`) %>%
#     rbind(tibble(Term = c("Sentrix_position", "DateOfSampling"),
#                  Effect = c(effect_sentrix, effect_day),
#                  P_value = c(p_sentrix, p_day))) %>%
#     as_tibble()
# }

### fit_mixed_model with batch effects as random effect
fit_mixed_model <- function(y, db, test_var) {
  db$y <- y
  test_var_term <- paste0("(1 | ", test_var, ")")
  fm <- y ~ X_CD3pos_of_total.panel1 + Age + Sex + CMVPositiveSerology + Smoking + MCHC + (1 | DateOfSampling)
  fm <- update(fm, paste0(". ~ . + " , test_var_term))
  m <- lmer(fm, data = db, REML = FALSE)
  
  var_comps <- as.data.frame(VarCorr(m))
  normalize <- sqrt(sum(var_comps$vcov))

  effect_test <- var_comps %>%
    dplyr::filter(grp == test_var) %>%
    pull(sdcor) / normalize
  
  effect_day <- var_comps %>%
    dplyr::filter(grp == "DateOfSampling") %>%
    pull(sdcor) / normalize

  p_test <- anova(m, update(m, paste0(". ~ . - ", test_var_term)))[["Pr(>Chisq)"]][2]
  p_day <- anova(m, update(m, . ~ . - (1 | DateOfSampling)))[["Pr(>Chisq)"]][2]

  summary(m)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column(var = "Term") %>%
    dplyr::filter(Term != "(Intercept)") %>%
    dplyr::select(Term, Effect = Estimate, P_value = `Pr(>|t|)`) %>%
    rbind(tibble(Term = c(test_var, "DateOfSampling"),
                 Effect = c(effect_test, effect_day),
                 P_value = c(p_test, p_day))) %>%
    as_tibble()
}


#
# Y: dataframe of numerical columns that will be response variables. In particular shouldn't contain ID's
# Y_IDs: identifiers for the rows of Y
# db: database containing the predictors X_CD3pos_of_total.panel1, Age, Sex, CMVPositiveSerology, Smoking, Sentrix_position, DateOfSampling
#
#
col_fit <- function(Y, Y_IDs, db, test_var) {
  db <- right_join(db, tibble(SUBJID = Y_IDs))
  map_dfr(Y, fit_mixed_model, db = db, test_var = test_var, .id = "Response")
}

get_subjid <- function(sample_name) {
  gsub("_.*", "", sample_name)
}

strongest_signals <- function(tib) {
  tib %>% 
    mutate(Response = as.numeric(str_sub(Response, start = 3, end = 4))) %>% 
    group_by(Response) %>% 
    summarize(strongest_signal = Term[which.min(P_value)], P_value = min(P_value)) %>% 
    print(n = Inf)
}
 


# Create RGSet object -----------------------------------------------------
# Read Illumina .csv file and import names of all .idat files
mimeth.sheet <- read.metharray.sheet("/pasteur/projets/policy01/LabExMI/Data/Methylation/2016_Milieu/")
mimeth.sheet$SentrixID <- paste0(mimeth.sheet$Slide, "_", mimeth.sheet$Array)
mimeth.sheet$SUBJID <- get_subjid(mimeth.sheet$Sample_Name)
mimeth.sheet$Sentrix_position <- str_sub(mimeth.sheet$Array, start = 1, end = 3)


# Read raw IDAT files
mimeth.RGset <- read.metharray.exp(targets = mimeth.sheet)
saveRDS(mimeth.RGset, file = "./Data/RData/Methylation/MIMETH.minfi.RGset.rds")
#mimeth.RGset <- readRDS("Data/RData/Methylation/MIMETH.minfi.RGset.rds")



# Raw: Control quality filters of samples ---------------------------------
# Sample quality based on detection P-values
mimeth.detP <- detectionP(mimeth.RGset)
mimeth.sheet$mean.detP <- colMeans(mimeth.detP)
samples2exclude <- as.vector(mimeth.sheet$Sample_Name[mimeth.sheet$mean.detP > 0.005]) #only '_BAD' samples

# Generate raw MethylSet for QC filters
mimeth.MSet.raw <- preprocessRaw(mimeth.RGset)

# Sample quality based on medians of Meth and Unmeth channels
mimeth.qcMed <- getQC(mimeth.MSet.raw)
mimeth.sheet$mMed <- mimeth.qcMed$mMed
mimeth.sheet$uMed <- mimeth.qcMed$uMed
samples2exclude <- c(samples2exclude, as.vector(mimeth.sheet$Sample_Name[mimeth.sheet$mMed < mean(mimeth.sheet$mMed)-3*sd(mimeth.sheet$mMed)])) #only '_BAD' samples

# Predict sex based on X/Y methylation
mimeth.GMSet.raw <- mapToGenome(mimeth.MSet.raw)
mimeth.getSex <- getSex(mimeth.GMSet.raw)
mimeth.getSex <- rownames_to_column(as.data.frame(mimeth.getSex))
colnames(mimeth.getSex)[1] <- "SentrixID"
mimeth.sheet <- dplyr::left_join(mimeth.sheet, mimeth.getSex, by = "SentrixID")
mimeth.sheet$predictedSex <- recode(mimeth.sheet$predictedSex, M = "Male", F = "Female")

# Compare predicted sex to actual sex
ecrf <- readRDS("./Data/RData/Environment/ecrf_curated.rds")
ecrf <- select(ecrf$V1, "SUBJID", "Age", "Sex")
mimeth.sheet <- dplyr::left_join(mimeth.sheet, ecrf, by = "SUBJID")
samples2exclude <- c(samples2exclude, mimeth.sheet$Sample_Name[mimeth.sheet$predictedSex != mimeth.sheet$Sex]) #only 1 sample: 962_V2

# Get EPIC genotypes for the 59 control SNPs
mimeth.controlSNPs <- getSnpBeta(mimeth.RGset)
write.table(rownames(mimeth.controlSNPs), file = "./Data/Quality_filters/MIMETH.59_controlSNPs.txt", quote = F, row.names = F, col.names = F)
mimeth.controlSNPs <- t(mimeth.controlSNPs)
mimeth.controlSNPs <- apply(mimeth.controlSNPs, 2, function(x) x * 2)
mimeth.controlSNPs <- rownames_to_column(as.data.frame(mimeth.controlSNPs))
colnames(mimeth.controlSNPs)[1] <- "SentrixID"
mimeth.controlSNPs <- dplyr::left_join(mimeth.controlSNPs, mimeth.sheet[, c("SentrixID", "SUBJID")], by = "SentrixID")

# Get Omni genotypes for the 59 control SNPs
omni.controlSNPs <- read.table("./Data/LabExMI_29_controlSNPs.raw", header = T, sep = " ")
colnames(omni.controlSNPs)[2] <- "SUBJID"
omni.controlSNPs$SUBJID <- as.character(omni.controlSNPs$SUBJID)
omni.controlSNPs <- omni.controlSNPs[, -c(1,3:6)]
SNPlist <- gsub("_.*", "", colnames(omni.controlSNPs)[grepl("rs", colnames(omni.controlSNPs))])

# Merge EPIC and Omni genotypes for 29 common SNPs
mimeth.controlSNPs <- mimeth.controlSNPs[, c("SentrixID", "SUBJID", SNPlist)]
mimeth.controlSNPs <- dplyr::left_join(mimeth.controlSNPs, omni.controlSNPs, by = "SUBJID")

# Match allele codes between datasets
SNPlist.corr <- data.frame(corr=rep(0, length(SNPlist)), row.names = SNPlist)
for (i in 3:31) {
	SNPlist.corr[i-2,] <- cor(mimeth.controlSNPs[,i], mimeth.controlSNPs[,i+29], use = "complete.obs")
	if (SNPlist.corr[i-2,] < 0) {
		mimeth.controlSNPs[,i+29] <- abs(mimeth.controlSNPs[,i+29] - 2)
	}
	SNPlist.corr[i-2,] <- cor(mimeth.controlSNPs[,i], mimeth.controlSNPs[,i+29], use = "complete.obs")
}
#SNPlist.corr #All correlations are > 99%

# Calculate mean discordance per sample
mimeth.controlSNPs$discord <- 0
for (i in 1:length(mimeth.controlSNPs$SUBJID)) {
	for (j in 3:31) {
		if (!is.na(mimeth.controlSNPs[i,j+29])) {
			mimeth.controlSNPs$discord[i] <-  mimeth.controlSNPs$discord[i] + abs(mimeth.controlSNPs[i,j] - mimeth.controlSNPs[i,j+29])
		}
	}
	mimeth.controlSNPs$discord[i] <- mimeth.controlSNPs$discord[i] / sum(!is.na(mimeth.controlSNPs[i,3:31]))
}
mimeth.sheet <- dplyr::left_join(mimeth.sheet, mimeth.controlSNPs[, c("SentrixID", "discord")], by = "SentrixID")
samples2exclude <- c(samples2exclude, as.vector(mimeth.sheet$Sample_Name[mimeth.sheet$discord > mean(mimeth.sheet$discord)+3*sd(mimeth.sheet$discord)])) #2 samples: 529_3260590139, 962_V2

#Add _BAD, PBMC, V2 samples to samples2exclude
samples2exclude <- c(samples2exclude, as.vector(mimeth.sheet$Sample_Name[grepl("BAD|V2|PBMC|v2", mimeth.sheet$Sample_Name)]))

# Add replicates to samples2exclude (chosen manually based on higher mean.detP, higher discord, and external raw.PCA positions)
samples2exclude <- c(samples2exclude, "11_3260590265_Rep2", "310_3260590430", "488_3260587777_Rep1", "513_3260588503", "60_3260589989_Rep1", "811_3260589730_Rep1", "893_3260587692_Rep2", "893_3260587692_Rep3", "961_3260586195")
samples2exclude <- samples2exclude[!duplicated(samples2exclude)]

# Exclude samples listed in samples2exclude: 16 BAD, 27 V2, 57 PBMC, 9 replicates, 1 discordant
mimeth.RGset.filt <- mimeth.RGset[, ! mimeth.RGset$Sample_Name %in% samples2exclude]

# Generate raw MSet, detP and MMatrix after initial sample filter
mimeth.MSet.raw <- preprocessRaw(mimeth.RGset.filt)
mimeth.MMatrix.raw <- getM(mimeth.MSet.raw)
mimeth.detP <- detectionP(mimeth.RGset.filt)
saveRDS(mimeth.detP, file = "./Data/RData/Methylation/MIMETH.minfi.detP_978.rds")

# Remove probes based on #detP>0.05, #NA or #inf >10 and impute remaining NAs by rowMeans
# Excluded 3311 probes listed in probes2exclude: 3250 with #detP>=10 (including #missing >=10), 61 with #inf>10
mimeth.MMatrix.raw.filt <- filter_impute_MMatrix(mimeth.MMatrix.raw, mimeth.detP)
saveRDS(mimeth.MMatrix.raw.filt, file = "./Data/RData/Methylation/MIMETH.minfi.MMatrix.raw_978.rds")

# Perform PCA of filtered raw M matrix
mimeth.raw.PCA <- perform_pca_of_MMatrix(mimeth.MMatrix.raw.filt)
saveRDS(mimeth.raw.PCA, file = "./Data/RData/Methylation/MIMETH.PCA.raw.rds")

# Write all QC statistics to RDS file
saveRDS(mimeth.sheet, file = "./Data/RData/Methylation/Annotation/MIMETH.sample_sheet.rds")
#mimeth.sheet <- readRDS("Data/RData/Methylation/Annotation/MIMETH.sample_sheet.rds")
qcReport(mimeth.RGset.filt, pdf= "./Plots/Methylation_processing/MIMETH.minfi.qcReport.pdf")



# Noob: normalize RGset ---------------------------------------------------
mimeth.MSet.noob <- preprocessNoob(mimeth.RGset.filt, dyeMethod = "reference")
saveRDS(mimeth.MSet.noob, file = "./Data/RData/Methylation/MIMETH.minfi.MSet.noob_978.rds")
mimeth.MMatrix.noob <- getM(mimeth.MSet.noob)
mimeth.MMatrix.noob.filt <- filter_impute_MMatrix(mimeth.MMatrix.noob, mimeth.detP)

# Perform PCA of noob M matrix
mimeth.noob.PCA <- perform_pca_of_MMatrix(mimeth.MMatrix.noob.filt)
saveRDS(mimeth.noob.PCA, file = "./Data/RData/Methylation/MIMETH.PCA.noob_978.rds")

# List PCA outliers
samples2exclude <- as.vector(mimeth.sheet$Sample_Name[mimeth.sheet$Slide == "200298840149"])
samples2exclude <- c(samples2exclude, as.vector(mimeth.noob.PCA$Sample_Name[mimeth.noob.PCA$PC1 > 1500]), 
                     as.vector(mimeth.noob.PCA$Sample_Name[mimeth.noob.PCA$PC2 > 1000])) # 309_3260590429
samples2exclude <- samples2exclude[!duplicated(samples2exclude)]
write.table(samples2exclude, file = "/pasteur/projets/policy01/LabExMI/Methylation/MIMETH/Data/MIMETH.9_excluded_samples.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Exclude samples in samples2exclude: 9 PCA outliers
mimeth.RGset.filt <- mimeth.RGset.filt[, !mimeth.RGset.filt$Sample_Name %in% samples2exclude]
saveRDS(mimeth.RGset.filt, file = "./Data/RData/Methylation/MIMETH.minfi.RGset_969.rds")
#mimeth.RGset.filt <- readRDS("Data/RData/Methylation/MIMETH.minfi.RGset_969.rds")

# Sheet with bad samples removed
mimeth.sheet.filtered <- mimeth.sheet %>% 
  filter(SentrixID %in% colnames(mimeth.RGset.filt))
saveRDS(mimeth.sheet.filtered, file = "./Data/RData/Methylation/Annotation/MIMETH.sample_sheet_969.rds")
#mimeth.sheet.filtered <- readRDS("Data/RData/Methylation/Annotation/MIMETH.sample_sheet_969.rds")



# Noob: normalize filtered RGset ------------------------------------------
mimeth.MSet.noob <- preprocessNoob(mimeth.RGset.filt, dyeMethod = "reference")
saveRDS(mimeth.MSet.noob, file = "./Data/RData/Methylation/MIMETH.minfi.MSet.noob_969.rds")

mimeth.detP <- detectionP(mimeth.RGset.filt)
saveRDS(mimeth.detP, file = "./Data/RData/Methylation/MIMETH.minfi.detP_969.rds")

mimeth.MMatrix.noob <- getM(mimeth.MSet.noob)
mimeth.MMatrix.noob.filt <- filter_impute_MMatrix(mimeth.MMatrix.noob, mimeth.detP)
# 286289 values are missing out of 866836*969=839964084 (0.034%)
# Excluded 2930 probes
# Imputed 31660 values out of 863906*969=837124914 (0.0038%)
saveRDS(mimeth.MMatrix.noob.filt, file = "./Data/RData/Methylation/MIMETH.minfi.MMatrix.noob_969.rds")



# Noob: Perform LMM on PCs ------------------------------------------------
mimeth.noob.PCA <- perform_pca_of_MMatrix(mimeth.MMatrix.noob.filt)
saveRDS(mimeth.noob.PCA, file = "./Data/RData/Methylation/MIMETH.PCA.noob_969.rds")

ecrf <- readRDS("./Data/RData/Environment/ecrf_curated.rds")$V1
cells <- readRDS("./Data/RData/Cells/proportions_of_total_outrm.rds") %>% 
  mutate(SUBJID = as.character(SUBJID)) %>% 
  select(SUBJID, X_CD3pos_of_total.panel1)

db <- full_join(ecrf, cells) %>% 
  inner_join(mimeth.sheet.filtered[c("SUBJID", "Sentrix_position", "Sample_Plate")])

mimeth.noob.PCA.mixed.fit <- col_fit(select(mimeth.noob.PCA, -Sample_Name), 
                                            get_subjid(mimeth.noob.PCA$Sample_Name), 
		              db = db, 
		              test_var = "Sentrix_position")
strongest_signals(mimeth.noob.PCA.mixed.fit)



# Noob: ComBat correction for Sentrix_position ----------------------------
Sentrix_position <- substr(mimeth.MSet.noob$Array, 0, 3)
mimeth.MMatrix.noob.ComBat <- ComBat(mimeth.MMatrix.noob.filt,
                                        batch = Sentrix_position, 
                                        mod = NULL, 
                                        par.prior = TRUE, 
                                        prior.plots = FALSE)
saveRDS(mimeth.MMatrix.noob.ComBat, file = "./Data/RData/Methylation/MIMETH.minfi.MMatrix.noob_969.ComBat.rds")

mimeth.noob.ComBat.PCA <- perform_pca_of_MMatrix(mimeth.MMatrix.noob.ComBat)
saveRDS(mimeth.noob.ComBat.PCA, file = "./Data/RData/Methylation/MIMETH.PCA.noob_969.ComBat.rds")
mimeth.noob.ComBat.PCA.mixed.fit <- col_fit(select(mimeth.noob.ComBat.PCA, -Sample_Name), 
                                            get_subjid(mimeth.noob.ComBat.PCA$Sample_Name), 
		              db = db, 
		              test_var = "Sentrix_position")
strongest_signals(mimeth.noob.ComBat.PCA.mixed.fit)



# Noob: ComBat correction for Sample_Plate ------------------------------
mimeth.MMatrix.noob.ComBat2 <- ComBat(mimeth.MMatrix.noob.ComBat,
                                        batch = mimeth.sheet.filtered$Sample_Plate, 
                                        mod = NULL, 
                                        par.prior = TRUE, 
                                        prior.plots = FALSE)
saveRDS(mimeth.MMatrix.noob.ComBat2, file = "./Data/RData/Methylation/MIMETH.minfi.MMatrix.noob_969.ComBat2.rds")

mimeth.noob.ComBat2.PCA <- perform_pca_of_MMatrix(mimeth.MMatrix.noob.ComBat2)
saveRDS(mimeth.noob.ComBat2.PCA, file = "./Data/RData/Methylation/MIMETH.PCA.noob_969.ComBat2.rds")
mimeth.noob.ComBat2.PCA.mixed.fit <- col_fit(select(mimeth.noob.ComBat2.PCA, -Sample_Name), 
                                            get_subjid(mimeth.noob.ComBat2.PCA$Sample_Name), 
		              db = db, 
		              test_var = "Sample_Plate")
strongest_signals(mimeth.noob.ComBat2.PCA.mixed.fit)



# Funnorm: normalize filtered RGset ---------------------------------------
mimeth.GRSet.funnorm <- preprocessFunnorm(mimeth.RGset.filt, nPCs = 2, sex = NULL, bgCorr = TRUE, dyeCorr = TRUE)
saveRDS(mimeth.GRSet.funnorm, file = "./Data/RData/Methylation/MIMETH.minfi.GRSet.funnorm_969.rds")

mimeth.MMatrix.funnorm <- getM(mimeth.GRSet.funnorm)
mimeth.detP.funnorm <- mimeth.detP[rownames(mimeth.MMatrix.funnorm), ]
mimeth.MMatrix.funnorm.filt <- filter_impute_MMatrix(mimeth.MMatrix.funnorm, mimeth.detP.funnorm)
saveRDS(mimeth.MMatrix.funnorm.filt, file = "./Data/RData/Methylation/MIMETH.minfi.MMatrix.funnorm_969.rds")



# Funnorm: Perform LMM on PCs ---------------------------------------------
mimeth.funnorm.PCA <- perform_pca_of_MMatrix(mimeth.MMatrix.funnorm.filt)
saveRDS(mimeth.funnorm.PCA, file = "./Data/RData/Methylation/MIMETH.PCA.funnorm_969.rds")

mimeth.funnorm.PCA.mixed.fit <- col_fit(select(mimeth.funnorm.PCA, -Sample_Name), 
                                            get_subjid(mimeth.funnorm.PCA$Sample_Name), 
		              db = db, 
		              test_var = "Sentrix_position")
strongest_signals(mimeth.funnorm.PCA.mixed.fit)



# Funnorm: ComBat correction for sentrix position -------------------------
Sentrix_position <- substr(mimeth.GRSet.funnorm$Array, 0, 3)
mimeth.MMatrix.funnorm.ComBat <- ComBat(mimeth.MMatrix.funnorm.filt, 
                                        batch = Sentrix_position, 
                                        mod = NULL, 
                                        par.prior = TRUE, 
                                        prior.plots = FALSE)
saveRDS(mimeth.MMatrix.funnorm.ComBat, file = "./Data/RData/Methylation/MIMETH.minfi.MMatrix.funnorm_969.ComBat.rds")

mimeth.funnorm.ComBat.PCA <- perform_pca_of_MMatrix(mimeth.MMatrix.funnorm.ComBat)
saveRDS(mimeth.funnorm.ComBat.PCA, file = "./Data/RData/Methylation/MIMETH.PCA.funnorm_969.ComBat.rds")

mimeth.funnorm.ComBat.PCA.mixed.fit <- col_fit(select(mimeth.funnorm.ComBat.PCA, -Sample_Name), 
                                            get_subjid(mimeth.funnorm.ComBat.PCA$Sample_Name), 
		              db = db, 
		              test_var = "Sample_Plate")
strongest_signals(mimeth.funnorm.ComBat.PCA.mixed.fit)



# Additional quality filters of probes -----------------------------------
# Check that manifest is correct
getManifest(mimeth.RGset.filt)

# List probes based on multiple matches to the genome (done by Lucas based on Maud's scripts)
multmatch <- scan(file = "./Data/Quality_filters/MIMETH.83635_multiple_match_probes.txt", what = character())

# List probes based on multiple matches to the genome (done by Lucas based on Maud's script: Lucas/Normalisation/Scripts/Scripts_Maud/find_snp_in_probes_1000G.pl)
overlapSNP <- scan(file = "/pasteur/projets/policy01/LabExMI/Methylation/Lucas/Normalisation/MMprobes/Probes_wCpG_1000G.txt", what = character())          # SNP in CpG
overlapSNP10 <- scan(file = "/pasteur/projets/policy01/LabExMI/Methylation/Lucas/Normalisation/MMprobes/Probes_overlap10_1000G.txt", what = character())   # SNP in [-1,-10] bp of CpG 
overlapSNP40 <- scan(file = "/pasteur/projets/policy01/LabExMI/Methylation/Lucas/Normalisation/MMprobes/Probes_overlap40_1000G.txt", what = character())   # SNP in [-11,-50] bp of CpG 
overlapSNP <- c(overlapSNP, overlapSNP10, overlapSNP40)
overlapSNP <- overlapSNP[!duplicated(overlapSNP)]

# Remove 969 probes absent from the newest version of the manifest
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
probes2exclude <- c(multmatch, overlapSNP, setdiff(rownames(mimeth.MMatrix.noob.ComBat2), rownames(Locations)))
probes2exclude <- probes2exclude[!duplicated(probes2exclude)]



# Generate the final M and beta matrices  ----------------------------------
mimeth.MMatrix.noob.ComBat2.filt <- mimeth.MMatrix.noob.ComBat2[!rownames(mimeth.MMatrix.noob.ComBat2) %in% probes2exclude, ]
final.MMatrix <- as_tibble(t(mimeth.MMatrix.noob.ComBat2.filt))
final.MMatrix <- add_column(final.MMatrix, mimeth.sheet.filtered$SUBJID, .before = 1)
colnames(final.MMatrix)[1] <- "SUBJID"
saveRDS(final.MMatrix, file = "./Data/RData/Methylation/MIMETH.minfi.final.MMatrix.rds")

final.betaMatrix <- bind_cols(final.MMatrix[1], map(final.MMatrix[-1], M2beta))
saveRDS(final.betaMatrix, file = "./Data/RData/Methylation/MIMETH.minfi.final.betaMatrix.rds")

# Out of 866836 probes, 661393 remain
# 2930 removed based on missingness
# 83380 removed based on multiple matches
# 118575 removed based on SNP overlap
# 558 removed based on location annotation



# Remove outlier values  ---------------------------------------------------
final.MMatrix.no_outliers <- bind_cols(final.MMatrix[1], map(final.MMatrix[-1], ~ fix_outliers(., 5)))
sum(as.numeric(map(final.MMatrix.no_outliers, ~ sum(is.na(.))))) / (dim(final.MMatrix.no_outliers)[1] *  dim(final.MMatrix.no_outliers)[2])  ## = 0.03435%
saveRDS(final.MMatrix.no_outliers, file = "./Data/RData/Methylation/MIMETH.minfi.final.MMatrix.no_outliers.rds")

final.betaMatrix.no_outliers <- bind_cols(final.betaMatrix[1], map(final.betaMatrix[-1], ~ fix_outliers(., 5)))
sum(as.numeric(map(final.betaMatrix.no_outliers, ~ sum(is.na(.))))) / (dim(final.betaMatrix.no_outliers)[1] *  dim(final.betaMatrix.no_outliers)[2]) ## = 0.09437%
saveRDS(final.betaMatrix.no_outliers, file = "./Data/RData/Methylation/MIMETH.minfi.final.betaMatrix.no_outliers.rds")


