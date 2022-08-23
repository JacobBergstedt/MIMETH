#!/opt/gensoft/exe/R/3.6.2/scripts/Rscript

library(optparse)
library(parallel)
library(tidyverse)
library(MatrixEQTL)
library(splines)
source("./Scripts/R_scripts/Libraries/functions_for_inference.R")
source("./Scripts/R_scripts/Libraries/functions_for_logistics.R")
source("./Scripts/R_scripts/Libraries/functions_for_meQTL_mapping2.R")

# Parse options
option_list = list(
  make_option("--result_path", action = "store", type = "character"),
  make_option("--methylation_path", action = "store", type = "character"),
  make_option("--chr", action = "store", type = "character")
)

# External options
opt <- parse_args(OptionParser(option_list=option_list))
result_path <- opt$result_path
chr <- opt$chr
methylation_path <- opt$methylation_path

# Load data ---------------------------------------------------------------
cell_list_15_cells <- get_15_props()
meth_locations <- get_anno_meth_location()
probes_trans <- readRDS("./Data/RData/Methylation/probes_trans_50K.rds")

mf <- as.data.frame(get_covs_884())
rownames(mf) <- mf$SUBJID

meth <- readRDS(methylation_path)
meth_id <- meth$SUBJID
meth <- meth[probes_trans]
meth <- as.matrix(meth)
rownames(meth) <- meth_id
meth <- meth[rownames(mf), ]

fm_control <- get_control_fm_15_cells()
mf <- model.matrix(fm_control, mf)[, -1]

snp_locations <- readRDS(glue::glue("./Data/Chunk_data/Genotype/Chromosome_{chr}/LabExMI_imputation_1000x5699237_chr{chr}_map.rds"))
snps <- get_snp_matrix(chr)[rownames(mf), ]

meqtl_meth <- setup_methylation_for_meQTL_mapping(meth)
meqtl_snps <- setup_genotypes_for_meQTL_mapping(snps)
meqtl_covariates <- setup_covariates_for_meQTL_mapping(mf)

# Call the main analysis function
res <- Matrix_eQTL_main(snps = meqtl_snps,
                        gene = meqtl_meth,
                        cvrt = meqtl_covariates,
                        snpspos = as.data.frame(snp_locations[c("SNP_ID", "SNP_chr", "SNP_position")]),
                        genepos = as.data.frame(meth_locations[c("Probe", "Probe_chr", "Probe_position", "Probe_position")]),
                        pvOutputThreshold = 1e-3,
                        pvOutputThreshold.cis = 1e-3,
                        output_file_name = NULL,
                        output_file_name.cis = NULL,
                        cisDist = 5 * 1e5,
                        useModel = modelLINEAR,
                        verbose = TRUE)$trans$eqtls

res <- res %>% 
  clean_matrix_eQTL_results(meth_locations, snp_locations = snp_locations)

saveRDS(res, result_path)