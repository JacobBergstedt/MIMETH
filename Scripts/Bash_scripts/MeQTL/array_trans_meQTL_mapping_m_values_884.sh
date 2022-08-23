#!/bin/bash
#SBATCH --partition=geh
#SBATCH --mem=30G
#SBATCH --output ./Cluster/MeQTL/Trans_meQTL/Outputs/jarray_%A_%a_outputs
#SBATCH --error ./Cluster/MeQTL/Trans_meQTL/Errors/jarray_%A_%a_errors
#SBATCH --nodes=1                   ### Node count required for the job
#SBATCH --verbose                   ### Increase informational messages
#SBATCH --array=1-22

# Load loader
source /opt/gensoft/adm/etc/profile.d/modules.sh

# Load R
module load R/3.6.2

# snp input data path
chr="${SLURM_ARRAY_TASK_ID}"
methylation_path="./Data/RData/Methylation/MIMETH.minfi.final.MMatrix.autosomes.no_outliers.rds"
result_path="./Data/Chunk_data/Results/MeQTL/Trans_meQTL/M_values/Normalized_cells_884/chromosome_""${chr}"".rds"

echo "${result_path}"

./Scripts/R_scripts/Run_scripts/MeQTL/run_trans_meQTL_mapping_884.R --result_path="${result_path}" \
                                                                    --methylation_path="${methylation_path}" \
                                                                    --chr="${chr}"
