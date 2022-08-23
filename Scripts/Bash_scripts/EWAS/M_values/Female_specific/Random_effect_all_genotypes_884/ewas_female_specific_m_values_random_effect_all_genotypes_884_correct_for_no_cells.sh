#!/bin/bash

# Bash functions that stores all the paths
source ./Scripts/Bash_scripts/bash_logistics.sh

# Construct the formula
null_formula="$(cat ./Data/methylation_control_formula)"

echo "${null_formula}"

# Path to save results for each covariate
chunk_path="./Data/Chunk_data/Results/EWAS/M_values/Female_specific/Random_effect_all_genotypes_884/Correct_for_no_cells/"

jobid=$(sbatch \
        --error="./Cluster/EWAS/M_values/Female_specific/Random_effect_all_genotypes_884/Errors/correct_for_no_cells_%A_%a_errors" \
        --output="./Cluster/EWAS/M_values/Female_specific/Random_effect_all_genotypes_884/Outputs/correct_for_no_cells_%A_%a_outputs" \
        --parsable \
        ./Scripts/Bash_scripts/EWAS/array_ewas_female_specific.sh \
        --chunk_path="${chunk_path}" \
        --null_formula="${null_formula}")

echo "${jobid}"
