#!/bin/bash

# Bash functions that stores all the paths
source ./Scripts/Bash_scripts/bash_logistics.sh

# Construct the formula
null_formula="$(cat ./Data/methylation_control_formula)"
echo "${null_formula}"

# Path to save results for each covariate
chunk_path="$(get_m_value_env_random_effect_normalized_props_884_chunk_path)""Correct_for_no_cells/"

jobid=$(sbatch \
        --error="$(get_m_value_env_random_effect_normalized_props_884_cluster_path)""Errors/correct_for_no_cells_%A_%a_errors" \
        --output="$(get_m_value_env_random_effect_normalized_props_884_cluster_path)""Outputs/correct_for_no_cells_%A_%a_outputs" \
        --parsable \
        ./Scripts/Bash_scripts/EWAS/array_ewas_random_effect_normalized_props_884.sh \
        --chunk_path="${chunk_path}" \
        --null_formula="${null_formula}")

echo "${jobid}"
