#!/bin/bash

# Define variables

# Chunk paths
get_m_value_env_const_chunk_path () {
  local path="./Data/Chunk_data/Results/EWAS/M_values/Environment/Constant_variance/"
  echo "${path}"
}

get_m_value_env_sandwich_chunk_path () {
  local path="./Data/Chunk_data/Results/EWAS/M_values/Environment/Sandwich/"
  echo "${path}"
}

get_m_value_female_specific_sandwich_chunk_path () {
  local path="./Data/Chunk_data/Results/EWAS/M_values/Female_specific/Sandwich/"
  echo "${path}"
}

get_m_value_female_specific_random_effect_chunk_path () {
  local path="./Data/Chunk_data/Results/EWAS/M_values/Female_specific/Random_effect/"
  echo "${path}"
}

get_m_value_env_random_effect_chunk_path () {
  local path="./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect/"
  echo "${path}"
}

get_m_value_env_random_effect_all_genotypes_chunk_path () {
  local path="./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes/"
  echo "${path}"
}

get_m_value_env_random_effect_normalized_props_884_chunk_path () {
  local path="./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_all_genotypes_normalized_props_884/"
  echo "${path}"
}

get_m_value_env_random_effect_no_genotypes_chunk_path () {
  local path="./Data/Chunk_data/Results/EWAS/M_values/Environment/Random_effect_no_genotypes/"
  echo "${path}"
}

# get_m_value_immune_const_chunk_path () {
#  local path="./Data/Chunk_data/Results/EWAS/M_values/Constant_variance/Immune_cells/"
#  echo "${path}"
# }
#
# get_m_value_immune_sandwich_chunk_path () {
#  local path="./Data/Chunk_data/Results/EWAS/M_values/Sandwich/Immune_cells/"
#  echo "${path}"
# }

####

get_beta_value_env_sandwich_chunk_path () {
  local path="./Data/Chunk_data/Results/EWAS/Beta_values/Environment/Sandwich/"
  echo "${path}"
}

get_beta_value_env_random_effect_chunk_path () {
  local path="./Data/Chunk_data/Results/EWAS/Beta_values/Environment/Random_effect/"
  echo "${path}"
}

# get_beta_value_props_sandwich_chunk_path () {
#   local path="./Data/Chunk_data/Results/EWAS/Beta_values/Sandwich/Props/"
#   echo "${path}"
# }
#
# get_beta_value_counts_sandwich_chunk_path () {
#   local path="./Data/Chunk_data/Results/EWAS/Beta_values/Sandwich/Counts/"
#   echo "${path}"
# }

# Result paths

get_m_value_env_const_result_path () {
  local path="./Data/RData/Results/EWAS/M_values/Environment/Constant_variance/"
  echo "${path}"
}

get_m_value_env_random_effect_result_path () {
  local path="./Data/RData/Results/EWAS/M_values/Environment/Random_effect/"
  echo "${path}"
}

get_m_value_env_sandwich_result_path () {
  local path="./Data/RData/Results/EWAS/M_values/Environment/Sandwich/"
  echo "${path}"
}

get_m_value_female_specific_sandwich_result_path () {
  local path="./Data/RData/Results/EWAS/M_values/Female_specific/Sandwich/"
  echo "${path}"
}


get_beta_value_env_result_path () {
  local path="./Data/RData/Results/EWAS/Beta_values/Environment/Sandwich/"
  echo "${path}"
}

get_beta_value_env_random_effect_result_path () {
  local path="./Data/RData/Results/EWAS/Beta_values/Environment/Random_effect/"
  echo "${path}"
}

get_m_value_immune_const_result_path () {
  local path="./Data/RData/Results/EWAS/M_values/Constant_variance/Immune_cells/"
  echo "${path}"
}

get_m_value_immune_sandwich_result_path () {
  local path="./Data/RData/Results/EWAS/M_values/Sandwich/Immune_cells/"
  echo "${path}"
}

# Cluster paths

get_m_value_env_const_cluster_path () {
  local path="./Cluster/EWAS/M_values/Environment/Constant_variance/"
  echo "${path}"
}

get_m_value_env_random_effect_cluster_path () {
  local path="./Cluster/EWAS/M_values/Environment/Random_effect/"
  echo "${path}"
}

get_m_value_env_random_effect_all_genotypes_cluster_path () {
  local path="./Cluster/EWAS/M_values/Environment/Random_effect_all_genotypes/"
  echo "${path}"
}

get_m_value_env_random_effect_normalized_props_884_cluster_path () {
  local path="./Cluster/EWAS/M_values/Environment/Random_effect_normalized_props_884/"
  echo "${path}"
}

get_m_value_env_random_effect_no_genotypes_cluster_path () {
  local path="./Cluster/EWAS/M_values/Environment/Random_effect_no_genotypes/"
  echo "${path}"
}

get_m_value_env_sandwich_cluster_path () {
  local path="./Cluster/EWAS/M_values/Environment/Sandwich/"
  echo "${path}"
}

get_m_value_female_specific_sandwich_cluster_path () {
  local path="./Cluster/EWAS/M_values/Female_specific/Sandwich/"
  echo "${path}"
}

get_m_value_female_specific_random_effect_cluster_path () {
  local path="./Cluster/EWAS/M_values/Female_specific/Random_effect/"
  echo "${path}"
}

get_beta_value_env_cluster_path () {
  local path="./Cluster/EWAS/Beta_values/Environment/Sandwich/"
  echo "${path}"
}

get_beta_value_env_random_effect_cluster_path () {
  local path="./Cluster/EWAS/Beta_values/Environment/Random_effect/"
  echo "${path}"
}


get_facs_gwas_path () {
  local path="../../../GWAS/LMM/FACS_GWAS/"
  echo "${path}"
}

clean_log_files () {
  find ./Cluster -type f -delete
}
