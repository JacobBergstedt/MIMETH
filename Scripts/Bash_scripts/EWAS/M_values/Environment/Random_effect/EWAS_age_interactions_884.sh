#!/bin/bash
sbatch --output="./Cluster/EWAS/M_values/Environment/Random_effect/Age_interaction/Outputs/age_interaction_sex.output" \
       --error="./Cluster/EWAS/M_values/Environment/Random_effect/Age_interaction/Errors/age_interaction_sex.error" \
       ./Scripts/R_scripts/Run_scripts/EWAS/run_ewas_age_interaction_884.R \
       --interacting_var="Sex" \
       --n_cores=6

sbatch --output="./Cluster/EWAS/M_values/Environment/Random_effect/Age_interaction/Outputs/age_interaction_smoking.output" \
       --error="./Cluster/EWAS/M_values/Environment/Random_effect/Age_interaction/Errors/age_interaction_smoking.error" \
       ./Scripts/R_scripts/Run_scripts/EWAS/run_ewas_age_interaction_884.R \
       --interacting_var="Smoking_status" \
       --n_cores=6

sbatch --output="./Cluster/EWAS/M_values/Environment/Random_effect/Age_interaction/Outputs/age_interaction_cmv.output" \
       --error="./Cluster/EWAS/M_values/Environment/Random_effect/Age_interaction/Errors/age_interaction_cmv.error" \
       ./Scripts/R_scripts/Run_scripts/EWAS/run_ewas_age_interaction_884.R \
       --interacting_var="CMV_serostatus" \
       --n_cores=6
