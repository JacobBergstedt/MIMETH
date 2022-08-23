#!/bin/bash
#SBATCH --mem=120G
#SBATCH --cpus-per-task=6
#SBATCH --partition=geh
#SBATCH --qos=geh
#SBATCH --nodes=1
#SBATCH --verbose
#SBATCH --array=56,63,64,66,67,68,74,102,103,110

if (( $# != 2 ))
then
  echo "Not all parameters supplied"
  exit
fi

for arg in "$@"
do

  flag="${arg%%=*}"
  value="${arg#*=}"

  case "${flag}" in

    --chunk_path)
      chunk_path="${value}"
    ;;
    --null_formula)
      null_formula="${value}"
    ;;

    *)
      echo "Unknown parameter supplied"
      exit
    ;;

  esac

done

# Load loader
source /opt/gensoft/adm/etc/profile.d/modules.sh

# Load R
module load R/3.6.2

# Get Array ID
i=${SLURM_ARRAY_TASK_ID}

# Get i:th covariate
cov=$(cat "./Data/covariate_list.txt" | awk "NR == "${i}"")

# Print variables
echo "Covariate: "
echo "${cov}"
echo "Path to store intermediate results: "
echo "${chunk_path}"
echo "Statistical model: "
echo "${null_formula}"

./Scripts/R_scripts/Run_scripts/EWAS/run_ewas_random_effect_snps_884_samples.R --cov="${cov}" \
                                                                               --n_cores="${SLURM_CPUS_PER_TASK}" \
                                                                               --null_formula="${null_formula}" \
                                                                               --chunk_path="${chunk_path}"
