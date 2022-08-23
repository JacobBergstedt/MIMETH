#!/bin/bash
#SBATCH --mem=120G
#SBATCH --qos=geh
#SBATCH --partition=geh
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --verbose
#SBATCH --array=1-12


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
source /local/gensoft2/adm/etc/profile.d/modules.sh

# Load R
module load R/3.6.2

# Get Array ID
i=${SLURM_ARRAY_TASK_ID}

# Get i:th covariate
cov=$(cat "./Data/female_variables_list.txt" | awk "NR == "${i}"")


./Scripts/R_scripts/Run_scripts/EWAS/run_female_specific_ewas.R --cov="${cov}" \
                                                                --n_cores="${SLURM_CPUS_PER_TASK}" \
                                                                --null_formula="${null_formula}" \
                                                                --chunk_path="${chunk_path}"
