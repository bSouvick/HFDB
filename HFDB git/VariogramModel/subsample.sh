#!/bin/bash
#SBATCH -p cpu
#SBATCH --account=bfss-delta-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=0-14:00:00
#SBATCH -o output.%j
#SBATCH -e error.%j
#SBATCH --array=1-6

source ~/.bashrc
conda activate renv

mkdir -p logs

CACHE_DIR="$HOME/variogram_model_ng_cache"
SCRIPT="specMeanCIhpc_variogramModel.R"

blocks=(5 6 7 8 9 10)
BETA=${blocks[$SLURM_ARRAY_TASK_ID-1]}

echo "Starting subsample for variogram model non-Gaussian case, block size $BETA"
echo "Using cache dir: $CACHE_DIR"
date

Rscript ${SCRIPT} subsample "$BETA" "$CACHE_DIR"

date
echo "Finished subsample for variogram model non-Gaussian case, block size $BETA"
