#!/bin/bash
#SBATCH -A mth260002p
#SBATCH -p RM
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

CACHE_DIR="$HOME/specmeans_cache_linear_filter"
SCRIPT="specMeanCIhpc_linearFilter.R"

blocks=(6 7 8 9 10 11)
BETA=${blocks[$SLURM_ARRAY_TASK_ID-1]}

echo "Starting subsample for non-Gaussian linear filter case, block size $BETA"
echo "Using cache dir: $CACHE_DIR"
date

Rscript ${SCRIPT} subsample "$BETA" "$CACHE_DIR"

date
echo "Finished subsample for non-Gaussian linear filter case, block size $BETA"
