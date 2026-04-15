#!/bin/bash
#SBATCH -p cpu
#SBATCH --account=bfss-delta-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=0-12:00:00
#SBATCH -o output.%j
#SBATCH -e error.%j

source ~/.bashrc
conda activate renv

mkdir -p logs

CACHE_DIR="$HOME/variogram_model_ng_cache"
mkdir -p "$CACHE_DIR"

SCRIPT="specMeanCIhpc_variogramModel.R"

echo "Starting precompute for variogram model non-Gaussian case"
echo "Using cache dir: $CACHE_DIR"
date

Rscript ${SCRIPT} precompute 7 "$CACHE_DIR"

date
echo "Finished precompute for variogram model non-Gaussian case"
