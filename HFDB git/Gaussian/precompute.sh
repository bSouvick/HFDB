#!/bin/bash
#SBATCH -A mth260002p
#SBATCH -p RM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=0-12:00:00
#SBATCH -o output.%j
#SBATCH -e error.%j

source ~/.bashrc
conda activate renv

mkdir -p logs

CACHE_DIR="$HOME/specmeans_cache_g95_gausscov2"
mkdir -p "$CACHE_DIR"

SCRIPT="specMeanCIhpc_gaussian95_gausscov.R"

echo "Starting precompute for Gaussian field (95% CI setting)"
echo "Using cache dir: $CACHE_DIR"
date

Rscript ${SCRIPT} precompute 5 "$CACHE_DIR"

date
echo "Finished precompute for Gaussian field (95% CI setting)"

