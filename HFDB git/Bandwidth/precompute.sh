#!/bin/bash
#SBATCH -A mth250064
#SBATCH -p wholenode
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=1-12:00:00
#SBATCH -o logs/precompute_bwscan_%j.out
#SBATCH -e logs/precompute_bwscan_%j.err

source ~/.bashrc
conda activate renv

mkdir -p logs

CACHE_DIR="$HOME/specmeans_cache_linear_filter_bwscan"
SCRIPT="specMeanCIhpc_linearFilter_bwScan.R"

echo "Starting precompute for 50x50 non-Gaussian linear filter bandwidth scan"
echo "Fixed block size in downstream subsampling: 10x10"
echo "Bandwidth grid: 0.05, 0.10, 0.15, 0.20, 0.25"
echo "Using cache dir: $CACHE_DIR"
date

Rscript ${SCRIPT} precompute "$CACHE_DIR"

date
echo "Finished precompute for 50x50 non-Gaussian linear filter bandwidth scan"

