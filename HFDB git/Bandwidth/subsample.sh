#!/bin/bash
#SBATCH -A mth250064
#SBATCH -p wholenode
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=0-18:00:00
#SBATCH -o logs/subsample_bwscan_%j.out
#SBATCH -e logs/subsample_bwscan_%j.err

source ~/.bashrc
conda activate renv

mkdir -p logs

CACHE_DIR="$HOME/specmeans_cache_linear_filter_bwscan"
SCRIPT="specMeanCIhpc_linearFilter_bwScan.R"

echo "Starting fixed-block subsampling for 50x50 non-Gaussian linear filter bandwidth scan"
echo "Fixed block size: 10x10"
echo "Lag: (1,1)"
echo "Bandwidth grid: 0.05, 0.10, 0.15, 0.20, 0.25"
echo "Using cache dir: $CACHE_DIR"
date

Rscript ${SCRIPT} subsample "$CACHE_DIR"

date
echo "Finished fixed-block subsampling and bandwidth summary"

