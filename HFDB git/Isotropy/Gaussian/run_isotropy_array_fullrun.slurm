#!/bin/bash
#SBATCH -A mth260002p
#SBATCH -p RM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=0-24:00:00
#SBATCH --array=0-5

source ~/.bashrc
conda activate renv

mkdir -p logs

SCRIPT="Isotropy_Test_Gaussian_HPC_array_fullrun.R"
BASE_OUTDIR="${HOME}/isotropy_gaussian_array_results"
psiR_values=(1 1.1 1.2 1.3 1.4 1.5)

psiR=${psiR_values[$SLURM_ARRAY_TASK_ID]}
OUTDIR="${BASE_OUTDIR}/psiR_${psiR}"

mkdir -p "$OUTDIR"

echo "Starting isotropy Gaussian full run"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "psiR: ${psiR}"
echo "Output directory: ${OUTDIR}"
date

Rscript "$SCRIPT" "$psiR" "$OUTDIR"

status=$?
date

echo "Finished isotropy Gaussian full run for psiR=${psiR}"
exit $status
