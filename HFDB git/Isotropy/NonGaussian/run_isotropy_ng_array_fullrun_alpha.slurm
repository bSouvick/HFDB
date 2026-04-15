#!/bin/bash
#SBATCH -A mth250064
#SBATCH -p wholenode
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=0-24:00:00
#SBATCH -o logs/isotropy_ng_array_%A_%a.out
#SBATCH -e logs/isotropy_ng_array_%A_%a.err
#SBATCH --array=0-7

source ~/.bashrc
conda activate renv

mkdir -p logs

SCRIPT="IsotropyTest_ng_HPC_array_fullrun_spherical_quadratic.R"
BASE_OUTDIR="${HOME}/isotropy_nongaussian_array_results"
psiR_values=(1 1.2 1.4 1.6 1.8 2.0 2.2 2.3)

# Set the non-Gaussian strength here, or override at submission time:
# sbatch --export=ALL,ALPHA_NG=0.7 run_isotropy_ng_array_fullrun_alpha.slurm
ALPHA_NG=${ALPHA_NG:-0.15}

psiR=${psiR_values[$SLURM_ARRAY_TASK_ID]}
OUTDIR="${BASE_OUTDIR}/alpha_${ALPHA_NG}/psiR_${psiR}"

mkdir -p "$OUTDIR"

echo "Starting isotropy non-Gaussian full run"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "psiR: ${psiR}"
echo "ALPHA_NG: ${ALPHA_NG}"
echo "Output directory: ${OUTDIR}"
date

ALPHA_NG="$ALPHA_NG" Rscript "$SCRIPT" "$psiR" "$OUTDIR"

status=$?
date

echo "Finished isotropy non-Gaussian full run for psiR=${psiR}, ALPHA_NG=${ALPHA_NG}"
exit $status

