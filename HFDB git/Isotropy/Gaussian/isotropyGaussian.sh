#!/bin/bash
set -e

mkdir -p logs
sbatch run_isotropy_array_fullrun.slurm
