#!/bin/bash
#SBATCH --job-name=iter003_c4b_susie
#SBATCH --partition=debug
#SBATCH --time=03:55:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/slurm_c4b_susie_%j.out
#SBATCH --error=logs/slurm_c4b_susie_%j.err
set -eo pipefail
cd /arf/scratch/ycicek/iter_003_mtdna_psy_mr
export PATH=$HOME/bin:$PATH; export MAMBA_ROOT_PREFIX=$HOME/.conda
eval "$(micromamba shell hook --shell bash)"; micromamba activate iter003_mr
Rscript scripts/08d_coloc_susie_v2.R 2>&1 | tee logs/c4b_susie_${SLURM_JOB_ID}.log
