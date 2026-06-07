#!/bin/bash
#SBATCH --job-name=iter003_c3_3way
#SBATCH --partition=debug
#SBATCH --time=03:55:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/slurm_c3_3way_%j.out
#SBATCH --error=logs/slurm_c3_3way_%j.err
set -eo pipefail
cd /arf/scratch/ycicek/iter_003_mtdna_psy_mr
export PATH=$HOME/bin:$PATH; export MAMBA_ROOT_PREFIX=$HOME/.conda
eval "$(micromamba shell hook --shell bash)"; micromamba activate iter003_mr
echo "[$(date)] three-way coloc (eQTL x mtDNA-CN x disorder) start"
Rscript scripts/08c_threeway_coloc.R 2>&1 | tee logs/c3_3way_${SLURM_JOB_ID}.log
echo "[$(date)] done"
