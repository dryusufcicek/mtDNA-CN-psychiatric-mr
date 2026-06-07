#!/bin/bash
#SBATCH --job-name=iter003_qc
#SBATCH --partition=orfoz
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=64G
#SBATCH --output=logs/slurm_qc_%j.out
#SBATCH --error=logs/slurm_qc_%j.err

set -eo pipefail

PROJECT=/arf/scratch/ycicek/iter_003_mtdna_psy_mr
cd "$PROJECT"
export ITER003_HOME="$PROJECT"

# Activate R env
export PATH=$HOME/bin:$PATH
export MAMBA_ROOT_PREFIX=$HOME/.conda
eval "$(micromamba shell hook --shell bash)"
micromamba activate iter003_mr

echo "[$(date)] === START QC ==="
echo "Host: $(hostname)"
echo "Mem available: $(free -g | grep Mem)"

Rscript scripts/01b_data_qc.R

echo "[$(date)] === DONE ==="
