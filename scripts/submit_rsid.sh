#!/bin/bash
#SBATCH --job-name=iter003_rsid
#SBATCH --partition=orfoz
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=64G
#SBATCH --output=logs/slurm_rsid_%j.out
#SBATCH --error=logs/slurm_rsid_%j.err

set -eo pipefail
PROJECT=/arf/scratch/ycicek/iter_003_mtdna_psy_mr
cd "$PROJECT"
export ITER003_HOME="$PROJECT"
export PATH=$HOME/bin:$PATH
export MAMBA_ROOT_PREFIX=$HOME/.conda
eval "$(micromamba shell hook --shell bash)"
micromamba activate iter003_mr

echo "[$(date)] === rsID lookup build START ==="
Rscript scripts/01c_build_rsid_lookup.R
echo "[$(date)] === DONE ==="
