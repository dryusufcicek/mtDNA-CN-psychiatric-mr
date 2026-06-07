#!/bin/bash
#SBATCH --job-name=iter003_coloc2
#SBATCH --partition=orfoz
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=64G
#SBATCH --output=logs/slurm_coloc2_%j.out
#SBATCH --error=logs/slurm_coloc2_%j.err

set -eo pipefail
PROJECT=/arf/scratch/ycicek/iter_003_mtdna_psy_mr
cd "$PROJECT"; export ITER003_HOME="$PROJECT"
export PATH=$HOME/bin:$PATH; export MAMBA_ROOT_PREFIX=$HOME/.conda
eval "$(micromamba shell hook --shell bash)"

echo "[$(date)] === C2 COLOC (liftOver-fixed) START === job $SLURM_JOB_ID on $(hostname)"
echo "[$(date)] Step A: liftOver eQTL GRCh38->GRCh37 (env v11)"
micromamba activate v11
python scripts/08b_liftover_eqtl.py
micromamba deactivate

echo "[$(date)] Step B: coloc.abf chr:pos-matched, MitoCarta-flagged (env iter003_mr)"
micromamba activate iter003_mr
Rscript scripts/08b_coloc_fixed.R 2>&1 | tee logs/coloc2_${SLURM_JOB_ID}.log
echo "[$(date)] === C2 COLOC END ==="
