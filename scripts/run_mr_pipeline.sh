#!/bin/bash
#SBATCH --job-name=iter003_mr
#SBATCH --partition=orfoz
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err

# iter_003 — full MR pipeline on TRUBA
# Submit: sbatch scripts/run_mr_pipeline.sh
# Date: 2026-05-22

set -euo pipefail

PROJECT=$HOME/scz-research/iter_003_mtdna_psy_mr
cd "$PROJECT"

# Activate R env
export PATH=$HOME/bin:$PATH
eval "$(micromamba shell hook --shell bash)"
micromamba activate iter003_mr

# Verify R packages
echo "[$(date)] === R version ==="
R --version | head -1

echo "[$(date)] === Package check ==="
Rscript -e 'for (p in c("TwoSampleMR","MRPRESSO","RadialMR","metafor","data.table","dplyr","ggplot2","MendelianRandomization")) cat(sprintf("%-25s %s\n", p, ifelse(requireNamespace(p,quietly=TRUE),"OK","MISSING")))'

# Pipeline steps
echo "[$(date)] === Step 01: Header inspection ==="
Rscript scripts/01_inspect_headers.R

echo "[$(date)] === Step 02: Harmonize exposures ==="
Rscript scripts/02_harmonize_exposures.R

echo "[$(date)] === Step 03: Forward MR (mtDNA-CN -> psy disorders) ==="
Rscript scripts/03_forward_mr.R

echo "[$(date)] === Step 04: Cross-panel synthesis (meta + design-effect correction) ==="
Rscript scripts/04_meta_mr.R

echo "[$(date)] === DONE ==="
