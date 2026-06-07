#!/bin/bash
#SBATCH --job-name=iter003_mrlap
#SBATCH --partition=orfoz
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=64G
#SBATCH --output=logs/slurm_mrlap_%j.out
#SBATCH --error=logs/slurm_mrlap_%j.err

set -eo pipefail

PROJECT=/arf/scratch/ycicek/iter_003_mtdna_psy_mr
cd "$PROJECT"
export ITER003_HOME="$PROJECT"

# Activate iter003_mr R env
export PATH=$HOME/bin:$PATH
export MAMBA_ROOT_PREFIX=$HOME/.conda
eval "$(micromamba shell hook --shell bash)"
micromamba activate iter003_mr

echo "[$(date)] === MRlap sample-overlap correction START ==="
echo "Host: $(hostname)"
echo "Job: $SLURM_JOB_ID"
echo "Mem: $(free -g | grep Mem)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Disk: $(df -h /arf/scratch/ycicek | tail -1)"

# Sanity check packages
echo ""
echo "[$(date)] === Package check ==="
Rscript -e 'for (p in c("MRlap","TwoSampleMR","data.table","dplyr","GenomicSEM")) {
  cat(sprintf("  %-15s %s\n", p, ifelse(requireNamespace(p, quietly=TRUE), "OK", "MISSING")))
}'

# Sanity check LD reference
echo ""
echo "[$(date)] === LD reference check ==="
LD_DIR=$PROJECT/data/ldsc_mrlap/eur_w_ld_chr
n_ldscore=$(ls $LD_DIR/*.l2.ldscore.gz 2>/dev/null | wc -l)
n_m=$(ls $LD_DIR/*.l2.M 2>/dev/null | wc -l)
n_m550=$(ls $LD_DIR/*.l2.M_5_50 2>/dev/null | wc -l)
echo "  ldscore.gz files: $n_ldscore (expect 22)"
echo "  M files:           $n_m (expect 22)"
echo "  M_5_50 files:      $n_m550 (expect 22)"

echo ""
echo "[$(date)] === Running 06_mrlap.R ==="
Rscript scripts/06_mrlap.R 2>&1 | tee logs/mrlap_${SLURM_JOB_ID}.log

echo ""
echo "[$(date)] === MRlap PIPELINE COMPLETE ==="
echo "Disk after: $(df -h /arf/scratch/ycicek | tail -1)"
ls -la results/mrlap/ 2>/dev/null || echo "results/mrlap not found"
