#!/bin/bash
#SBATCH --job-name=iter003_rev_mr
#SBATCH --partition=orfoz
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=64G
#SBATCH --output=logs/slurm_reverse_mr_%j.out
#SBATCH --error=logs/slurm_reverse_mr_%j.err

set -eo pipefail

PROJECT=/arf/scratch/ycicek/iter_003_mtdna_psy_mr
cd "$PROJECT"
export ITER003_HOME="$PROJECT"

# Activate iter003_mr R env
export PATH=$HOME/bin:$PATH
export MAMBA_ROOT_PREFIX=$HOME/.conda
eval "$(micromamba shell hook --shell bash)"
micromamba activate iter003_mr

echo "[$(date)] === Reverse MR START ==="
echo "Host: $(hostname)"
echo "Job: $SLURM_JOB_ID"
echo "Mem: $(free -g | grep Mem)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Disk: $(df -h /arf/scratch/ycicek | tail -1)"

# Sanity check packages
echo ""
echo "[$(date)] === Package check ==="
Rscript -e 'for (p in c("TwoSampleMR","MRPRESSO","RadialMR","data.table","dplyr")) {
  cat(sprintf("  %-15s %s\n", p, ifelse(requireNamespace(p, quietly=TRUE), "OK", "MISSING")))
}'

# Sanity check PLINK + 1KG reference
echo ""
echo "[$(date)] === PLINK + 1KG check ==="
PLINK=/arf/home/ycicek/bin/plink
LD_REF_DIR=/arf/home/ycicek/EVOSCZ/data/ldsc/sldsc_ref/1000G_EUR_Phase3_plink
echo "  PLINK: $PLINK ($([ -x $PLINK ] && echo OK || echo MISSING))"
echo "  1KG ref dir: $LD_REF_DIR ($([ -d $LD_REF_DIR ] && echo OK || echo MISSING))"
echo "  Merged BED exists: $([ -e $LD_REF_DIR/1000G.EUR.QC.merged.bed ] && echo YES || echo NO_will_build)"

echo ""
echo "[$(date)] === Running 07_reverse_mr.R ==="
Rscript scripts/07_reverse_mr.R 2>&1 | tee logs/reverse_mr_${SLURM_JOB_ID}.log

echo ""
echo "[$(date)] === Reverse MR PIPELINE COMPLETE ==="
echo "Disk after: $(df -h /arf/scratch/ycicek | tail -1)"
ls -la results/reverse_mr/ 2>/dev/null || echo "results/reverse_mr not found"
