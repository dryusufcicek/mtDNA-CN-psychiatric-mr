#!/bin/bash
#SBATCH --job-name=iter003_coloc
#SBATCH --partition=orfoz
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=64G
#SBATCH --output=logs/slurm_coloc_%j.out
#SBATCH --error=logs/slurm_coloc_%j.err

set -eo pipefail

PROJECT=/arf/scratch/ycicek/iter_003_mtdna_psy_mr
cd "$PROJECT"
export ITER003_HOME="$PROJECT"

export PATH=$HOME/bin:$PATH
export MAMBA_ROOT_PREFIX=$HOME/.conda
eval "$(micromamba shell hook --shell bash)"
micromamba activate iter003_mr

echo "[$(date)] === COLOC START ==="
echo "Host: $(hostname)"
echo "Job: $SLURM_JOB_ID"
echo "Mem: $(free -g | grep Mem)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Disk: $(df -h /arf/scratch/ycicek | tail -1)"

echo ""
echo "[$(date)] === Package check ==="
Rscript -e 'for (p in c("coloc","data.table","dplyr","TwoSampleMR")) {
  cat(sprintf("  %-15s %s\n", p, ifelse(requireNamespace(p, quietly=TRUE), "OK", "MISSING")))
}'

echo ""
echo "[$(date)] === Pyarrow env (v11) check ==="
micromamba deactivate
micromamba activate v11
python -c "import pyarrow, pandas; print('pyarrow', pyarrow.__version__, 'pandas', pandas.__version__)"
micromamba deactivate
micromamba activate iter003_mr

echo ""
echo "[$(date)] === Running 08_coloc_brain.R ==="
Rscript scripts/08_coloc_brain.R 2>&1 | tee logs/coloc_${SLURM_JOB_ID}.log

echo ""
echo "[$(date)] === Final disk ==="
df -h /arf/scratch/ycicek | tail -1

echo ""
echo "[$(date)] === COLOC END ==="
