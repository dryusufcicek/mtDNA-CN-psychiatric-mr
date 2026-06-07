#!/bin/bash
#SBATCH --job-name=iter003_magma
#SBATCH --partition=orfoz
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=64G
#SBATCH --output=logs/slurm_magma_%j.out
#SBATCH --error=logs/slurm_magma_%j.err

set -eo pipefail

PROJECT=/arf/scratch/ycicek/iter_003_mtdna_psy_mr
cd "$PROJECT"
export ITER003_HOME="$PROJECT"

export PATH=$HOME/bin:$PATH
export MAMBA_ROOT_PREFIX=$HOME/.conda

echo "[$(date)] === MAGMA START ==="
echo "Host: $(hostname)"
echo "Job: $SLURM_JOB_ID"
echo "Mem: $(free -g | grep Mem)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Disk: $(df -h /arf/scratch/ycicek | tail -1)"

# Sanity check binaries
echo ""
echo "[$(date)] === Binary check ==="
MAGMA=/arf/home/ycicek/bin/magma_v1.10
PLINK_REF=/arf/home/ycicek/EVOSCZ/data/ldsc/sldsc_ref/1000G_EUR_Phase3_plink/1000G.EUR.QC.merged
GENE_LOC=$PROJECT/data/magma_aux/NCBI37.3.gene.loc
MITOC=$PROJECT/data/magma_aux/Human.MitoCarta3.0.xls
for f in $MAGMA $GENE_LOC $MITOC ${PLINK_REF}.bed; do
  if [ -e "$f" ]; then
    echo "  OK: $f"
  else
    echo "  MISSING: $f"
    exit 1
  fi
done

echo ""
echo "[$(date)] === MAGMA version ==="
$MAGMA --version 2>&1 | head -3

echo ""
echo "[$(date)] === Running 09_magma_mitocarta.sh ==="
bash $PROJECT/scripts/09_magma_mitocarta.sh 2>&1 | tee logs/magma_${SLURM_JOB_ID}.log

echo ""
echo "[$(date)] === Final disk ==="
df -h /arf/scratch/ycicek | tail -1

echo ""
echo "[$(date)] === MAGMA END ==="
