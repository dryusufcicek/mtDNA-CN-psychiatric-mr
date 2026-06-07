#!/bin/bash
#SBATCH --job-name=iter003_c1_ldsc
#SBATCH --partition=orfoz
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=96G
#SBATCH --output=logs/slurm_c1_ldsc_%j.out
#SBATCH --error=logs/slurm_c1_ldsc_%j.err

set -eo pipefail
PROJECT=/arf/scratch/ycicek/iter_003_mtdna_psy_mr
cd "$PROJECT"
export PATH=$HOME/bin:$PATH; export MAMBA_ROOT_PREFIX=$HOME/.conda
eval "$(micromamba shell hook --shell bash)"
micromamba activate iter003_mr
mkdir -p data/ldsc_rho results/ldsc_rho

echo "[$(date)] build HM3 snplist (w_hm3.justrs ∩ g1000_eur.bim, with alleles)"
JUSTRS=/arf/scratch/ycicek/h2_phase1_strengthening/aim2_mixer/data/ref/w_hm3.justrs
BIM=/arf/scratch/ycicek/h2_paper2/data/g1000_eur/g1000_eur.bim
awk 'BEGIN{print "SNP\tA1\tA2"} NR==FNR{a[$1]=1; next} ($2 in a){print $2"\t"$5"\t"$6}' "$JUSTRS" "$BIM" > data/ldsc_rho/w_hm3.snplist
echo "snplist lines: $(wc -l < data/ldsc_rho/w_hm3.snplist)"

echo "[$(date)] run C1_ldsc_rho.R"
Rscript scripts/C1_ldsc_rho.R 2>&1 | tee logs/c1_ldsc_${SLURM_JOB_ID}.log
echo "[$(date)] DONE C1"
