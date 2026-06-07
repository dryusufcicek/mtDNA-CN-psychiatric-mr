#!/bin/bash
#SBATCH --job-name=supp_t1
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/arf/scratch/ycicek/iter_003_mtdna_psy_mr/scripts/supp_t1_%j.out
#SBATCH --error=/arf/scratch/ycicek/iter_003_mtdna_psy_mr/scripts/supp_t1_%j.err

set -euo pipefail
python3 -u /arf/scratch/ycicek/iter_003_mtdna_psy_mr/scripts/build_supp_table_1.py
echo "RETCODE=$?"
ls -la /arf/scratch/ycicek/iter_003_mtdna_psy_mr/manuscript_draft/supplementary/
