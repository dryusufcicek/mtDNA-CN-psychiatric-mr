#!/bin/bash
# ============================================================================
# 09_magma_mitocarta.sh — Task 3 (MAGMA gene-level + MitoCarta gene-set)
#
# Pipeline:
#  1. Build SNP-P inputs (python helper 09a)
#  2. SNP-to-gene annotation (MAGMA --annotate, NCBI37 gene_loc, 1KG EUR ref)
#  3. Gene-level analysis (MAGMA --gene-results) per trait
#  4. Gene-set analysis with MitoCarta3 (--gene-results + --set-annot)
#
# 7 traits total: Longchamps, Chong, Gupta_raw, Gupta_adj, ASD, ADHD, F3
# ============================================================================
set -eo pipefail

PROJECT=/arf/scratch/ycicek/iter_003_mtdna_psy_mr
MAGMA=/arf/home/ycicek/bin/magma_v1.10
GENE_LOC=$PROJECT/data/magma_aux/NCBI37.3.gene.loc
PLINK_REF=/arf/home/ycicek/EVOSCZ/data/ldsc/sldsc_ref/1000G_EUR_Phase3_plink/1000G.EUR.QC.merged
INPUT_DIR=$PROJECT/data/magma_inputs
OUT_DIR=$PROJECT/results/magma
SET_ANNOT=$INPUT_DIR/mitocarta3_set_annot.tsv

mkdir -p $OUT_DIR

ts() { date +"%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(ts)] $*"; }

log "=== 09_magma_mitocarta.sh START ==="
log "MAGMA: $MAGMA"
log "Gene_loc: $GENE_LOC"
log "PLINK ref: $PLINK_REF"

# ---------------------------------------------------------------------
# Step 0: prep inputs via Python helper (v11 env has pandas + openpyxl)
# ---------------------------------------------------------------------
log "--- Step 0: prep MAGMA inputs (python helper) ---"
eval "$(micromamba shell hook --shell bash)"
micromamba activate v11
pip install xlrd --quiet >/dev/null 2>&1 || true
python $PROJECT/scripts/09a_prep_magma_inputs.py 2>&1 | tee $PROJECT/logs/magma_prep_${SLURM_JOB_ID:-local}.log
micromamba deactivate

# ---------------------------------------------------------------------
# Step 1: SNP -> gene annotation (one-time, shared across traits)
# Window: 35kb upstream + 10kb downstream (common GWAS gene-mapping)
# ---------------------------------------------------------------------
ANNOT_PREFIX=$OUT_DIR/snp2gene_NCBI37
if [ ! -f "${ANNOT_PREFIX}.genes.annot" ]; then
  log "--- Step 1: SNP -> gene annotation ---"
  # Need a SNPLOC file (one trait's SNP-CHR-BP suffices as reference SNP map)
  # Use 1KG bim as canonical SNP location (covers EUR Phase3 SNPs)
  TMP_SNPLOC=$OUT_DIR/1kg_snploc.tsv
  if [ ! -f "$TMP_SNPLOC" ]; then
    awk 'BEGIN{OFS="\t"; print "SNP","CHR","BP"} {print $2,$1,$4}' \
      ${PLINK_REF}.bim > $TMP_SNPLOC
    log "  Built SNPLOC from 1KG merged bim: $(wc -l < $TMP_SNPLOC) SNPs"
  fi
  $MAGMA --annotate window=35,10 \
    --snp-loc $TMP_SNPLOC \
    --gene-loc $GENE_LOC \
    --out $ANNOT_PREFIX 2>&1 | tail -10
  log "  Annot file: ${ANNOT_PREFIX}.genes.annot"
else
  log "--- Step 1: annot already exists, skip ---"
fi

# ---------------------------------------------------------------------
# Step 2 + 3: gene-level + gene-set per trait
# ---------------------------------------------------------------------
TRAITS=(Longchamps2022 Chong2022 Gupta2023_raw Gupta2023_adj \
        ASD_Grove2019 ADHD_Demontis2023 CDG3_F3_Neurodev)

# Sample sizes per trait (used for gene analysis when --gene-model snp-wise)
declare -A NTRAIT
NTRAIT[Longchamps2022]=465809
NTRAIT[Chong2022]=295150
NTRAIT[Gupta2023_raw]=395000
NTRAIT[Gupta2023_adj]=395000
NTRAIT[ASD_Grove2019]=46350
NTRAIT[ADHD_Demontis2023]=225534
NTRAIT[CDG3_F3_Neurodev]=500000

for trait in "${TRAITS[@]}"; do
  log "--- Step 2 ($trait): gene-level analysis ---"
  GENE_OUT=$OUT_DIR/${trait}.gene
  SNP_P=$INPUT_DIR/${trait}.snp_p.tsv
  if [ ! -f "$SNP_P" ]; then
    log "  [WARN] SNP_P not found: $SNP_P -- skip"
    continue
  fi
  if [ -f "${GENE_OUT}.genes.raw" ]; then
    log "  Skip (exists): ${GENE_OUT}.genes.raw"
  else
    $MAGMA \
      --bfile $PLINK_REF \
      --gene-annot ${ANNOT_PREFIX}.genes.annot \
      --pval $SNP_P use=SNP,P ncol=N \
      --gene-model snp-wise=mean \
      --out $GENE_OUT 2>&1 | tail -20
  fi

  log "--- Step 3 ($trait): MitoCarta set enrichment ---"
  SET_OUT=$OUT_DIR/${trait}.mitoset
  if [ -f "${SET_OUT}.gsa.out" ]; then
    log "  Skip (exists): ${SET_OUT}.gsa.out"
  else
    $MAGMA \
      --gene-results ${GENE_OUT}.genes.raw \
      --set-annot $SET_ANNOT \
      --out $SET_OUT 2>&1 | tail -10
  fi
done

# ---------------------------------------------------------------------
# Step 4: combine results into one TSV
# ---------------------------------------------------------------------
log "--- Step 4: aggregate set enrichment results ---"
SUMM=$OUT_DIR/MAGMA_mitocarta_summary.tsv
echo -e "trait\tfull_set_name\tn_genes\tbeta\tse\tp\tq_bonf" > $SUMM
for trait in "${TRAITS[@]}"; do
  RESULT=$OUT_DIR/${trait}.mitoset.gsa.out
  if [ ! -f "$RESULT" ]; then
    log "  [WARN] gsa.out missing for $trait"
    continue
  fi
  # gsa.out format: VARIABLE TYPE NGENES BETA BETA_STD SE P
  # Skip comment lines (#) and header; emit trait + name + ngenes + beta + se + p
  awk -v t="$trait" 'BEGIN{OFS="\t"} !/^#/ && NR>1 && $1 != "VARIABLE" {
    print t, $1, $3, $4, $6, $7
  }' $RESULT > /tmp/.row_${trait}_$$
  # Sort by p, apply Bonferroni based on number of sets per trait
  Nsets=$(wc -l < /tmp/.row_${trait}_$$)
  awk -v ns=$Nsets 'BEGIN{OFS="\t"} {
    qb = $6 * ns; if (qb > 1) qb = 1;
    print $1, $2, $3, $4, $5, $6, qb
  }' /tmp/.row_${trait}_$$ >> $SUMM
  rm /tmp/.row_${trait}_$$
done
log "  Summary -> $SUMM"
log "  Lines (incl header): $(wc -l < $SUMM)"
# Top 5 hits per trait
log ""
log "=== Top hits (p<0.01) per trait ==="
awk -F'\t' 'NR>1 && $6 < 0.01' $SUMM | sort -t$'\t' -k1,1 -k6,6g | head -30 || true

log "=== 09_magma_mitocarta.sh END ==="
