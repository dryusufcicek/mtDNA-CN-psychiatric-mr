#!/usr/bin/env Rscript
# ============================================================================
# 08_coloc_brain.R — Task 2 (Coloc with GTEx v10 brain eQTL, 13 tissues)
#
# Goal: For protective MR signals (ASD, ADHD, F3), test whether the lead
# instrument SNPs colocalize with brain cis-eQTL signals -> mechanistic
# anchoring at gene/tissue level.
#
# Inputs:
#   results/forward_mr/all_forward_results.rds (read SNP-level dat)
#   data/outcomes/{ASD_Grove2019.gz, ADHD_Demontis2023.meta.gz, CDG3_F3_Neurodev.tsv.gz}
#   data/coloc_inputs/cis_eqtl_<TISSUE>.tsv  (built by 08a python helper)
#
# Outputs:
#   results/coloc/lead_loci.tsv                — lead SNPs +/- 500kb regions
#   results/coloc/COLOC_per_pair.tsv           — long: region x tissue x gene x PP*
#   results/coloc/COLOC_summary.tsv            — strong (PP4>0.75) / suggestive (PP4>0.5)
# ============================================================================

suppressMessages({
  library(data.table)
  library(dplyr)
  library(coloc)
})

# -------- robustness lib (log_step, log_info, log_warn) --------
PROJECT <- Sys.getenv("ITER003_HOME",
                      "/arf/scratch/ycicek/iter_003_mtdna_psy_mr")
setwd(PROJECT)
source(file.path(PROJECT, "scripts", "lib_robust.R"))

OUT_DIR    <- file.path(PROJECT, "results", "coloc")
COLOC_IN   <- file.path(PROJECT, "data", "coloc_inputs")
FORWARD    <- file.path(PROJECT, "results", "forward_mr",
                        "all_forward_results.rds")
OUT_LOCI   <- file.path(OUT_DIR, "lead_loci.tsv")
OUT_PAIR   <- file.path(OUT_DIR, "COLOC_per_pair.tsv")
OUT_SUMM   <- file.path(OUT_DIR, "COLOC_summary.tsv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(COLOC_IN, showWarnings = FALSE, recursive = TRUE)

# Outcome metadata: case/control N + sample type
# (needed for coloc.abf type='cc' or 'quant')
OUTCOME_META <- list(
  ASD  = list(file = "data/outcomes/ASD_Grove2019.gz",            n_case = 18381,  n_control = 27969,
              type = "cc", s = 18381/(18381+27969)),
  ADHD = list(file = "data/outcomes/ADHD_Demontis2023.meta.gz",   n_case = 38691,  n_control = 186843,
              type = "cc", s = 38691/(38691+186843)),
  F3   = list(file = "data/outcomes/CDG3_F3_Neurodev.tsv.gz",     n_case = NA,     n_control = NA,
              type = "quant", N = 500000)  # CDG3 factor -- approximate effective N
)

# GTEx v10 brain tissues (13)
BRAIN_TISSUES <- c(
  "Brain_Amygdala",
  "Brain_Anterior_cingulate_cortex_BA24",
  "Brain_Caudate_basal_ganglia",
  "Brain_Cerebellar_Hemisphere",
  "Brain_Cerebellum",
  "Brain_Cortex",
  "Brain_Frontal_Cortex_BA9",
  "Brain_Hippocampus",
  "Brain_Hypothalamus",
  "Brain_Nucleus_accumbens_basal_ganglia",
  "Brain_Putamen_basal_ganglia",
  "Brain_Spinal_cord_cervical_c-1",
  "Brain_Substantia_nigra"
)

# GTEx v10 brain median N (https://gtexportal.org -- v10 release notes)
# Values approximate, used as N for the eQTL track
GTEX_N <- c(
  Brain_Amygdala = 152,
  Brain_Anterior_cingulate_cortex_BA24 = 176,
  Brain_Caudate_basal_ganglia = 246,
  Brain_Cerebellar_Hemisphere = 215,
  Brain_Cerebellum = 250,
  Brain_Cortex = 255,
  Brain_Frontal_Cortex_BA9 = 209,
  Brain_Hippocampus = 197,
  Brain_Hypothalamus = 202,
  Brain_Nucleus_accumbens_basal_ganglia = 246,
  Brain_Putamen_basal_ganglia = 205,
  `Brain_Spinal_cord_cervical_c-1` = 159,
  Brain_Substantia_nigra = 140
)

REGION_HALF <- 500000  # +/- 500 kb

# ============================================================================
# Step 1 — Build lead loci from forward MR significant exposure-outcome pairs
# ============================================================================
log_step("Step 1: build lead loci from forward MR")
res <- readRDS(FORWARD)
log_info(sprintf("Loaded %d forward MR pair results", length(res)))

# Significant outcomes: ASD, ADHD, F3 (across all 4 exposures)
EXPOSURES <- c("Longchamps_2022", "Chong_2022", "Gupta_2023_raw", "Gupta_2023_adj")
OUTCOMES  <- c("ASD", "ADHD", "F3")

target_pairs <- as.vector(outer(EXPOSURES, OUTCOMES,
                                FUN = function(e, o) paste0(e, "__", o)))
target_pairs <- target_pairs[target_pairs %in% names(res)]
log_info(sprintf("Target pairs: %d", length(target_pairs)))

# Aggregate harmonized SNPs across all target pairs
all_dat <- bind_rows(lapply(target_pairs, function(p) {
  d <- res[[p]]$dat
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d$pair <- p
  d
}))

if (is.null(all_dat) || nrow(all_dat) == 0) {
  log_error("No harmonized data found for target pairs")
  quit(status = 1)
}

log_info(sprintf("Aggregated harmonized SNPs: %d rows (%d unique SNPs)",
                 nrow(all_dat), length(unique(all_dat$SNP))))

# Per-SNP minimum P across the 3 outcomes (used for lead selection)
lead_snps_raw <- all_dat %>%
  group_by(SNP) %>%
  summarise(
    chr     = first(chr.outcome),
    pos     = first(pos.outcome),
    min_p_outcome = suppressWarnings(min(pval.outcome, na.rm = TRUE)),
    n_pairs = n_distinct(pair),
    .groups = "drop"
  ) %>%
  filter(!is.na(chr), !is.na(pos))

log_info(sprintf("Unique candidate SNPs with position: %d", nrow(lead_snps_raw)))

# Lead selection: SNPs nominally significant in at least one psy outcome
# AND a top set per chromosome to limit region count
lead_snps <- lead_snps_raw %>%
  filter(min_p_outcome <= 0.05) %>%   # nominal sig in 1+ outcome
  arrange(min_p_outcome) %>%
  slice_head(n = 25)                  # cap at top 25 to keep run tractable

# If <5 SNPs survive (under-powered selection), fall back to top 15 by outcome P
if (nrow(lead_snps) < 5) {
  log_warn(sprintf("Only %d SNPs at outcome P<=0.05; falling back to top-15 by outcome P",
                   nrow(lead_snps)))
  lead_snps <- lead_snps_raw %>%
    arrange(min_p_outcome) %>%
    slice_head(n = 15)
}

log_info(sprintf("Selected %d lead SNPs for coloc", nrow(lead_snps)))

# Collapse nearby leads (<1 Mb same chr) into single regions
lead_snps <- lead_snps %>% arrange(chr, pos)
lead_snps$region_id <- NA_character_
cur_id <- 0L
prev_chr <- -1; prev_end <- -1L
for (i in seq_len(nrow(lead_snps))) {
  if (lead_snps$chr[i] != prev_chr ||
      (lead_snps$pos[i] - REGION_HALF) > prev_end) {
    cur_id <- cur_id + 1L
    prev_chr <- lead_snps$chr[i]
  }
  lead_snps$region_id[i] <- sprintf("R%03d", cur_id)
  prev_end <- max(prev_end, lead_snps$pos[i] + REGION_HALF)
}

# Build region boundaries for parquet preprocess
regions <- lead_snps %>%
  group_by(region_id, CHR = chr) %>%
  summarise(
    POS_START = min(pos) - REGION_HALF,
    POS_END   = max(pos) + REGION_HALF,
    n_leads   = n(),
    lead_snp  = paste(SNP, collapse = ";"),
    .groups   = "drop"
  )
fwrite(regions, OUT_LOCI, sep = "\t")
log_info(sprintf("Wrote %d regions -> %s", nrow(regions), OUT_LOCI))

# Also save mapping for downstream
fwrite(lead_snps, file.path(OUT_DIR, "lead_snps_full.tsv"), sep = "\t")

# ============================================================================
# Step 2 — Trigger python extraction of cis-eQTL parquet rows for each region
# ============================================================================
log_step("Step 2: extract cis-eQTL parquet rows (Python pyarrow)")
PY_HELPER <- file.path(PROJECT, "scripts", "08a_extract_cis_eqtl.py")
if (!file.exists(PY_HELPER)) {
  log_error(sprintf("Python helper not found: %s", PY_HELPER))
  quit(status = 1)
}

# Use the v11 micromamba env (has pyarrow). Spawn shell-isolated session.
cmd <- sprintf(
  "bash -c 'export PATH=$HOME/bin:$PATH; export MAMBA_ROOT_PREFIX=$HOME/.conda; eval \"$(micromamba shell hook --shell bash)\"; micromamba activate v11; python %s %s %s'",
  PY_HELPER, OUT_LOCI, COLOC_IN
)
log_info(sprintf("Running: %s", cmd))
ret <- system(cmd)
if (ret != 0) {
  log_error(sprintf("Python extraction failed (exit %d)", ret))
  quit(status = 1)
}

# ============================================================================
# Step 3 — Build OUTCOME GWAS subset per region (for coloc dataset2)
# ============================================================================
log_step("Step 3: build outcome GWAS subsets per region")

read_outcome <- function(meta) {
  log_info(sprintf("Reading outcome: %s", meta$file))
  d <- fread(meta$file)
  # Standardize columns -- different headers per outcome
  if ("OR" %in% colnames(d) && !"BETA" %in% colnames(d)) {
    d[, BETA := log(OR)]
  }
  if (!"BP" %in% colnames(d)) {
    pos_col <- intersect(c("POS", "BP", "POSITION"), colnames(d))[1]
    if (!is.na(pos_col)) setnames(d, pos_col, "BP")
  }
  if (!"CHR" %in% colnames(d)) {
    chr_col <- intersect(c("CHROM", "Chromosome", "chr"), colnames(d))[1]
    if (!is.na(chr_col)) setnames(d, chr_col, "CHR")
  }
  d$CHR <- as.integer(gsub("chr", "", as.character(d$CHR)))
  d
}

outcome_data <- list()
for (oname in names(OUTCOME_META)) {
  outcome_data[[oname]] <- read_outcome(OUTCOME_META[[oname]])
  log_info(sprintf("  %s: %d SNPs loaded", oname, nrow(outcome_data[[oname]])))
}

# ============================================================================
# Step 4 — Run coloc.abf per (region x tissue x gene) triplet x outcome
# ============================================================================
log_step("Step 4: run coloc.abf")

coloc_results <- list()
slot_i <- 0L

for (tissue in BRAIN_TISSUES) {
  eqtl_path <- file.path(COLOC_IN, sprintf("cis_eqtl_%s.tsv", tissue))
  if (!file.exists(eqtl_path)) {
    log_warn(sprintf("eQTL file missing: %s -- skip", eqtl_path))
    next
  }
  eqtl <- fread(eqtl_path)
  if (nrow(eqtl) == 0) {
    log_warn(sprintf("Empty eQTL table for %s", tissue))
    next
  }
  log_info(sprintf("[%s] eQTL rows=%d, genes=%d",
                   tissue, nrow(eqtl), eqtl[, uniqueN(gene_id)]))
  # GTEx v10 has variant_pos in GRCh38 -- outcome GWAS uses GRCh37/38 depending
  # We rely on rs ID where available. For coloc we'll use chr+variant_pos
  # matching to outcome BP. WARNING: outcome may be GRCh37; we'll match by rs ID
  # via rs_id_dbSNP155_GRCh38p13 column if present, else by chr+pos (lossy).
  # Default safer path: match by SNP rsid.

  # Sample size for eQTL track (population N for GTEx tissue)
  n_eqtl <- GTEX_N[tissue]
  if (is.null(n_eqtl) || is.na(n_eqtl)) n_eqtl <- 200

  for (rid in unique(eqtl$region_id)) {
    sub_eqtl <- eqtl[region_id == rid]
    region_row <- regions[regions$region_id == rid, ]
    chr_i <- region_row$CHR[1]
    pos_lo <- region_row$POS_START[1]
    pos_hi <- region_row$POS_END[1]

    for (gene in unique(sub_eqtl$gene_id)) {
      e <- sub_eqtl[gene_id == gene]
      # Need >= 50 SNPs in region for sensible coloc; otherwise skip
      if (nrow(e) < 30) next

      # eQTL coloc dataset
      d1 <- list(
        type     = "quant",
        beta     = e$slope,
        varbeta  = e$slope_se^2,
        snp      = if ("rs_id_dbSNP155_GRCh38p13" %in% colnames(e))
                     e$rs_id_dbSNP155_GRCh38p13 else e$variant_id,
        position = e$variant_pos,
        MAF      = pmin(e$af, 1 - e$af),
        N        = n_eqtl
      )
      # Drop SNPs with NA rsid or NA stats
      keep1 <- !is.na(d1$snp) & d1$snp != "." & !is.na(d1$beta) &
               !is.na(d1$varbeta) & d1$varbeta > 0 & !is.na(d1$MAF) &
               d1$MAF > 0 & d1$MAF < 1
      if (sum(keep1) < 30) next
      d1 <- lapply(d1, function(x) if (length(x) == 1) x else x[keep1])

      for (oname in names(OUTCOME_META)) {
        omd <- outcome_data[[oname]]
        meta <- OUTCOME_META[[oname]]
        # Outcome subset matching by SNP (rs ID)
        og <- omd[SNP %in% d1$snp]
        if (nrow(og) < 30) next

        # Filter to common SNPs
        common <- intersect(d1$snp, og$SNP)
        if (length(common) < 30) next
        idx1 <- match(common, d1$snp)
        og2 <- og[match(common, og$SNP)]

        # Build outcome dataset for coloc
        if (meta$type == "cc") {
          d2 <- list(
            type    = "cc",
            beta    = if ("BETA" %in% colnames(og2)) og2$BETA else log(og2$OR),
            varbeta = og2$SE^2,
            snp     = og2$SNP,
            position = og2$BP,
            s       = meta$s,
            N       = meta$n_case + meta$n_control,
            MAF     = if ("FRQ_U" %in% colnames(og2)) og2$FRQ_U else 0.2  # placeholder when missing
          )
        } else {
          d2 <- list(
            type    = "quant",
            beta    = og2$BETA,
            varbeta = og2$SE^2,
            snp     = og2$SNP,
            position = og2$BP,
            N       = meta$N,
            MAF     = if ("MAF" %in% colnames(og2)) og2$MAF else 0.2
          )
        }
        keep2 <- !is.na(d2$beta) & !is.na(d2$varbeta) & d2$varbeta > 0
        if (sum(keep2) < 30) next
        d2 <- lapply(d2, function(x) if (length(x) <= 1) x else x[keep2])
        # Recompute common SNPs after keep2
        common2 <- intersect(d1$snp, d2$snp)
        if (length(common2) < 30) next
        idx1b <- match(common2, d1$snp)
        idx2b <- match(common2, d2$snp)
        d1b <- lapply(d1, function(x) if (length(x) <= 1) x else x[idx1b])
        d2b <- lapply(d2, function(x) if (length(x) <= 1) x else x[idx2b])

        cres <- tryCatch(
          coloc::coloc.abf(d1b, d2b),
          error = function(e) { log_warn(sprintf("coloc fail %s/%s/%s: %s",
                                                 rid, tissue, gene,
                                                 conditionMessage(e)));
                                NULL }
        )
        if (is.null(cres)) next
        slot_i <- slot_i + 1L
        coloc_results[[slot_i]] <- data.frame(
          region_id = rid,
          chr       = chr_i,
          tissue    = tissue,
          gene_id   = gene,
          outcome   = oname,
          n_snps    = length(common2),
          PP0       = cres$summary[["PP.H0.abf"]],
          PP1       = cres$summary[["PP.H1.abf"]],
          PP2       = cres$summary[["PP.H2.abf"]],
          PP3       = cres$summary[["PP.H3.abf"]],
          PP4       = cres$summary[["PP.H4.abf"]],
          PP4_PP3_ratio = cres$summary[["PP.H4.abf"]] / max(cres$summary[["PP.H3.abf"]], 1e-30),
          stringsAsFactors = FALSE
        )
        if (slot_i %% 50 == 0) {
          log_info(sprintf("  ... ran %d coloc tests so far", slot_i))
        }
      }
    }
  }
}

if (length(coloc_results) == 0) {
  log_error("No coloc tests succeeded -- aborting")
  quit(status = 1)
}

coloc_df <- bind_rows(coloc_results)
fwrite(coloc_df, OUT_PAIR, sep = "\t")
log_info(sprintf("Wrote per-pair coloc -> %s (%d rows)",
                 OUT_PAIR, nrow(coloc_df)))

# ============================================================================
# Step 5 — Summary: strong (PP4>0.75) + suggestive (PP4>0.5)
# ============================================================================
log_step("Step 5: summary tables")

# Add gene symbol from GTEx eGenes file (gene_id -> gene_name)
eg_path <- "/arf/scratch/ycicek/h2_paper2/data/gtex_v10_eqtl_brain/Brain_Cortex.v10.eGenes.txt.gz"
egene <- fread(eg_path, select = c("gene_id", "gene_name"))
coloc_df <- coloc_df %>%
  left_join(distinct(egene), by = "gene_id")

summary_strong <- coloc_df %>%
  filter(PP4 >= 0.75) %>%
  arrange(desc(PP4))
summary_sugg <- coloc_df %>%
  filter(PP4 >= 0.5, PP4 < 0.75) %>%
  arrange(desc(PP4))

OUT_STRONG <- file.path(OUT_DIR, "COLOC_strong_PP4_075.tsv")
OUT_SUGG <- file.path(OUT_DIR, "COLOC_suggestive_PP4_050.tsv")
fwrite(summary_strong, OUT_STRONG, sep = "\t")
fwrite(summary_sugg,   OUT_SUGG, sep = "\t")

# Top per outcome
top_per_outcome <- coloc_df %>%
  group_by(outcome) %>%
  arrange(desc(PP4)) %>%
  slice_head(n = 20) %>%
  ungroup()
fwrite(top_per_outcome, OUT_SUMM, sep = "\t")

log_info(sprintf("Strong coloc (PP4>=0.75): %d -> %s",
                 nrow(summary_strong), OUT_STRONG))
log_info(sprintf("Suggestive (0.5<=PP4<0.75): %d -> %s",
                 nrow(summary_sugg), OUT_SUGG))
log_info(sprintf("Top-20 per outcome -> %s", OUT_SUMM))

# Print top 5 strong if any
if (nrow(summary_strong) > 0) {
  cat("\n=== Top strong coloc (PP4>=0.75) ===\n")
  print(head(summary_strong %>% select(outcome, tissue, gene_name, gene_id,
                                       PP4, PP3, n_snps), 5))
}
cat("\n=== Per-outcome best PP4 ===\n")
print(coloc_df %>% group_by(outcome) %>%
        summarise(max_pp4 = max(PP4),
                  n_strong = sum(PP4 >= 0.75),
                  n_sugg = sum(PP4 >= 0.5 & PP4 < 0.75),
                  .groups = "drop"))

log_step("DONE -- 08_coloc_brain.R")
