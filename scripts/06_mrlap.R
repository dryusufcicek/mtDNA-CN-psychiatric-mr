#!/usr/bin/env Rscript
# iter_003 — Step 06: MRlap sample-overlap corrected MR
# Reviewer #2 robustness: 4 mtDNA-CN exposures × 6 outcomes
# - Longchamps_2022 (CHARGE qPCR, UKB ~half of sample)
# - Chong_2022      (UKB array, full UKB)
# - Gupta_2023_raw  (UKB WGS, full UKB)
# - Gupta_2023_adj  (UKB WGS blood-cell adjusted, full UKB)
# Outcomes (5 primary + F3): SCZ, BD, MDD, ASD, ADHD, F3_NDev
# Date: 2026-05-26
#
# Reference: Mounier & Kutalik (2023) Bioinformatics 39:btac903
# Method: ld score regression intercept gives the correlation between Z-statistics
#         that share controls/cases; this lambda is used to deflate the MR estimate

suppressMessages({
  library(data.table)
  library(dplyr)
  library(MRlap)
})

PROJECT <- Sys.getenv("ITER003_HOME", unset = "/arf/scratch/ycicek/iter_003_mtdna_psy_mr")
setwd(PROJECT)

source(file.path(PROJECT, "scripts", "lib_robust.R"))

log_step("MRlap SAMPLE-OVERLAP CORRECTION START")
log_info(sprintf("Project: %s", PROJECT))
log_info(sprintf("R version: %s", R.version.string))

# ============================================================================
# Configuration: LD reference + HM3
# ============================================================================
# MRlap expects:
#   ld:  directory of <chr>.l2.ldscore.gz + <chr>.l2.M + <chr>.l2.M_5_50 (per-chromosome)
#   hm3: w_hm3.snplist (canonical 1.2M HM3 SNPs)
LD_DIR   <- "data/ldsc_mrlap/eur_w_ld_chr"
HM3_FILE <- "/arf/home/ycicek/EVOSCZ/data/ldsc/sldsc_ref/w_hm3.snplist"

stopifnot(dir.exists(LD_DIR))
stopifnot(file.exists(HM3_FILE))

log_info(sprintf("LD ref dir: %s", LD_DIR))
ld_files <- list.files(LD_DIR, pattern = "\\.l2\\.(ldscore\\.gz|M|M_5_50)$")
log_info(sprintf("LD ref files: %d (expect 66 = 22 chr x 3 files)", length(ld_files)))
log_info(sprintf("HM3 snplist: %s (%d MB)", HM3_FILE, round(file.info(HM3_FILE)$size / 1e6)))

# ============================================================================
# Output directory + checkpoint dir
# ============================================================================
dir.create("results/mrlap", showWarnings = FALSE, recursive = TRUE)
ckpt_dir <- "results/checkpoints/mrlap"
dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Exposures (harmonized): 4 mtDNA-CN GWAS
# ============================================================================
EXPOSURES <- list(
  list(name = "Longchamps_2022", path = "data/exposures_harmonized/Longchamps2022.tsv",     N = 465809),
  list(name = "Chong_2022",      path = "data/exposures_harmonized/Chong2022.tsv",          N = 383476),
  list(name = "Gupta_2023_raw",  path = "data/exposures_harmonized/Gupta2023_raw.tsv",      N = 155998),
  list(name = "Gupta_2023_adj",  path = "data/exposures_harmonized/Gupta2023_adjusted.tsv", N = 155998)
)

# ============================================================================
# Outcomes: 5 primary + F3
# ============================================================================
OUTCOMES <- list(
  list(name = "SCZ",     path = "data/outcomes/SCZ_PGC3_UKBdedup.gz",      n_cases = 67323,  n_controls = 93456),
  list(name = "BD",      path = "data/outcomes/BD_Mullins2024_EUR.gz",     n_cases = 59287,  n_controls = 781022),
  list(name = "MDD",     path = "data/outcomes/MDD_Adams2025_EUR.gz",      n_cases = 412305, n_controls = 1588397),
  list(name = "ADHD",    path = "data/outcomes/ADHD_Demontis2023.meta.gz", n_cases = 38691,  n_controls = 186843),
  list(name = "ASD",     path = "data/outcomes/ASD_Grove2019.gz",          n_cases = 18381,  n_controls = 27969),
  list(name = "F3_NDev", path = "data/outcomes/CDG3_F3_Neurodev.tsv.gz",   n_total = 84760)
)

# ============================================================================
# Helper: read outcome file with format auto-detection
# ============================================================================
read_outcome_file <- function(path) {
  con <- gzfile(path, "r")
  n_skip <- 0
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "##")) { n_skip <- n_skip + 1; next }
    break
  }
  close(con)
  if (n_skip > 0) log_info(sprintf("  skipping %d ## comment lines", n_skip))
  dt <- fread(path, skip = n_skip, fill = TRUE)
  if (startsWith(colnames(dt)[1], "#")) {
    setnames(dt, 1, sub("^#", "", colnames(dt)[1]))
  }
  dt
}

.detect_col <- function(cn, exact_set, regex = NULL) {
  cu <- toupper(cn)
  hit <- match(TRUE, cu %in% toupper(exact_set))
  if (!is.na(hit)) return(cn[hit])
  if (!is.null(regex)) {
    hit <- which(grepl(regex, cu, ignore.case = TRUE))[1]
    if (!is.na(hit)) return(cn[hit])
  }
  NA_character_
}

# ============================================================================
# Build a MRlap-format data.frame from outcome path
# Required columns (per MRlap docs): rsid (or snpid), chr (optional), pos (optional),
#   a1 (effect_allele), a2 (other_allele), beta, se, N
# We provide: snpid, chr, pos, a1, a2, beta, se, N (+ pval is optional)
# ============================================================================
build_mrlap_outcome_df <- function(out_def) {
  log_info(sprintf("[outcome %s] reading %s", out_def$name, out_def$path))
  dt <- read_outcome_file(out_def$path)
  cn <- colnames(dt)

  snp_col   <- .detect_col(cn, c("SNP","RSID","RS_ID","VARIANT_ID","MARKERNAME","ID","SNPID"))
  chr_col   <- .detect_col(cn, c("CHR","CHROM","CHROMOSOME","#CHR"))
  bp_col    <- .detect_col(cn, c("BP","POS","POSITION","BASE_PAIR_LOCATION"))
  a1_col    <- .detect_col(cn, c("A1","EFFECT_ALLELE","EA","ALT","ALLELE1","TESTED_ALLELE"))
  a2_col    <- .detect_col(cn, c("A2","OTHER_ALLELE","NEA","REF","ALLELE0"))
  beta_col  <- .detect_col(cn, c("BETA","B","EFFECT"))
  or_col    <- .detect_col(cn, c("OR","ODDS_RATIO"))
  se_col    <- .detect_col(cn, c("SE","STANDARD_ERROR","STDERR","STDDEV"))
  z_col     <- .detect_col(cn, c("Z","ZSCORE","Z_SCORE","Z_STAT","ZSTAT"))
  p_col     <- .detect_col(cn, c("P","PVAL","PVALUE","P_VALUE"), regex = "^P[\\-_]?VALUE$|^P[\\-_]VAL$")
  n_col     <- .detect_col(cn, c("N","N_TOTAL","SAMPLE_SIZE","NEFF","N_EFF","NEFFDIV2","TOTALSAMPLESIZE","NEFF_HALF"))
  eaf_col   <- .detect_col(cn, c("EAF","FRQ","FREQ","MAF","A1FREQ","EFFECT_ALLELE_FREQUENCY",
                                  "FREQ_TESTED_ALLELE_IN_HRS","FREQ1","FREQ_A1"),
                            regex = "^FRQ_A_\\d+$|^FREQ_A_\\d+$")

  # Z-score format detection
  can_derive_z <- (is.na(beta_col) && is.na(or_col) && is.na(se_col) &&
                   !is.na(z_col) && !is.na(eaf_col) && !is.na(n_col))

  # Sanity: must have SNP, A1, A2, and either (BETA+SE), or OR+SE (→ logOR), or Z-derive
  if (is.na(snp_col) || is.na(a1_col) || is.na(a2_col) ||
      (is.na(se_col) && !can_derive_z) ||
      (is.na(beta_col) && is.na(or_col) && !can_derive_z)) {
    log_warn(sprintf("[outcome %s] missing required column; cn=%s", out_def$name, paste(cn, collapse=",")))
    return(NULL)
  }

  # Z derivation (PTSD case — not used in current outcome list but kept for safety)
  if (can_derive_z) {
    eaf_vec <- as.numeric(dt[[eaf_col]])
    z_vec   <- as.numeric(dt[[z_col]])
    n_vec   <- as.numeric(dt[[n_col]])
    maf_vec <- pmin(eaf_vec, 1 - eaf_vec)
    dt[, SE_derived   := 1 / sqrt(2 * maf_vec * (1 - maf_vec) * n_vec)]
    dt[, BETA_derived := z_vec * SE_derived]
    se_col   <- "SE_derived"
    beta_col <- "BETA_derived"
  }

  # OR → log(OR)
  beta_col_to_use <- beta_col
  if (is.na(beta_col) && !is.na(or_col)) {
    dt[, BETA_log := log(suppressWarnings(as.numeric(get(or_col))))]
    beta_col_to_use <- "BETA_log"
  }

  # N resolution
  n_col_to_use <- n_col
  if (is.na(n_col)) {
    if (all(c("Nca", "Nco") %in% cn)) {
      dt[, N_total := Nca + Nco]
      n_col_to_use <- "N_total"
    } else {
      n_default <- if (!is.null(out_def$n_total)) out_def$n_total else (out_def$n_cases + out_def$n_controls)
      dt[, N_total := n_default]
      n_col_to_use <- "N_total"
    }
  }

  # Numeric coercion
  for (col in c(beta_col_to_use, se_col, p_col, n_col_to_use)) {
    if (!is.na(col) && col %in% colnames(dt) && !is.numeric(dt[[col]])) {
      dt[, (col) := suppressWarnings(as.numeric(get(col)))]
    }
  }
  # Drop NA
  dt <- dt[!is.na(get(beta_col_to_use)) & !is.na(get(se_col)) & get(se_col) > 0 &
            !is.na(get(n_col_to_use))]
  if (nrow(dt) == 0) {
    log_warn(sprintf("[outcome %s] zero rows after NA filtering", out_def$name))
    return(NULL)
  }

  # MRlap-required minimal cols: rsid, chr, pos, a1, a2, beta, se, N
  out_df <- data.frame(
    rsid = as.character(dt[[snp_col]]),
    chr  = if (!is.na(chr_col)) as.integer(dt[[chr_col]]) else NA_integer_,
    pos  = if (!is.na(bp_col))  as.integer(dt[[bp_col]])  else NA_integer_,
    a1   = toupper(as.character(dt[[a1_col]])),
    a2   = toupper(as.character(dt[[a2_col]])),
    beta = as.numeric(dt[[beta_col_to_use]]),
    se   = as.numeric(dt[[se_col]]),
    N    = as.numeric(dt[[n_col_to_use]]),
    stringsAsFactors = FALSE
  )
  log_info(sprintf("[outcome %s] %d SNPs ready for MRlap", out_def$name, nrow(out_df)))
  out_df
}

# ============================================================================
# Build MRlap-format exposure df from harmonized TSV
# ============================================================================
build_mrlap_exposure_df <- function(exp_def) {
  log_info(sprintf("[exposure %s] reading %s", exp_def$name, exp_def$path))
  dt <- fread(exp_def$path)
  # Standard harmonized cols: SNP, CHR, BP, A1, A2, BETA, SE, P, EAF, N
  # N may be missing in harmonized exposures — fallback to exp_def$N
  if (!"N" %in% colnames(dt)) {
    dt[, N := exp_def$N]
  } else if (any(is.na(dt$N))) {
    dt[is.na(N), N := exp_def$N]
  }
  out_df <- data.frame(
    rsid = as.character(dt$SNP),
    chr  = as.integer(dt$CHR),
    pos  = as.integer(dt$BP),
    a1   = toupper(as.character(dt$A1)),
    a2   = toupper(as.character(dt$A2)),
    beta = as.numeric(dt$BETA),
    se   = as.numeric(dt$SE),
    N    = as.numeric(dt$N),
    stringsAsFactors = FALSE
  )
  # Drop NA
  out_df <- out_df[complete.cases(out_df[, c("rsid","a1","a2","beta","se","N")]) & out_df$se > 0, ]
  log_info(sprintf("[exposure %s] %d SNPs ready for MRlap", exp_def$name, nrow(out_df)))
  out_df
}

# ============================================================================
# MRlap per (exposure, outcome) pair
# ============================================================================
run_mrlap_pair <- function(exp_def, out_def, exp_df, out_df) {
  pair_id <- sprintf("%s__%s", exp_def$name, out_def$name)
  ckpt <- file.path(ckpt_dir, sprintf("mrlap_%s.rds", pair_id))
  if (file.exists(ckpt)) {
    log_info(sprintf("[MRlap %s] checkpoint exists, skip", pair_id))
    return(invisible(readRDS(ckpt)))
  }
  log_info(sprintf("[MRlap %s] starting MRlap()", pair_id))
  t0 <- Sys.time()

  res <- tryCatch({
    r <- MRlap::MRlap(
      exposure        = exp_df,
      exposure_name   = exp_def$name,
      outcome         = out_df,
      outcome_name    = out_def$name,
      ld              = LD_DIR,
      hm3             = HM3_FILE,
      do_pruning      = TRUE,
      MR_threshold    = 5e-08,
      MR_pruning_dist = 500,
      MR_pruning_LD   = 0,
      MR_reverse      = 0.001,
      save_logfiles   = FALSE,
      verbose         = FALSE
    )
    list(pair_id = pair_id, exposure = exp_def$name, outcome = out_def$name,
         result = r, status = "OK",
         runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }, error = function(e) {
    log_error(sprintf("[MRlap %s] FAILED: %s", pair_id, conditionMessage(e)))
    list(pair_id = pair_id, exposure = exp_def$name, outcome = out_def$name,
         result = NULL, status = "FAILED", error = conditionMessage(e),
         runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  })

  saveRDS(res, ckpt)
  log_info(sprintf("[MRlap %s] done (%.1fs) status=%s",
                   pair_id, res$runtime_sec, res$status))
  invisible(res)
}

# ============================================================================
# Run grid: 4 exposures × 6 outcomes = 24 pairs
# ============================================================================
log_step("RUN GRID: 24 PAIRS")
all_results <- list()

# Cache outcome df reads — one per outcome (re-used across 4 exposures)
out_df_cache <- list()
for (out_def in OUTCOMES) {
  log_info(sprintf("[cache outcome] building df for %s", out_def$name))
  out_df_cache[[out_def$name]] <- build_mrlap_outcome_df(out_def)
}

for (exp_def in EXPOSURES) {
  log_step(sprintf("Exposure: %s", exp_def$name))
  exp_df <- build_mrlap_exposure_df(exp_def)
  for (out_def in OUTCOMES) {
    out_df <- out_df_cache[[out_def$name]]
    if (is.null(out_df)) {
      log_warn(sprintf("[%s -> %s] outcome df NULL, skip", exp_def$name, out_def$name))
      next
    }
    res <- run_mrlap_pair(exp_def, out_def, exp_df, out_df)
    all_results[[res$pair_id]] <- res
  }
  log_mem(sprintf("after exposure %s", exp_def$name))
}

saveRDS(all_results, "results/mrlap/all_mrlap_results.rds")

# ============================================================================
# Build summary TSV (defensive extraction — MRlap output structure)
# Expected output structure (per MRlap source):
#   $MRcorrection: list with components
#     - observed_effect, observed_effect_se, observed_effect_p
#     - corrected_effect, corrected_effect_se, corrected_effect_p
#     - test_difference, p_difference
#     - m_IVs (number of instruments used)
#   $LDSC: list with components
#     - h2_exp, h2_exp_se, h2_out, h2_out_se
#     - rg, rg_se, gencov_int (i.e., cross-trait intercept = lambda)
#     - h2_int_exp, h2_int_out (single-trait intercepts)
#   $GeneticArchitecture: list (h2 + intercept info)
# ============================================================================
log_step("BUILD SUMMARY TSV")

safe_get <- function(lst, key, default = NA_real_) {
  if (is.null(lst) || !key %in% names(lst)) return(default)
  v <- lst[[key]]
  if (is.null(v) || length(v) == 0) return(default)
  as.numeric(v[1])
}

summary_rows <- list()
for (key in names(all_results)) {
  r <- all_results[[key]]
  if (is.null(r) || r$status != "OK" || is.null(r$result)) {
    summary_rows[[key]] <- data.frame(
      pair = key, exposure = if (!is.null(r$exposure)) r$exposure else NA,
      outcome = if (!is.null(r$outcome)) r$outcome else NA,
      status = if (!is.null(r$status)) r$status else "MISSING",
      n_iv = NA_integer_,
      ivw_obs_b = NA_real_, ivw_obs_se = NA_real_, ivw_obs_p = NA_real_,
      ivw_corrected_b = NA_real_, ivw_corrected_se = NA_real_, ivw_corrected_p = NA_real_,
      test_difference = NA_real_, p_difference = NA_real_,
      h2_exp = NA_real_, h2_exp_se = NA_real_, h2_int_exp = NA_real_,
      h2_out = NA_real_, h2_out_se = NA_real_, h2_int_out = NA_real_,
      rg = NA_real_, rg_se = NA_real_, gencov_int = NA_real_,
      error = if (!is.null(r$error)) substr(r$error, 1, 200) else NA_character_,
      runtime_sec = if (!is.null(r$runtime_sec)) r$runtime_sec else NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }
  R <- r$result
  mrc  <- R$MRcorrection
  ldsc <- R$LDSC

  summary_rows[[key]] <- data.frame(
    pair = key, exposure = r$exposure, outcome = r$outcome,
    status = "OK",
    n_iv = as.integer(safe_get(mrc, "m_IVs", NA)),
    ivw_obs_b       = safe_get(mrc, "observed_effect"),
    ivw_obs_se      = safe_get(mrc, "observed_effect_se"),
    ivw_obs_p       = safe_get(mrc, "observed_effect_p"),
    ivw_corrected_b = safe_get(mrc, "corrected_effect"),
    ivw_corrected_se= safe_get(mrc, "corrected_effect_se"),
    ivw_corrected_p = safe_get(mrc, "corrected_effect_p"),
    test_difference = safe_get(mrc, "test_difference"),
    p_difference    = safe_get(mrc, "p_difference"),
    h2_exp          = safe_get(ldsc, "h2_exp"),
    h2_exp_se       = safe_get(ldsc, "h2_exp_se"),
    h2_int_exp      = safe_get(ldsc, "h2_int_exp"),
    h2_out          = safe_get(ldsc, "h2_out"),
    h2_out_se       = safe_get(ldsc, "h2_out_se"),
    h2_int_out      = safe_get(ldsc, "h2_int_out"),
    rg              = safe_get(ldsc, "rg"),
    rg_se           = safe_get(ldsc, "rg_se"),
    gencov_int      = safe_get(ldsc, "gencov_int"),
    error           = NA_character_,
    runtime_sec     = r$runtime_sec,
    stringsAsFactors = FALSE
  )
}
summary_df <- bind_rows(summary_rows)
fwrite(summary_df, "results/mrlap/MRLAP_summary.tsv", sep = "\t")
log_info(sprintf("MRlap summary saved: %d rows", nrow(summary_df)))

log_step("MRlap PIPELINE COMPLETE")
log_info(sprintf("Summary: %s/results/mrlap/MRLAP_summary.tsv", PROJECT))
