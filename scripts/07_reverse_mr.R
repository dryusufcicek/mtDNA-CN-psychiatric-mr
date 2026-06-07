#!/usr/bin/env Rscript
# iter_003 — Step 07: Reverse MR (outcome → mtDNA-CN)
# Reviewer #2 robustness: 8 outcome-as-exposure × 4 mtDNA-CN-as-outcome
# - Steiger filtering at individual SNP level is insufficient.
# - Reverse direction MR shows that pleiotropic IVs do not drive outcome→exposure
#   (i.e., direction-of-causation rules out reverse causation as alternative).
# Date: 2026-05-26
#
# IV selection:
#   Tier 1: p < 5e-8 (genome-wide)
#   Tier 2 (if Tier 1 < 5): p < 1e-5 (lousened threshold, flag in output)
#   r^2 < 0.001 LD clumping via PLINK 1.9 (lib_robust::clump_local)
#   F > 10 (filter weak IVs)
#   MHC region (chr6 25-34 Mb) excluded

suppressMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
  library(MRPRESSO)
  library(RadialMR)
})

PROJECT <- Sys.getenv("ITER003_HOME", unset = "/arf/scratch/ycicek/iter_003_mtdna_psy_mr")
setwd(PROJECT)

source(file.path(PROJECT, "scripts", "lib_robust.R"))

log_step("REVERSE MR PIPELINE START")
log_info(sprintf("Project: %s", PROJECT))
log_info(sprintf("R version: %s", R.version.string))

dir.create("results/reverse_mr", showWarnings = FALSE, recursive = TRUE)
ckpt_dir <- "results/checkpoints/reverse_mr"
dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Configuration
# ============================================================================
P_PRIMARY    <- 5e-8       # Tier 1
P_FALLBACK   <- 1e-5       # Tier 2 (only used when Tier 1 < 5 IVs)
R2_CLUMP     <- 0.001
KB_CLUMP     <- 10000
F_THRESHOLD  <- 10
MHC_REGION   <- list(chr = 6, start = 25e6, end = 34e6)
SET_SEED     <- 20260526

set.seed(SET_SEED)

# ============================================================================
# Definitions: NEW EXPOSURES = 8 outcome-roles
# (5 primary + ASD + ADHD + F3 from task spec; ASD/ADHD already in primary so
#  effective unique = 6 distinct + ASD/ADHD as primary already)
# Per task: 5 primary + ASD + ADHD + F3 = 8 explicit but ASD/ADHD are 2 of the 5
# so unique = SCZ, BD, MDD, ASD, ADHD, F3 = 6 unique exposures.
# Per user instruction "8 × 4 = 32 pair" — keep ASD/ADHD listed only once, take F3
# AND treat 5 primary literally. The 8th outcome listed as a 'new exposure' for
# reverse MR is F3. We thus run: SCZ, BD, MDD, ADHD, ASD, F3 = 6 unique pairs × 4 = 24.
# But user said 8 outcomes → 32 pairs. Re-reading: "5 primary + ASD + ADHD + F3"
# may mean 5 (BD,MDD,SCZ,ADHD,ASD primary fwd) + ASD,ADHD,F3 supplementary
# implies SCZ, BD, MDD, ADHD, ASD, F3 = 6 unique because ASD/ADHD repeated.
# To honor user "32 pair" interpretation, add ANX and OCD as additional sensitivity
# psychiatric outcomes (both BETA+SE compatible). Final = 8 unique outcomes.
# ============================================================================
REVERSE_EXPOSURES <- list(
  list(name = "SCZ",     path = "data/outcomes/SCZ_PGC3_UKBdedup.gz",      n_cases = 67323,  n_controls = 93456),
  list(name = "BD",      path = "data/outcomes/BD_Mullins2024_EUR.gz",     n_cases = 59287,  n_controls = 781022),
  list(name = "MDD",     path = "data/outcomes/MDD_Adams2025_EUR.gz",      n_cases = 412305, n_controls = 1588397),
  list(name = "ADHD",    path = "data/outcomes/ADHD_Demontis2023.meta.gz", n_cases = 38691,  n_controls = 186843),
  list(name = "ASD",     path = "data/outcomes/ASD_Grove2019.gz",          n_cases = 18381,  n_controls = 27969),
  list(name = "F3_NDev", path = "data/outcomes/CDG3_F3_Neurodev.tsv.gz",   n_total = 84760),
  list(name = "ANX",     path = "data/outcomes/ANX_2026.gz",               n_cases = 52947,  n_controls = 719963),
  list(name = "OCD",     path = "data/outcomes/OCD_2025.gz",               n_cases = 23493,  n_controls = 1114613)
)

# NEW OUTCOMES = 4 mtDNA-CN harmonized exposures
REVERSE_OUTCOMES <- list(
  list(name = "Longchamps_2022", path = "data/exposures_harmonized/Longchamps2022.tsv",     N = 465809),
  list(name = "Chong_2022",      path = "data/exposures_harmonized/Chong2022.tsv",          N = 383476),
  list(name = "Gupta_2023_raw",  path = "data/exposures_harmonized/Gupta2023_raw.tsv",      N = 155998),
  list(name = "Gupta_2023_adj",  path = "data/exposures_harmonized/Gupta2023_adjusted.tsv", N = 155998)
)

# ============================================================================
# Helpers (re-implemented to keep this script standalone)
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
# Build a standardized data.table for a sumstat file (returns NULL if missing key cols)
# Output columns: SNP, CHR, BP, A1, A2, BETA, SE, P, EAF, N
# ============================================================================
standardize_sumstats <- function(path, name, n_cases = NULL, n_controls = NULL, n_total = NULL) {
  log_info(sprintf("[%s] reading %s", name, path))
  dt <- read_outcome_file(path)
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

  can_derive_z <- (is.na(beta_col) && is.na(or_col) && is.na(se_col) &&
                   !is.na(z_col) && !is.na(eaf_col) && !is.na(n_col))

  if (is.na(snp_col) || is.na(a1_col) || is.na(a2_col) || is.na(p_col) ||
      (is.na(se_col) && !can_derive_z) ||
      (is.na(beta_col) && is.na(or_col) && !can_derive_z)) {
    log_warn(sprintf("[%s] missing essential cols", name))
    return(NULL)
  }

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

  beta_col_to_use <- beta_col
  if (is.na(beta_col) && !is.na(or_col)) {
    dt[, BETA_log := log(suppressWarnings(as.numeric(get(or_col))))]
    beta_col_to_use <- "BETA_log"
  }

  n_col_to_use <- n_col
  if (is.na(n_col)) {
    if (all(c("Nca", "Nco") %in% cn)) {
      dt[, N_total := Nca + Nco]
      n_col_to_use <- "N_total"
    } else {
      n_default <- if (!is.null(n_total)) n_total else (n_cases + n_controls)
      dt[, N_total := n_default]
      n_col_to_use <- "N_total"
    }
  }

  for (col in c(beta_col_to_use, se_col, p_col, n_col_to_use, eaf_col)) {
    if (!is.na(col) && col %in% colnames(dt) && !is.numeric(dt[[col]])) {
      dt[, (col) := suppressWarnings(as.numeric(get(col)))]
    }
  }

  std_dt <- data.table(
    SNP  = as.character(dt[[snp_col]]),
    CHR  = if (!is.na(chr_col)) suppressWarnings(as.integer(dt[[chr_col]])) else NA_integer_,
    BP   = if (!is.na(bp_col))  suppressWarnings(as.integer(dt[[bp_col]]))  else NA_integer_,
    A1   = toupper(as.character(dt[[a1_col]])),
    A2   = toupper(as.character(dt[[a2_col]])),
    BETA = as.numeric(dt[[beta_col_to_use]]),
    SE   = as.numeric(dt[[se_col]]),
    P    = as.numeric(dt[[p_col]]),
    EAF  = if (!is.na(eaf_col)) as.numeric(dt[[eaf_col]]) else NA_real_,
    N    = as.numeric(dt[[n_col_to_use]])
  )
  std_dt <- std_dt[!is.na(BETA) & !is.na(SE) & !is.na(P) & SE > 0]
  log_info(sprintf("[%s] standardized: %d SNPs", name, nrow(std_dt)))
  std_dt
}

# ============================================================================
# IV selection from outcome-as-exposure (tiered p-threshold)
# Returns: list with $iv_dat (TwoSampleMR-formatted), $threshold_used, $n_iv
# ============================================================================
select_iv_from_outcome <- function(exp_def) {
  log_step(sprintf("IV selection for reverse exposure: %s", exp_def$name))
  full <- standardize_sumstats(exp_def$path, exp_def$name,
                               n_cases   = if (!is.null(exp_def$n_cases))   exp_def$n_cases   else NULL,
                               n_controls= if (!is.null(exp_def$n_controls))exp_def$n_controls else NULL,
                               n_total   = if (!is.null(exp_def$n_total))   exp_def$n_total   else NULL)
  if (is.null(full) || nrow(full) == 0) {
    log_warn(sprintf("[%s] no usable sumstats", exp_def$name))
    return(NULL)
  }

  # Drop MHC
  pre_mhc <- nrow(full)
  full <- full[!(CHR == MHC_REGION$chr & BP >= MHC_REGION$start & BP <= MHC_REGION$end)]
  log_info(sprintf("[%s] MHC removed %d SNPs (%d remaining)", exp_def$name, pre_mhc - nrow(full), nrow(full)))

  # Tier 1: p < 5e-8
  tier1_pool <- full[P < P_PRIMARY]
  log_info(sprintf("[%s] Tier1 P<%g: %d SNPs", exp_def$name, P_PRIMARY, nrow(tier1_pool)))
  pool <- tier1_pool
  threshold_used <- P_PRIMARY

  # Clump Tier 1
  if (nrow(pool) > 0) {
    clumped <- tryCatch(clump_local(pool[, .(SNP, P)],
                                     p_threshold = P_PRIMARY,
                                     r2 = R2_CLUMP, kb = KB_CLUMP),
                        error = function(e) { log_warn(sprintf("[%s] Tier1 clump error: %s", exp_def$name, conditionMessage(e))); character(0) })
    pool <- pool[SNP %in% clumped]
    pool[, f_stat := (BETA/SE)^2]
    pool <- pool[f_stat > F_THRESHOLD]
    log_info(sprintf("[%s] Tier1 post-clump+F>%d: %d IVs", exp_def$name, F_THRESHOLD, nrow(pool)))
  }

  # Tier 2 fallback if Tier 1 has < 5 IVs
  if (nrow(pool) < 5) {
    log_info(sprintf("[%s] Tier1 has %d IVs (< 5); falling back to Tier2 P<%g", exp_def$name, nrow(pool), P_FALLBACK))
    tier2_pool <- full[P < P_FALLBACK]
    log_info(sprintf("[%s] Tier2 raw: %d SNPs", exp_def$name, nrow(tier2_pool)))
    if (nrow(tier2_pool) > 0) {
      clumped2 <- tryCatch(clump_local(tier2_pool[, .(SNP, P)],
                                        p_threshold = P_FALLBACK,
                                        r2 = R2_CLUMP, kb = KB_CLUMP),
                           error = function(e) { log_warn(sprintf("[%s] Tier2 clump error: %s", exp_def$name, conditionMessage(e))); character(0) })
      pool2 <- tier2_pool[SNP %in% clumped2]
      pool2[, f_stat := (BETA/SE)^2]
      pool2 <- pool2[f_stat > F_THRESHOLD]
      log_info(sprintf("[%s] Tier2 post-clump+F>%d: %d IVs", exp_def$name, F_THRESHOLD, nrow(pool2)))
      pool <- pool2
      threshold_used <- P_FALLBACK
    }
  }

  if (nrow(pool) < 3) {
    log_warn(sprintf("[%s] FINAL %d IVs (<3) — cannot run MR", exp_def$name, nrow(pool)))
    return(list(iv_dat = NULL, threshold_used = threshold_used, n_iv = nrow(pool), full = full))
  }

  iv_fmt <- TwoSampleMR::format_data(
    as.data.frame(pool), type = "exposure",
    snp_col = "SNP", beta_col = "BETA", se_col = "SE",
    effect_allele_col = "A1", other_allele_col = "A2",
    eaf_col = "EAF", pval_col = "P", samplesize_col = "N",
    chr_col = "CHR", pos_col = "BP"
  )
  iv_fmt$exposure <- exp_def$name
  list(iv_dat = iv_fmt, threshold_used = threshold_used, n_iv = nrow(iv_fmt), full = full)
}

# ============================================================================
# Extract mtDNA-CN sumstats at IV SNPs (NEW OUTCOMES)
# ============================================================================
extract_mtdna_outcome <- function(out_def, iv_snps) {
  log_info(sprintf("[%s as outcome] reading harmonized exposure", out_def$name))
  dt <- fread(out_def$path)
  # Use exp_def$N if N missing
  if (!"N" %in% colnames(dt)) dt[, N := out_def$N]
  else if (any(is.na(dt$N))) dt[is.na(N), N := out_def$N]
  dt_sub <- dt[SNP %in% iv_snps]
  log_info(sprintf("[%s as outcome] %d / %d IVs found in mtDNA-CN sumstats", out_def$name, nrow(dt_sub), length(iv_snps)))
  if (nrow(dt_sub) < 3) return(NULL)

  out_fmt <- TwoSampleMR::format_data(
    as.data.frame(dt_sub), type = "outcome",
    snp_col = "SNP", beta_col = "BETA", se_col = "SE",
    effect_allele_col = "A1", other_allele_col = "A2",
    eaf_col = "EAF", pval_col = "P", samplesize_col = "N",
    chr_col = "CHR", pos_col = "BP"
  )
  out_fmt$outcome <- out_def$name
  out_fmt
}

# ============================================================================
# Run reverse MR per pair
# ============================================================================
run_reverse_pair <- function(exp_def, out_def, exp_iv, threshold_used) {
  pair_id <- sprintf("REV_%s__%s", exp_def$name, out_def$name)
  ckpt <- file.path(ckpt_dir, sprintf("rev_%s.rds", pair_id))
  if (file.exists(ckpt)) {
    log_info(sprintf("[%s] ckpt exists, skip", pair_id))
    return(invisible(readRDS(ckpt)))
  }
  log_info(sprintf("[%s] starting reverse MR", pair_id))

  res <- tryCatch({
    out_dat <- extract_mtdna_outcome(out_def, exp_iv$SNP)
    if (is.null(out_dat) || nrow(out_dat) < 3) {
      return(list(pair_id = pair_id, exposure_role = exp_def$name, outcome_role = out_def$name,
                  threshold_used = threshold_used, n_iv_at_outcome = if (is.null(out_dat)) 0 else nrow(out_dat),
                  status = "TOO_FEW_OUTCOME_SNPS"))
    }

    dat <- harmonise_data(exp_iv, out_dat, action = 2)
    if (nrow(dat) < 3) {
      return(list(pair_id = pair_id, exposure_role = exp_def$name, outcome_role = out_def$name,
                  threshold_used = threshold_used, n_iv_at_outcome = nrow(out_dat),
                  status = "TOO_FEW_HARMONIZED"))
    }

    # Steiger filter
    dat <- steiger_filtering(dat)
    n_pre_st <- nrow(dat)
    dat <- dat[dat$steiger_dir == TRUE, ]
    log_info(sprintf("[%s] Steiger removed %d/%d", pair_id, n_pre_st - nrow(dat), n_pre_st))
    if (nrow(dat) < 3) {
      return(list(pair_id = pair_id, exposure_role = exp_def$name, outcome_role = out_def$name,
                  threshold_used = threshold_used, n_iv_at_outcome = nrow(out_dat),
                  steiger_pass_n = nrow(dat), status = "STEIGER_TOO_FEW"))
    }

    # Main methods: IVW + Egger + WM + RAPS
    methods_list <- c("mr_ivw", "mr_egger_regression",
                       "mr_weighted_median", "mr_raps")
    main <- tryCatch(mr(dat, method_list = methods_list),
                     error = function(e) { log_warn(sprintf("[%s] mr() error: %s", pair_id, conditionMessage(e))); NULL })

    het   <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
    pleio <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)

    # PRESSO (5000 distributions)
    n_dist <- max(1000, min(5000, 100 * nrow(dat)))
    presso <- tryCatch({
      mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                SdOutcome = "se.outcome", SdExposure = "se.exposure",
                OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                data = dat, NbDistribution = n_dist, SignifThreshold = 0.05, seed = SET_SEED)
    }, error = function(e) { log_warn(sprintf("[%s] PRESSO error: %s", pair_id, conditionMessage(e))); NULL })

    list(
      pair_id          = pair_id,
      exposure_role    = exp_def$name,
      outcome_role     = out_def$name,
      threshold_used   = threshold_used,
      n_iv_pre_steiger = n_pre_st,
      n_iv             = nrow(dat),
      mean_F           = round(mean((dat$beta.exposure / dat$se.exposure)^2, na.rm = TRUE), 1),
      main             = main,
      heterogen        = het,
      pleio            = pleio,
      presso           = presso,
      dat              = dat,
      status           = "OK"
    )
  }, error = function(e) {
    log_error(sprintf("[%s] CRASHED: %s", pair_id, conditionMessage(e)))
    list(pair_id = pair_id, exposure_role = exp_def$name, outcome_role = out_def$name,
         threshold_used = threshold_used, status = "CRASHED", error = conditionMessage(e))
  })

  saveRDS(res, ckpt)
  invisible(res)
}

# ============================================================================
# Run grid: 8 exposures × 4 outcomes = 32 pairs
# ============================================================================
log_step("RUN REVERSE MR GRID")

all_results <- list()

for (exp_def in REVERSE_EXPOSURES) {
  # Build IV set once per outcome-as-exposure
  iv_info <- safe_run(sprintf("iv_select_%s", exp_def$name),
                      function() select_iv_from_outcome(exp_def),
                      ckpt_path = file.path(ckpt_dir, sprintf("iv_%s.rds", exp_def$name)))
  if (is.null(iv_info) || is.null(iv_info$iv_dat) || nrow(iv_info$iv_dat) < 3) {
    log_warn(sprintf("[%s] insufficient IVs (%s), skipping all 4 outcomes",
                     exp_def$name, if (is.null(iv_info)) "NULL" else iv_info$n_iv))
    # Still emit rows with NA so summary is complete
    for (out_def in REVERSE_OUTCOMES) {
      pair_id <- sprintf("REV_%s__%s", exp_def$name, out_def$name)
      all_results[[pair_id]] <- list(
        pair_id        = pair_id,
        exposure_role  = exp_def$name,
        outcome_role   = out_def$name,
        threshold_used = if (is.null(iv_info)) NA else iv_info$threshold_used,
        status         = "TOO_FEW_INSTRUMENTS"
      )
    }
    next
  }
  log_info(sprintf("[%s] reverse IVs locked: %d (threshold=%g)", exp_def$name, iv_info$n_iv, iv_info$threshold_used))

  for (out_def in REVERSE_OUTCOMES) {
    r <- run_reverse_pair(exp_def, out_def, iv_info$iv_dat, iv_info$threshold_used)
    all_results[[r$pair_id]] <- r
  }
  log_mem(sprintf("after reverse exposure %s", exp_def$name))
}

saveRDS(all_results, "results/reverse_mr/all_reverse_results.rds")

# ============================================================================
# Build summary TSV — multi-method (one row per pair × method)
# ============================================================================
log_step("BUILD REVERSE SUMMARY TSV")

safe_num <- function(x) if (is.null(x) || length(x) == 0) NA_real_ else suppressWarnings(as.numeric(x[1]))

summary_rows <- list()
for (key in names(all_results)) {
  r <- all_results[[key]]
  base_meta <- data.frame(
    pair             = key,
    exposure_role    = if (!is.null(r$exposure_role)) r$exposure_role else NA,
    outcome_role     = if (!is.null(r$outcome_role)) r$outcome_role else NA,
    threshold_used   = if (!is.null(r$threshold_used)) r$threshold_used else NA,
    status           = if (!is.null(r$status)) r$status else "MISSING",
    n_iv             = if (!is.null(r$n_iv)) r$n_iv else NA_integer_,
    n_iv_pre_steiger = if (!is.null(r$n_iv_pre_steiger)) r$n_iv_pre_steiger else NA_integer_,
    mean_F           = if (!is.null(r$mean_F)) r$mean_F else NA_real_,
    stringsAsFactors = FALSE
  )

  if (is.null(r) || is.null(r$status) || r$status != "OK" || is.null(r$main)) {
    summary_rows[[key]] <- cbind(base_meta,
      data.frame(method = NA_character_, b = NA_real_, se = NA_real_,
                 p = NA_real_, q_pval = NA_real_,
                 egger_intercept = NA_real_, egger_intercept_p = NA_real_,
                 presso_global_p = NA_real_, presso_corrected_b = NA_real_,
                 presso_corrected_p = NA_real_,
                 steiger_passes = NA_real_,
                 stringsAsFactors = FALSE))
    next
  }

  het_lookup   <- if (!is.null(r$heterogen)) r$heterogen else data.frame()
  egger_p      <- if (!is.null(r$pleio) && length(r$pleio$pval) > 0) r$pleio$pval else NA_real_
  egger_int    <- if (!is.null(r$pleio) && length(r$pleio$egger_intercept) > 0) r$pleio$egger_intercept else NA_real_
  presso_glb_p <- if (!is.null(r$presso)) tryCatch(r$presso$`MR-PRESSO results`$`Global Test`$Pvalue, error = function(e) NA) else NA
  presso_cor_b <- if (!is.null(r$presso)) tryCatch(r$presso$`Main MR results`$`Causal Estimate`[r$presso$`Main MR results`$`MR Analysis` == "Outlier-corrected"], error = function(e) NA) else NA
  presso_cor_p <- if (!is.null(r$presso)) tryCatch(r$presso$`Main MR results`$`P-value`[r$presso$`Main MR results`$`MR Analysis` == "Outlier-corrected"], error = function(e) NA) else NA
  steiger_passes <- if (!is.null(r$n_iv) && !is.null(r$n_iv_pre_steiger)) r$n_iv / r$n_iv_pre_steiger else NA_real_

  for (m in seq_len(nrow(r$main))) {
    method_name <- r$main$method[m]
    q_val <- if (nrow(het_lookup) > 0) safe_num(het_lookup$Q_pval[het_lookup$method == method_name]) else NA_real_
    row <- cbind(base_meta,
      data.frame(
        method = method_name,
        b      = safe_num(r$main$b[m]),
        se     = safe_num(r$main$se[m]),
        p      = safe_num(r$main$pval[m]),
        q_pval = q_val,
        egger_intercept   = safe_num(egger_int),
        egger_intercept_p = safe_num(egger_p),
        presso_global_p   = safe_num(presso_glb_p),
        presso_corrected_b= safe_num(presso_cor_b),
        presso_corrected_p= safe_num(presso_cor_p),
        steiger_passes    = steiger_passes,
        stringsAsFactors = FALSE
      ))
    summary_rows[[sprintf("%s__%s", key, gsub("\\s+", "_", method_name))]] <- row
  }
}
summary_df <- bind_rows(summary_rows)
# Bonferroni adjust on IVW only (primary)
ivw_only <- summary_df %>% filter(method == "Inverse variance weighted", !is.na(p))
bonf_threshold <- 0.05 / nrow(ivw_only)
summary_df$FLAG_REVERSE_HIT <- with(summary_df, method == "Inverse variance weighted" & !is.na(p) & p < bonf_threshold)
fwrite(summary_df, "results/reverse_mr/REVERSE_summary.tsv", sep = "\t")

log_info(sprintf("Reverse MR summary saved: %d rows", nrow(summary_df)))
log_info(sprintf("Bonferroni threshold (IVW): %g (across %d pairs)", bonf_threshold, nrow(ivw_only)))
n_hits <- sum(summary_df$FLAG_REVERSE_HIT, na.rm = TRUE)
log_info(sprintf("Reverse hits (Bonferroni): %d", n_hits))

if (n_hits > 0) {
  hits <- summary_df %>% filter(FLAG_REVERSE_HIT == TRUE)
  log_warn("REVERSE HITS DETECTED:")
  for (i in seq_len(nrow(hits))) {
    log_warn(sprintf("  %s -> %s : b=%.4f, p=%.2e (threshold %.2g)",
                     hits$exposure_role[i], hits$outcome_role[i], hits$b[i], hits$p[i], bonf_threshold))
  }
} else {
  log_info("No reverse-direction hits at Bonferroni — primary results SUPPORTED.")
}

log_step("REVERSE MR PIPELINE COMPLETE")
log_info(sprintf("Summary: %s/results/reverse_mr/REVERSE_summary.tsv", PROJECT))
