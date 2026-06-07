# iter_003 — Robustness library
# Shared helpers for defensive execution + checkpointing

suppressMessages({
  library(data.table)
  library(dplyr)
})

# ============================================================================
# Logging
# ============================================================================
log_msg <- function(level, msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] [%s] %s\n", ts, level, msg))
  flush.console()
}
log_info  <- function(msg) log_msg("INFO ", msg)
log_warn  <- function(msg) log_msg("WARN ", msg)
log_error <- function(msg) log_msg("ERROR", msg)
log_step  <- function(msg) {
  log_msg("STEP ", paste0(strrep("=", 8), " ", msg, " ", strrep("=", 8)))
}

# ============================================================================
# Safe execution with checkpoint
# ============================================================================
safe_run <- function(label, fn, ckpt_path = NULL, force = FALSE) {
  if (!is.null(ckpt_path) && file.exists(ckpt_path) && !force) {
    log_info(sprintf("[%s] skip: checkpoint exists at %s", label, ckpt_path))
    return(invisible(readRDS(ckpt_path)))
  }
  log_info(sprintf("[%s] starting", label))
  t0 <- Sys.time()
  result <- tryCatch(fn(), error = function(e) {
    log_error(sprintf("[%s] FAILED: %s", label, conditionMessage(e)))
    NULL
  })
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(result)) {
    log_error(sprintf("[%s] returned NULL (after %.1fs)", label, dt))
  } else {
    log_info(sprintf("[%s] OK (%.1fs)", label, dt))
    if (!is.null(ckpt_path)) {
      dir.create(dirname(ckpt_path), showWarnings = FALSE, recursive = TRUE)
      saveRDS(result, ckpt_path)
      log_info(sprintf("[%s] checkpoint saved: %s", label, ckpt_path))
    }
  }
  invisible(result)
}

# ============================================================================
# Local LD clumping (uses PLINK + 1KG reference)
# ============================================================================
PLINK_BIN <- "/arf/home/ycicek/bin/plink"
LD_REF_DIR <- "/arf/home/ycicek/EVOSCZ/data/ldsc/sldsc_ref/1000G_EUR_Phase3_plink"
LD_REF_PREFIX <- file.path(LD_REF_DIR, "1000G.EUR.QC")  # per-chr files

# Concatenate all per-chr files into one merged set if not exists
ensure_ld_merged <- function() {
  merged_prefix <- file.path(LD_REF_DIR, "1000G.EUR.QC.merged")
  if (file.exists(paste0(merged_prefix, ".bed"))) return(merged_prefix)
  log_info("Merging per-chr 1KG EUR PLINK files into one set...")
  mergelist <- tempfile()
  writeLines(paste0(LD_REF_PREFIX, ".", 2:22), mergelist)
  cmd <- sprintf("%s --bfile %s.1 --merge-list %s --make-bed --out %s",
                 PLINK_BIN, LD_REF_PREFIX, mergelist, merged_prefix)
  log_info(paste("Running:", cmd))
  system(cmd, intern = FALSE)
  merged_prefix
}

# Local clumping using PLINK
clump_local <- function(snp_df, p_threshold = 5e-8, r2 = 0.001, kb = 10000) {
  # snp_df should have columns: SNP, P
  stopifnot(all(c("SNP", "P") %in% colnames(snp_df)))
  ref_prefix <- ensure_ld_merged()

  tmp <- tempfile(fileext = ".tsv")
  fwrite(snp_df %>% select(SNP, P), tmp, sep = "\t")
  out_prefix <- tempfile()
  cmd <- sprintf("%s --bfile %s --clump %s --clump-p1 %g --clump-p2 %g --clump-r2 %g --clump-kb %d --out %s 2>&1",
                 PLINK_BIN, ref_prefix, tmp, p_threshold, p_threshold, r2, kb, out_prefix)
  log_info(sprintf("Clumping %d SNPs (P<%g)...", sum(snp_df$P < p_threshold, na.rm = TRUE), p_threshold))
  system(cmd, intern = FALSE)
  clump_file <- paste0(out_prefix, ".clumped")
  if (!file.exists(clump_file)) {
    log_warn("PLINK clumping produced no .clumped file (all SNPs may be in LD with sentinels)")
    return(character(0))
  }
  cl <- fread(clump_file, header = TRUE, fill = TRUE)
  unlink(c(tmp, paste0(out_prefix, ".*")))
  cl$SNP
}

# ============================================================================
# Memory monitoring (simple)
# ============================================================================
mem_used_gb <- function() {
  round(sum(gc()[,2]) * 1e-3, 2)
}
log_mem <- function(label) {
  log_info(sprintf("[%s] mem ~%.1f GB used", label, mem_used_gb()))
}

# ============================================================================
# Safe MR pair execution
# ============================================================================
safe_mr_pair <- function(exp_dat, out_dat, exp_name, out_name, ckpt_dir, force = FALSE) {
  pair_id <- sprintf("%s__%s", exp_name, out_name)
  ckpt <- file.path(ckpt_dir, sprintf("mr_%s.rds", pair_id))
  if (file.exists(ckpt) && !force) {
    log_info(sprintf("[MR %s] checkpoint exists, skip", pair_id))
    return(invisible(readRDS(ckpt)))
  }

  res <- tryCatch({
    library(TwoSampleMR)
    # Harmonize
    dat <- harmonise_data(exp_dat, out_dat, action = 2)
    if (nrow(dat) < 3) {
      log_warn(sprintf("[MR %s] <3 SNPs after harmonization (%d)", pair_id, nrow(dat)))
      return(list(pair_id = pair_id, n_ivs = nrow(dat), main = NULL, status = "TOO_FEW_IVS"))
    }
    # Steiger
    dat <- steiger_filtering(dat)
    dat <- dat[dat$steiger_dir == TRUE, ]
    if (nrow(dat) < 3) {
      log_warn(sprintf("[MR %s] <3 SNPs after Steiger (%d)", pair_id, nrow(dat)))
      return(list(pair_id = pair_id, n_ivs = nrow(dat), main = NULL, status = "STEIGER_TOO_FEW"))
    }

    methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_raps")
    main <- mr(dat, method_list = methods)

    het <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
    pleio <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)

    # PRESSO with adaptive distribution count
    n_dist <- max(1000, min(10000, 100 * nrow(dat)))
    presso <- tryCatch({
      suppressMessages(library(MRPRESSO))
      mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                SdOutcome = "se.outcome", SdExposure = "se.exposure",
                OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat,
                NbDistribution = n_dist, SignifThreshold = 0.05, seed = 20260522)
    }, error = function(e) { log_warn(sprintf("[MR %s] PRESSO failed: %s", pair_id, conditionMessage(e))); NULL })

    list(
      pair_id = pair_id,
      exposure = exp_name,
      outcome  = out_name,
      n_ivs    = nrow(dat),
      main     = main,
      heterogen= het,
      pleio    = pleio,
      presso   = presso,
      dat      = dat,
      status   = "OK"
    )
  }, error = function(e) {
    log_error(sprintf("[MR %s] CRASHED: %s", pair_id, conditionMessage(e)))
    list(pair_id = pair_id, status = "CRASHED", error = conditionMessage(e))
  })

  dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(res, ckpt)
  invisible(res)
}
