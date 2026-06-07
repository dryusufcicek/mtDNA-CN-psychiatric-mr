#!/usr/bin/env Rscript
# iter_003 — Step 03: Forward two-sample MR (mtDNA-CN → psychiatric outcomes)
# Per exposure × outcome: IVW + Egger + WM + RAPS + MR-PRESSO + sensitivity
# Date: 2026-05-22

suppressMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)        # devtools::install_github("MRCIEU/TwoSampleMR")
  library(MendelianRandomization)
  library(MRPRESSO)           # devtools::install_github("rondolab/MR-PRESSO")
  library(RadialMR)
  library(ggplot2)
})

PROJECT <- Sys.getenv("ITER003_HOME", unset = getwd())
setwd(PROJECT)
dir.create("results/forward_mr", showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Configuration
# ============================================================================
P_THRESHOLD     <- 5e-8
R2_CLUMP        <- 0.001
KB_CLUMP        <- 10000
F_THRESHOLD     <- 10
STEIGER_FILTER  <- TRUE
MHC_REGION      <- list(chr = 6, start = 25e6, end = 34e6)
SET_SEED        <- 20260522

set.seed(SET_SEED)

# Exposure list (4 harmonized mtDNA-CN GWAS)
EXPOSURES <- list(
  Longchamps_CHARGE_qPCR     = "data/exposures_harmonized/Longchamps2022.tsv",
  Chong_UKB_array            = "data/exposures_harmonized/Chong2022.tsv",
  Gupta_UKB_WGS_adjusted     = "data/exposures_harmonized/Gupta2023_adjusted.tsv",
  Gupta_UKB_WGS_raw          = "data/exposures_harmonized/Gupta2023_raw.tsv"
)

# Outcome list (5 primary + 3 bonus + 2 neg controls)
OUTCOMES <- list(
  SCZ      = "data/outcomes_harmonized/SCZ.tsv",
  BD       = "data/outcomes_harmonized/BD.tsv",
  MDD      = "data/outcomes_harmonized/MDD.tsv",
  ADHD     = "data/outcomes_harmonized/ADHD.tsv",
  ASD      = "data/outcomes_harmonized/ASD.tsv",
  PFactor  = "data/outcomes_harmonized/CDG3_PFactor.tsv",
  F3_NDev  = "data/outcomes_harmonized/CDG3_F3.tsv",
  F4_Int   = "data/outcomes_harmonized/CDG3_F4.tsv",
  BMI      = "data/outcomes_harmonized/NegCtrl_BMI.tsv",
  Height   = "data/outcomes_harmonized/NegCtrl_height.tsv"
)

# ============================================================================
# IV selection function
# ============================================================================
select_instruments <- function(exposure_path, exposure_name) {
  cat("\n[IV selection] Exposure:", exposure_name, "\n")
  exp <- fread(exposure_path)

  # Genome-wide significant
  exp_gw <- exp %>% filter(P < P_THRESHOLD)
  cat("  Genome-wide sig:", nrow(exp_gw), "SNPs\n")

  # Remove MHC
  exp_gw <- exp_gw %>%
    filter(!(CHR == MHC_REGION$chr & BP >= MHC_REGION$start & BP <= MHC_REGION$end))
  cat("  After MHC removal:", nrow(exp_gw), "\n")

  # Format for TwoSampleMR + LD clump
  exp_fmt <- format_data(
    exp_gw,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "EAF",
    pval_col = "P",
    samplesize_col = "N",
    chr_col = "CHR",
    pos_col = "BP"
  )
  exp_fmt$exposure <- exposure_name

  # Clump (uses MRBase 1KG EUR by default)
  exp_clumped <- clump_data(
    exp_fmt,
    clump_kb = KB_CLUMP,
    clump_r2 = R2_CLUMP,
    pop = "EUR"
  )
  cat("  After clumping (r²<", R2_CLUMP, "):", nrow(exp_clumped), "IVs\n")

  # F-stat
  exp_clumped$f_stat <- (exp_clumped$beta.exposure / exp_clumped$se.exposure)^2
  exp_clumped <- exp_clumped %>% filter(f_stat > F_THRESHOLD)
  cat("  After F>", F_THRESHOLD, ":", nrow(exp_clumped), "IVs (mean F =", round(mean(exp_clumped$f_stat), 1), ")\n")

  invisible(exp_clumped)
}

# ============================================================================
# Outcome extraction
# ============================================================================
extract_outcome <- function(outcome_path, outcome_name, exposure_snps) {
  cat("\n[Outcome extract]", outcome_name, "from", outcome_path, "\n")
  out <- fread(outcome_path)
  out_sub <- out %>% filter(SNP %in% exposure_snps)
  cat("  Found", nrow(out_sub), "/", length(exposure_snps), "matching SNPs\n")

  out_fmt <- format_data(
    out_sub,
    type = "outcome",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "EAF",
    pval_col = "P",
    samplesize_col = "N",
    chr_col = "CHR",
    pos_col = "BP"
  )
  out_fmt$outcome <- outcome_name
  invisible(out_fmt)
}

# ============================================================================
# Core MR per (exposure, outcome) pair
# ============================================================================
run_mr_pair <- function(exp_dat, out_dat, exposure_name, outcome_name) {
  cat("\n========================================\n")
  cat("MR PAIR:", exposure_name, "→", outcome_name, "\n")
  cat("========================================\n")

  # Harmonize
  dat <- harmonise_data(exp_dat, out_dat, action = 2)  # action=2 drops palindromic
  cat("  After harmonization:", nrow(dat), "SNPs\n")

  # Steiger filter
  if (STEIGER_FILTER && nrow(dat) > 1) {
    dat <- steiger_filtering(dat) %>% filter(steiger_dir == TRUE)
    cat("  After Steiger filter:", nrow(dat), "SNPs\n")
  }

  if (nrow(dat) < 3) {
    cat("  WARNING: <3 SNPs, skipping\n")
    return(NULL)
  }

  # Radial-MR outlier removal
  rad <- tryCatch({
    ivw_radial(dat, alpha = 0.05/nrow(dat), weights = 3, tol = 1e-4)
  }, error = function(e) { cat("  Radial failed:", conditionMessage(e), "\n"); NULL })

  if (!is.null(rad) && !is.null(rad$outliers)) {
    outlier_snps <- rad$outliers$SNP
    dat <- dat %>% filter(!(SNP %in% outlier_snps))
    cat("  Radial-MR removed", length(outlier_snps), "outliers; remaining:", nrow(dat), "\n")
  }

  # Primary MR estimators
  res_main <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                                       "mr_weighted_median", "mr_two_sample_ml"))

  # Heterogeneity / pleiotropy diagnostics
  het <- mr_heterogeneity(dat)
  pleio <- mr_pleiotropy_test(dat)

  # MR-PRESSO
  presso <- tryCatch({
    mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
              SdOutcome = "se.outcome", SdExposure = "se.exposure",
              OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
              data = dat, NbDistribution = 5000, SignifThreshold = 0.05, seed = SET_SEED)
  }, error = function(e) { cat("  PRESSO failed:", conditionMessage(e), "\n"); NULL })

  # Leave-one-out
  loo <- mr_leaveoneout(dat)

  list(
    exposure   = exposure_name,
    outcome    = outcome_name,
    n_ivs      = nrow(dat),
    mean_F     = round(mean(dat$f_stat, na.rm = TRUE), 1),
    main       = res_main,
    heterogen  = het,
    pleiotropy = pleio,
    presso     = presso,
    loo        = loo,
    dat        = dat
  )
}

# ============================================================================
# Run grid
# ============================================================================
log_path <- file.path(PROJECT, "logs", sprintf("03_forward_mr_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(dirname(log_path), showWarnings = FALSE, recursive = TRUE)
sink(log_path, split = TRUE)
cat("iter_003 — Forward MR (mtDNA-CN → psychiatric outcomes)\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Seed:", SET_SEED, "\n\n")

all_results <- list()
for (exp_name in names(EXPOSURES)) {
  exp_dat <- select_instruments(EXPOSURES[[exp_name]], exp_name)
  saveRDS(exp_dat, sprintf("results/forward_mr/IVs_%s.rds", exp_name))

  for (out_name in names(OUTCOMES)) {
    out_dat <- extract_outcome(OUTCOMES[[out_name]], out_name, exp_dat$SNP)
    pair_res <- run_mr_pair(exp_dat, out_dat, exp_name, out_name)
    all_results[[sprintf("%s__%s", exp_name, out_name)]] <- pair_res
  }
}

saveRDS(all_results, "results/forward_mr/all_forward_results.rds")

# Build summary table
summary_df <- bind_rows(lapply(all_results, function(r) {
  if (is.null(r)) return(NULL)
  ivw <- r$main %>% filter(method == "Inverse variance weighted")
  data.frame(
    exposure = r$exposure,
    outcome  = r$outcome,
    n_ivs    = r$n_ivs,
    mean_F   = r$mean_F,
    beta     = ivw$b,
    se       = ivw$se,
    p        = ivw$pval,
    Q_p      = r$heterogen$Q_pval[r$heterogen$method == "Inverse variance weighted"],
    egger_intercept_p = r$pleiotropy$pval,
    presso_global_p   = if (!is.null(r$presso)) r$presso$`MR-PRESSO results`$`Global Test`$Pvalue else NA
  )
}))
summary_df <- summary_df %>% mutate(q_BH = p.adjust(p, method = "BH"))
fwrite(summary_df, "results/forward_mr/SUMMARY_forward_mr.tsv", sep = "\t")

cat("\n[Done] Summary saved to results/forward_mr/SUMMARY_forward_mr.tsv\n")
sink()
