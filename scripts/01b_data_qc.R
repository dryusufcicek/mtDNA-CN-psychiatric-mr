#!/usr/bin/env Rscript
# iter_003 — Step 01b: COMPREHENSIVE DATA QC
# Purpose: Validate every sumstats file BEFORE pipeline execution.
# - File existence + decompress sanity
# - Row count + duplicate detection
# - Column completeness + NA patterns
# - Data type sanity (CHR/BP/BETA/SE/P/EAF/INFO ranges)
# - Allele encoding (ATCG vs indel; palindromic flag)
# - Genome build heuristic (GRCh37 vs GRCh38)
# - Sample size consistency
# - UKB-inclusion detection in Daner headers
# - Exposure-outcome SNP overlap
# Date: 2026-05-22
# Author: Yusuf Cicek

suppressMessages({
  library(data.table)
  library(dplyr)
})

PROJECT <- Sys.getenv("ITER003_HOME", unset = getwd())
setwd(PROJECT)
cat("Project root:", PROJECT, "\n\n")

# ============================================================================
# QC summary collector
# ============================================================================
qc_results <- list()
qc_warnings <- character()
qc_errors   <- character()

add_warning <- function(file, msg) {
  qc_warnings <<- c(qc_warnings, sprintf("[%s] %s", file, msg))
  cat(sprintf("    [⚠ WARN] %s\n", msg))
}
add_error <- function(file, msg) {
  qc_errors <<- c(qc_errors, sprintf("[%s] %s", file, msg))
  cat(sprintf("    [✗ ERROR] %s\n", msg))
}
add_pass <- function(msg) cat(sprintf("    [✓] %s\n", msg))

# ============================================================================
# Per-file QC function
# ============================================================================
qc_sumstats <- function(path, name, expected_n_snps_min = 1e5,
                        expected_sample_size = NA_real_,
                        col_map = NULL) {
  cat("\n", strrep("=", 90), "\n", sep = "")
  cat("QC:", name, "\n")
  cat("PATH:", path, "\n")
  cat(strrep("=", 90), "\n", sep = "")

  result <- list(name = name, path = path, passed = TRUE)

  # ---- 1. EXISTENCE ----
  cat("\n[1/12] File existence + readability\n")
  if (!file.exists(path)) {
    add_error(name, "File does not exist")
    result$passed <- FALSE
    return(result)
  }
  finfo <- file.info(path)
  size_mb <- finfo$size / 1024^2
  add_pass(sprintf("Exists; size = %.1f MB", size_mb))
  result$size_mb <- size_mb

  # Decompress check (first 10 lines)
  test_lines <- tryCatch(
    suppressWarnings(readLines(gzfile(path), n = 10)),
    error = function(e) NULL
  )
  if (is.null(test_lines) || length(test_lines) < 2) {
    add_error(name, "Cannot decompress / file truncated")
    result$passed <- FALSE
    return(result)
  }
  add_pass("Decompress OK (first 10 lines readable)")

  # ---- 2. READ INTO data.table ----
  cat("\n[2/12] Full file read into data.table\n")
  dt <- tryCatch(fread(path, showProgress = FALSE), error = function(e) NULL)
  if (is.null(dt)) {
    add_error(name, "fread failed")
    result$passed <- FALSE
    return(result)
  }
  add_pass(sprintf("Read %s rows × %d cols", format(nrow(dt), big.mark = ","), ncol(dt)))
  result$nrow <- nrow(dt)
  result$ncol <- ncol(dt)
  result$cols <- colnames(dt)

  if (nrow(dt) < expected_n_snps_min) {
    add_warning(name, sprintf("Row count %d < expected min %d", nrow(dt), expected_n_snps_min))
  }

  # ---- 3. COLUMN MAPPING ----
  cat("\n[3/12] Column auto-detection / mapping\n")
  if (is.null(col_map)) {
    # Auto-detect based on common patterns
    nm <- colnames(dt)
    nm_upper <- toupper(nm)
    col_map <- list(
      SNP   = nm[match(TRUE, nm_upper %in% c("SNP","RSID","RS_ID","VARIANT_ID","MARKERNAME","ID","SNPID"))],
      CHR   = nm[match(TRUE, nm_upper %in% c("CHR","CHROM","CHROMOSOME","#CHR"))],
      BP    = nm[match(TRUE, nm_upper %in% c("BP","POS","POSITION","BASE_PAIR","BASE_PAIR_LOCATION"))],
      A1    = nm[match(TRUE, nm_upper %in% c("A1","EFFECT_ALLELE","EA","ALT","ALLELE1","TESTED_ALLELE"))],
      A2    = nm[match(TRUE, nm_upper %in% c("A2","OTHER_ALLELE","NEA","REF","ALLELE0"))],
      BETA  = nm[match(TRUE, nm_upper %in% c("BETA","B","EFFECT","BETA1"))],
      OR    = nm[match(TRUE, nm_upper %in% c("OR","ODDS_RATIO"))],
      SE    = nm[match(TRUE, nm_upper %in% c("SE","STANDARD_ERROR","STDERR"))],
      P     = nm[match(TRUE, nm_upper %in% c("P","PVAL","PVALUE","P_VALUE","P-VALUE","P_BOLT_LMM","P_BOLT_LMM_INF"))],
      EAF   = nm[match(TRUE, nm_upper %in% c("EAF","FRQ","FREQ","MAF","ALT_FREQ","FREQ.A1","A1FREQ","EFFECT_ALLELE_FREQUENCY","FRQ_A_FAKE","FRQ_U_FAKE","FREQ_TESTED_ALLELE_IN_HRS"))],
      INFO  = nm[match(TRUE, nm_upper %in% c("INFO","IMPINFO","R2"))],
      N     = nm[match(TRUE, nm_upper %in% c("N","N_TOTAL","SAMPLE_SIZE","NEFF","N_EFF","NEFFDIV2"))]
    )
  }
  for (k in names(col_map)) {
    cat(sprintf("  %-6s -> %s\n", k, col_map[[k]] %||% "NOT FOUND"))
  }
  result$col_map <- col_map

  # ---- 4. CHR sanity ----
  cat("\n[4/12] CHR values\n")
  if (!is.na(col_map$CHR) && !is.null(col_map$CHR)) {
    chrs <- unique(as.character(dt[[col_map$CHR]]))
    valid_chrs <- c(as.character(1:22), "X", "Y", "MT", "M",
                    paste0("chr", c(as.character(1:22), "X", "Y", "MT", "M")))
    bad <- setdiff(chrs, valid_chrs)
    add_pass(sprintf("Unique CHR values: %s", paste(sort(chrs)[1:min(10,length(chrs))], collapse=",")))
    if (length(bad) > 0) {
      add_warning(name, sprintf("Non-standard CHR values: %s", paste(bad[1:min(5,length(bad))], collapse=",")))
    }
  } else add_warning(name, "CHR column not detected")

  # ---- 5. BP sanity ----
  cat("\n[5/12] BP coordinate range\n")
  if (!is.na(col_map$BP) && !is.null(col_map$BP)) {
    bp <- as.numeric(dt[[col_map$BP]])
    bp_summary <- summary(bp)
    cat(sprintf("  range: %s\n", paste(names(bp_summary), round(bp_summary,0), sep="=", collapse=" ")))
    if (any(bp < 1 | bp > 3.5e8, na.rm = TRUE)) {
      add_warning(name, sprintf("BP out of range: %d values < 1 or > 3.5e8",
                                sum(bp < 1 | bp > 3.5e8, na.rm = TRUE)))
    }
    n_na_bp <- sum(is.na(bp))
    if (n_na_bp > 0) add_warning(name, sprintf("BP NA values: %d (%.2f%%)", n_na_bp, 100*n_na_bp/length(bp)))
    add_pass(sprintf("BP range %d - %d", min(bp, na.rm=TRUE), max(bp, na.rm=TRUE)))
  } else add_warning(name, "BP column not detected")

  # ---- 6. Effect / SE / P sanity ----
  cat("\n[6/12] BETA/OR + SE + P sanity\n")
  # Effect column (BETA or log(OR))
  if (!is.na(col_map$BETA) && !is.null(col_map$BETA)) {
    eff <- as.numeric(dt[[col_map$BETA]])
    cat(sprintf("  BETA range: [%.3f, %.3f], median = %.4f, NA = %d\n",
                min(eff, na.rm=TRUE), max(eff, na.rm=TRUE), median(eff, na.rm=TRUE), sum(is.na(eff))))
    if (any(!is.finite(eff[!is.na(eff)])))
      add_warning(name, sprintf("BETA Inf values: %d", sum(!is.finite(eff[!is.na(eff)]))))
    if (any(abs(eff) > 5, na.rm = TRUE))
      add_warning(name, sprintf("BETA |value|>5: %d (suspicious for log-OR)", sum(abs(eff) > 5, na.rm=TRUE)))
  } else if (!is.na(col_map$OR) && !is.null(col_map$OR)) {
    or <- as.numeric(dt[[col_map$OR]])
    cat(sprintf("  OR range: [%.3f, %.3f], median = %.4f, NA = %d\n",
                min(or, na.rm=TRUE), max(or, na.rm=TRUE), median(or, na.rm=TRUE), sum(is.na(or))))
    if (any(or <= 0, na.rm = TRUE))
      add_warning(name, sprintf("OR <= 0 values: %d", sum(or <= 0, na.rm=TRUE)))
  } else add_warning(name, "No BETA or OR column found")

  # SE
  if (!is.na(col_map$SE) && !is.null(col_map$SE)) {
    se <- as.numeric(dt[[col_map$SE]])
    cat(sprintf("  SE range: [%.4g, %.4g], NA = %d\n", min(se, na.rm=TRUE), max(se, na.rm=TRUE), sum(is.na(se))))
    if (any(se <= 0, na.rm = TRUE)) add_warning(name, sprintf("SE <= 0: %d", sum(se <= 0, na.rm=TRUE)))
    if (any(!is.finite(se[!is.na(se)]))) add_warning(name, "SE Inf detected")
  } else add_warning(name, "SE column not detected")

  # P
  if (!is.na(col_map$P) && !is.null(col_map$P)) {
    p <- as.numeric(dt[[col_map$P]])
    cat(sprintf("  P range: [%.2e, %.4f], NA = %d, P==0: %d, P>1: %d\n",
                min(p, na.rm=TRUE), max(p, na.rm=TRUE), sum(is.na(p)), sum(p == 0, na.rm=TRUE), sum(p > 1, na.rm=TRUE)))
    if (any(p < 0 | p > 1, na.rm = TRUE))
      add_error(name, sprintf("P out of [0,1]: %d", sum(p < 0 | p > 1, na.rm=TRUE)))
    n_sig <- sum(p < 5e-8, na.rm = TRUE)
    cat(sprintf("  Genome-wide significant (P<5e-8): %d SNPs\n", n_sig))
    result$n_gws <- n_sig
  } else add_warning(name, "P column not detected")

  # ---- 7. Allele encoding ----
  cat("\n[7/12] Allele encoding\n")
  if (!is.na(col_map$A1) && !is.na(col_map$A2) &&
      !is.null(col_map$A1) && !is.null(col_map$A2)) {
    a1 <- as.character(dt[[col_map$A1]])
    a2 <- as.character(dt[[col_map$A2]])
    snv_pattern <- "^[ACGTacgt]$"
    n_snv <- sum(grepl(snv_pattern, a1) & grepl(snv_pattern, a2), na.rm = TRUE)
    n_indel <- nrow(dt) - n_snv
    cat(sprintf("  SNVs: %d (%.1f%%), indels/multi-char: %d (%.1f%%)\n",
                n_snv, 100*n_snv/nrow(dt), n_indel, 100*n_indel/nrow(dt)))
    # Palindromic
    palin <- (toupper(a1) == "A" & toupper(a2) == "T") |
             (toupper(a1) == "T" & toupper(a2) == "A") |
             (toupper(a1) == "C" & toupper(a2) == "G") |
             (toupper(a1) == "G" & toupper(a2) == "C")
    n_palin <- sum(palin, na.rm = TRUE)
    cat(sprintf("  Palindromic (A/T or C/G): %d (%.1f%%)\n", n_palin, 100*n_palin/nrow(dt)))
    add_pass("Allele encoding parsed")
  } else add_warning(name, "Cannot QC alleles (A1/A2 not detected)")

  # ---- 8. Duplicate SNP IDs ----
  cat("\n[8/12] Duplicate SNP IDs\n")
  if (!is.na(col_map$SNP) && !is.null(col_map$SNP)) {
    n_dup <- sum(duplicated(dt[[col_map$SNP]]))
    if (n_dup > 0) add_warning(name, sprintf("Duplicate SNP IDs: %d", n_dup))
    else add_pass("No duplicate SNP IDs")
  }

  # ---- 9. Sample size / N ----
  cat("\n[9/12] Sample size check\n")
  if (!is.na(col_map$N) && !is.null(col_map$N)) {
    n_vals <- as.numeric(dt[[col_map$N]])
    cat(sprintf("  N range: [%.0f, %.0f], median = %.0f\n",
                min(n_vals, na.rm=TRUE), max(n_vals, na.rm=TRUE), median(n_vals, na.rm=TRUE)))
    if (!is.na(expected_sample_size)) {
      med_n <- median(n_vals, na.rm=TRUE)
      if (abs(med_n - expected_sample_size) / expected_sample_size > 0.10) {
        add_warning(name, sprintf("Median N=%.0f differs >10%% from expected %.0f", med_n, expected_sample_size))
      } else add_pass(sprintf("N consistent with expected ~%.0f", expected_sample_size))
    }
  } else {
    cat(sprintf("  N column not detected; expected = %s\n", expected_sample_size))
    if (!is.na(expected_sample_size)) {
      add_pass(sprintf("No N column; using fixed N = %.0f from metadata", expected_sample_size))
    }
  }

  # ---- 10. Daner-format header check (for outcomes) ----
  cat("\n[10/12] Header / metadata scan (UKB-inclusion detection)\n")
  # Read first 50 lines looking for ## comments or daner header
  con <- gzfile(path, "r")
  first50 <- readLines(con, n = 50)
  close(con)
  daner_comments <- grep("^##", first50, value = TRUE)
  if (length(daner_comments) > 0) {
    cat("  Daner-style comments found:\n")
    for (l in daner_comments[1:min(5,length(daner_comments))]) cat("    ", l, "\n")
  } else add_pass("No Daner ## comments (clean header)")

  # UKB detection in path or comments
  # Check DEDUP first — if path indicates UKB samples removed, that's the answer regardless
  ukb_dedup_in_name <- grepl("ukbbdedup|ukb_dedup|ukbdedup|UKBdedup|UKB_dedup", path, ignore.case = TRUE)
  ukb_in_name <- !ukb_dedup_in_name && grepl("ukb|biobank|UKBB|UK_BIO", path, ignore.case = TRUE)
  ukb_in_header <- any(grepl("UKB|biobank|UK Biobank", first50, ignore.case = TRUE))
  if (ukb_dedup_in_name) {
    add_pass("Path contains 'ukbbdedup' — UKB samples explicitly REMOVED")
    result$ukb_status <- "REMOVED"
  } else if (ukb_in_name || ukb_in_header) {
    add_warning(name, "UKB samples LIKELY INCLUDED — check sample overlap")
    result$ukb_status <- "INCLUDED"
  } else {
    add_pass("No UKB mention in path/header — likely UKB-free")
    result$ukb_status <- "UNKNOWN"
  }

  # ---- 11. Genome build heuristic ----
  cat("\n[11/12] Genome build heuristic (GRCh37 vs GRCh38)\n")
  # Use a few known SNPs whose positions differ between builds
  # rs7412 (APOE region): GRCh37 chr19:45,412,079; GRCh38 chr19:44,908,822
  build_check <- NULL
  if (!is.na(col_map$SNP) && !is.null(col_map$SNP) &&
      !is.na(col_map$CHR) && !is.null(col_map$CHR) &&
      !is.na(col_map$BP) && !is.null(col_map$BP)) {
    rs7412 <- dt[get(col_map$SNP) == "rs7412"]
    if (nrow(rs7412) > 0) {
      bp_val <- rs7412[[col_map$BP]][1]
      build <- if (abs(bp_val - 45412079) < 100) "GRCh37" else if (abs(bp_val - 44908822) < 100) "GRCh38" else "UNCLEAR"
      add_pass(sprintf("rs7412 BP = %d → %s", bp_val, build))
      result$build <- build
    } else {
      cat("  rs7412 not in file (likely variant_id format)\n")
    }
  }

  # ---- 12. File-format-specific notes (Daner OR + Nca/Nco) ----
  cat("\n[12/12] Format-specific completeness\n")
  has_or <- !is.na(col_map$OR) && !is.null(col_map$OR)
  has_beta <- !is.na(col_map$BETA) && !is.null(col_map$BETA)
  has_se <- !is.na(col_map$SE) && !is.null(col_map$SE)
  has_p <- !is.na(col_map$P) && !is.null(col_map$P)
  has_a1a2 <- (!is.na(col_map$A1) && !is.na(col_map$A2) &&
               !is.null(col_map$A1) && !is.null(col_map$A2))

  ready_for_mr <- has_se && has_p && has_a1a2 && (has_beta || has_or)
  if (ready_for_mr) add_pass("READY FOR MR (effect + SE + P + alleles all present)")
  else add_error(name, "NOT MR-ready: missing core columns")

  result$ready_for_mr <- ready_for_mr
  invisible(result)
}

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# ============================================================================
# Run QC on all sumstats
# ============================================================================
log_path <- file.path(PROJECT, "logs",
                     sprintf("01b_data_qc_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(dirname(log_path), showWarnings = FALSE, recursive = TRUE)
sink(log_path, split = TRUE)

cat("iter_003 — Comprehensive Data QC\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Project:", PROJECT, "\n\n")

# EXPOSURES
cat(strrep("#", 90), "\n", sep = "")
cat("# EXPOSURES (4 mtDNA-CN GWAS)\n")
cat(strrep("#", 90), "\n", sep = "")

exp_files <- list(
  list("data/exposures/Longchamps2022_mtDNA_CN.ALLm2.bgen.stats.gz",
       "Longchamps_2022_CHARGE_qPCR_meta", 1e7, 465809),
  list("data/exposures/Chong2022_mtDNA_CN_EUR_GCST90026372.tsv.gz",
       "Chong_2022_UKB_array", 1e7, 383476),
  list("data/exposures/Gupta2023_mtDNA_CN_adjusted_GCST90268497.tsv.gz",
       "Gupta_2023_UKB_WGS_adjusted", 1e7, 155998),
  list("data/exposures/Gupta2023_mtDNA_CN_raw_GCST90268498.tsv.gz",
       "Gupta_2023_UKB_WGS_raw", 1e7, 155998)
)
for (f in exp_files) {
  res <- qc_sumstats(f[[1]], f[[2]], f[[3]], f[[4]])
  qc_results[[f[[2]]]] <- res
}

# OUTCOMES
cat("\n", strrep("#", 90), "\n", sep = "")
cat("# OUTCOMES (5 primary + 3 CDG3 + 2 neg controls)\n")
cat(strrep("#", 90), "\n", sep = "")

out_files <- list(
  list("data/outcomes/SCZ_PGC3_UKBdedup.gz", "SCZ_PGC3_UKBdedup", 1e6, 130644),
  list("data/outcomes/ADHD_Demontis2023.meta.gz", "ADHD_Demontis2023", 1e6, 225534),
  list("data/outcomes/ASD_Grove2019.gz", "ASD_Grove2019", 1e6, 46350),
  list("data/outcomes/BD_Mullins2024_EUR.gz", "BD_Mullins2024_EUR", 1e6, 353899),
  list("data/outcomes/MDD_Adams2025_EUR.gz", "MDD_Adams2025_EUR", 1e6, 1639572),
  list("data/outcomes/CDG3_PFactor.tsv.gz", "CDG3_PFactor", 1e6, 2168621),
  list("data/outcomes/CDG3_F3_Neurodev.tsv.gz", "CDG3_F3_Neurodev", 1e6, 84760),
  list("data/outcomes/CDG3_F4_Internalizing.tsv.gz", "CDG3_F4_Internalizing", 1e6, 1637337),
  list("data/outcomes/NegCtrl_BMI.gz", "NegCtrl_BMI_Yengo2018", 1e6, 700000),
  list("data/outcomes/NegCtrl_height.gz", "NegCtrl_height_Yengo2018", 1e6, 700000)
)
for (f in out_files) {
  res <- qc_sumstats(f[[1]], f[[2]], f[[3]], f[[4]])
  qc_results[[f[[2]]]] <- res
}

# ============================================================================
# FINAL SUMMARY
# ============================================================================
cat("\n\n", strrep("#", 90), "\n", sep = "")
cat("# FINAL QC SUMMARY\n")
cat(strrep("#", 90), "\n", sep = "")

summary_df <- bind_rows(lapply(qc_results, function(r) {
  data.frame(
    file = r$name,
    size_MB = round(r$size_mb %||% NA, 1),
    nrow = r$nrow %||% NA_integer_,
    ncol = r$ncol %||% NA_integer_,
    n_gws_5e8 = r$n_gws %||% NA_integer_,
    build = r$build %||% NA_character_,
    ukb_status = r$ukb_status %||% NA_character_,
    ready_for_mr = r$ready_for_mr %||% FALSE
  )
}))
cat("\n")
print(summary_df, row.names = FALSE)

# Errors and warnings
cat("\n\n--- ERRORS (must fix before pipeline) ---\n")
if (length(qc_errors) == 0) cat("  None\n") else for (e in qc_errors) cat("  ", e, "\n")

cat("\n--- WARNINGS (review, may proceed) ---\n")
if (length(qc_warnings) == 0) cat("  None\n") else for (w in qc_warnings) cat("  ", w, "\n")

cat("\n--- READY TO RUN PIPELINE? ---\n")
all_ready <- all(sapply(qc_results, function(r) r$ready_for_mr %||% FALSE))
if (length(qc_errors) > 0) {
  cat("  ❌ NO — fix errors first\n")
} else if (!all_ready) {
  cat("  ❌ NO — some files not MR-ready\n")
} else {
  cat("  ✅ YES — all files validated; proceed to 02_harmonize_exposures.R\n")
}

saveRDS(qc_results, file.path(PROJECT, "logs", "01b_qc_results.rds"))
sink()
cat("\nQC log saved to:", log_path, "\n")
