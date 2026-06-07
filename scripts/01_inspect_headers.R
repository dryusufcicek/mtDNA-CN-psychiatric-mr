#!/usr/bin/env Rscript
# iter_003 — Step 01: Header inspection for all sumstats
# Purpose: Verify column structure, allele coding, effect direction, INFO availability,
#          sample size, genome build for all 4 exposures + 8 outcomes
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
# Helper function
# ============================================================================
inspect_sumstats <- function(path, name) {
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("FILE:", name, "\n")
  cat("PATH:", path, "\n")
  cat(strrep("=", 80), "\n", sep = "")

  # File size
  finfo <- file.info(path)
  cat(sprintf("Size: %.1f MB\n", finfo$size / 1024^2))

  # First line (potentially header or comment)
  con <- gzfile(path, "r")
  first_lines <- readLines(con, n = 5)
  close(con)
  cat("\nFirst 5 lines:\n")
  for (i in seq_along(first_lines)) cat(sprintf("  [%d] %s\n", i, substr(first_lines[i], 1, 200)))

  # Try fread
  tryCatch({
    dt <- fread(path, nrows = 1000)
    cat("\nColumn names (n =", ncol(dt), "):\n")
    cat(paste0("  ", colnames(dt), collapse = "\n"), "\n")
    cat("\nFirst row preview:\n")
    print(head(dt, 2))

    # Count total rows
    cat("\nCounting total rows (may take a moment for large files)...\n")
    total_rows <- length(count.fields(path, sep = "\t"))
    cat(sprintf("Total rows (incl. header): %d\n", total_rows))

    # Detect key columns
    nm <- toupper(colnames(dt))
    has_snp <- any(c("SNP", "RSID", "RS_ID", "VARIANT_ID", "MARKERNAME", "ID") %in% nm)
    has_chr <- any(c("CHR", "CHROM", "CHROMOSOME", "#CHR") %in% nm)
    has_pos <- any(c("BP", "POS", "POSITION", "BASE_PAIR") %in% nm)
    has_a1  <- any(c("A1", "EFFECT_ALLELE", "EA", "ALT") %in% nm)
    has_a2  <- any(c("A2", "OTHER_ALLELE", "NEA", "REF") %in% nm)
    has_beta <- any(c("BETA", "B", "EFFECT") %in% nm)
    has_or <- any(c("OR", "ODDS_RATIO") %in% nm)
    has_se <- any(c("SE", "STANDARD_ERROR", "STDERR") %in% nm)
    has_p <- any(c("P", "PVAL", "PVALUE", "P_VALUE", "P-VALUE") %in% nm)
    has_n <- any(c("N", "N_TOTAL", "SAMPLE_SIZE", "NEFF", "N_EFF") %in% nm)
    has_freq <- any(c("EAF", "FRQ", "FREQ", "MAF", "ALT_FREQ") %in% nm)
    has_info <- any(c("INFO", "IMPINFO", "R2") %in% nm)

    cat("\nKey-column detection:\n")
    cat(sprintf("  SNP/rsID: %s\n  CHR:     %s\n  POS:     %s\n  A1/EA:   %s\n  A2/NEA:  %s\n  BETA:    %s\n  OR:      %s\n  SE:      %s\n  P:       %s\n  N:       %s\n  FREQ:    %s\n  INFO:    %s\n",
                has_snp, has_chr, has_pos, has_a1, has_a2, has_beta, has_or, has_se, has_p, has_n, has_freq, has_info))

    # Check sample-size statements in headers (Daner-style)
    daner_comment <- first_lines[grepl("^#", first_lines)]
    if (length(daner_comment) > 0) {
      cat("\nDaner-style header comments:\n")
      for (l in daner_comment) cat("  ", l, "\n")
    }

    invisible(list(name = name, ncol = ncol(dt), nrow = total_rows,
                   has_snp = has_snp, has_chr = has_chr, has_pos = has_pos,
                   has_a1 = has_a1, has_a2 = has_a2, has_beta = has_beta,
                   has_or = has_or, has_se = has_se, has_p = has_p,
                   has_n = has_n, has_freq = has_freq, has_info = has_info))
  }, error = function(e) {
    cat("\nfread ERROR:", conditionMessage(e), "\n")
    invisible(NULL)
  })
}

# ============================================================================
# Inspection batch
# ============================================================================
log_path <- file.path(PROJECT, "logs", sprintf("format_inspection_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(dirname(log_path), showWarnings = FALSE, recursive = TRUE)

sink(log_path, split = TRUE)
cat("iter_003 — Sumstats header inspection\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# EXPOSURES (4 mtDNA-CN GWAS)
cat("\n", strrep("#", 80), "\n", sep = "")
cat("EXPOSURES (mtDNA-CN GWAS)\n")
cat(strrep("#", 80), "\n", sep = "")

exposures <- list(
  list("data/exposures/Longchamps2022_mtDNA_CN.ALLm2.bgen.stats.gz", "Longchamps_2022_CHARGE_qPCR_meta"),
  list("data/exposures/Chong2022_mtDNA_CN_EUR_GCST90026372.tsv.gz",  "Chong_2022_UKB_array"),
  list("data/exposures/Gupta2023_mtDNA_CN_adjusted_GCST90268497.tsv.gz", "Gupta_2023_UKB_WGS_adjusted"),
  list("data/exposures/Gupta2023_mtDNA_CN_raw_GCST90268498.tsv.gz",  "Gupta_2023_UKB_WGS_raw")
)
for (e in exposures) inspect_sumstats(e[[1]], e[[2]])

# OUTCOMES (5 + 3 + 2)
cat("\n", strrep("#", 80), "\n", sep = "")
cat("OUTCOMES (psychiatric disorders + bonus + neg controls)\n")
cat(strrep("#", 80), "\n", sep = "")

outcomes <- list(
  list("data/outcomes/SCZ_PGC3_UKBdedup.gz", "SCZ_PGC3_UKBdedup"),
  list("data/outcomes/ADHD_Demontis2023.meta.gz", "ADHD_Demontis2023"),
  list("data/outcomes/ASD_Grove2019.gz", "ASD_Grove2019"),
  list("data/outcomes/BD_Mullins2024_EUR.gz", "BD_Mullins2024_EUR"),
  list("data/outcomes/MDD_Adams2025_EUR.gz", "MDD_Adams2025_EUR"),
  list("data/outcomes/CDG3_PFactor.tsv.gz", "CDG3_PFactor"),
  list("data/outcomes/CDG3_F3_Neurodev.tsv.gz", "CDG3_F3_Neurodev"),
  list("data/outcomes/CDG3_F4_Internalizing.tsv.gz", "CDG3_F4_Internalizing"),
  list("data/outcomes/NegCtrl_BMI.gz", "NegCtrl_BMI_Yengo2018"),
  list("data/outcomes/NegCtrl_height.gz", "NegCtrl_height_Yengo2018")
)
for (o in outcomes) inspect_sumstats(o[[1]], o[[2]])

sink()

cat("\n\n[ Inspection log saved to:", log_path, "]\n")
