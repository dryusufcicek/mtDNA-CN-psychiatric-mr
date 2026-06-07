#!/usr/bin/env Rscript
# iter_003 — Step 02: Harmonize 4 mtDNA-CN exposures to common schema
# Output: data/exposures_harmonized/*.tsv (uniform columns: SNP, CHR, BP, A1, A2, BETA, SE, P, EAF, N, INFO)
# Date: 2026-05-22

suppressMessages({
  library(data.table)
  library(dplyr)
})

PROJECT <- Sys.getenv("ITER003_HOME", unset = getwd())
setwd(PROJECT)
cat("Project root:", PROJECT, "\n")
dir.create("data/exposures_harmonized", showWarnings = FALSE)

# ============================================================================
# Define expected schemas per exposure (FILL IN AFTER 01_inspect_headers.R)
# ============================================================================
# TEMPLATE — actual column names to be determined post-inspection

# Longchamps 2022 ALLm2.bgen.stats (PLINK/BOLT-LMM output likely)
# Expected cols: SNP, CHR, BP, ALLELE1, ALLELE0, A1FREQ, INFO, BETA, SE, P_BOLT_LMM
harmonize_longchamps <- function() {
  # VERIFIED COLUMNS (2026-05-22 header inspection):
  # SNP CHR BP GENPOS ALLELE1 ALLELE0 A1FREQ INFO CHISQ_LINREG P_LINREG
  # BETA SE CHISQ_BOLT_LMM_INF P_BOLT_LMM_INF CHISQ_BOLT_LMM P_BOLT_LMM
  # ALLELE1 = effect allele (BOLT-LMM convention)
  dt <- fread("data/exposures/Longchamps2022_mtDNA_CN.ALLm2.bgen.stats.gz")
  out <- dt %>%
    filter(INFO >= 0.8) %>%       # imputation quality filter
    filter(CHR %in% 1:22) %>%     # autosomal only (drop X/Y/MT)
    rename(
      A1 = ALLELE1,                # effect allele
      A2 = ALLELE0,
      P = P_BOLT_LMM,              # use BOLT-LMM mixed model p-value
      EAF = A1FREQ
    ) %>%
    distinct(SNP, .keep_all = TRUE) %>%  # drop duplicate SNP IDs (multi-allelic)
    mutate(N = 465809, exposure = "Longchamps2022_CHARGE_qPCR") %>%
    select(SNP, CHR, BP, A1, A2, BETA, SE, P, EAF, N, INFO, exposure)
  fwrite(out, "data/exposures_harmonized/Longchamps2022.tsv", sep = "\t")
  cat("Longchamps: ", nrow(out), "SNPs harmonized (after INFO>=0.8, autosomal, dedup)\n")
  invisible(out)
}

# Chong 2022 GWAS Catalog format
harmonize_chong <- function() {
  # VERIFIED COLUMNS (2026-05-22): variant_id chromosome base_pair_location
  # effect_allele other_allele beta standard_error p_value
  # NOTE: Chong's GWAS Catalog harmonized format does NOT include EAF or INFO
  dt <- fread("data/exposures/Chong2022_mtDNA_CN_EUR_GCST90026372.tsv.gz")
  out <- dt %>%
    rename(
      SNP = variant_id,
      CHR = chromosome,
      BP = base_pair_location,
      A1 = effect_allele,
      A2 = other_allele,
      BETA = beta,
      SE = standard_error,
      P = p_value
    ) %>%
    filter(CHR %in% 1:22) %>%
    distinct(SNP, .keep_all = TRUE) %>%
    mutate(EAF = NA_real_, N = 383476, INFO = NA_real_, exposure = "Chong2022_UKB_array") %>%
    select(SNP, CHR, BP, A1, A2, BETA, SE, P, EAF, N, INFO, exposure)
  fwrite(out, "data/exposures_harmonized/Chong2022.tsv", sep = "\t")
  cat("Chong: ", nrow(out), "SNPs harmonized (autosomal, dedup)\n")
  invisible(out)
}

# Gupta 2023 GWAS Catalog format (Pan-UKB-derived schema)
# Helper: translate Gupta variant_id (chr_bp_a1_a2) -> rsID using 1KG lookup
.translate_gupta_to_rsid <- function(dt) {
  lookup_file <- "data/rsid_lookup_1KG_EUR.rds"
  if (!file.exists(lookup_file)) {
    stop(sprintf("rsID lookup missing: %s — run 01c_build_rsid_lookup.R first", lookup_file))
  }
  cat("Loading rsID lookup...\n")
  lookup_vec <- readRDS(lookup_file)
  cat(sprintf("  %d lookup keys loaded\n", length(lookup_vec)))

  # Build key from existing CHR/BP/A1/A2 columns
  dt[, key := sprintf("%s_%s_%s_%s", CHR, BP, A1, A2)]
  dt[, rsID := lookup_vec[key]]
  n_mapped <- sum(!is.na(dt$rsID))
  cat(sprintf("  %d / %d (%.1f%%) Gupta variants mapped to rsID\n",
              n_mapped, nrow(dt), 100*n_mapped/nrow(dt)))
  # Keep only mapped variants
  dt <- dt[!is.na(rsID)]
  dt[, SNP := rsID]
  dt[, c("key", "rsID") := NULL]
  dt
}

harmonize_gupta_adjusted <- function() {
  # VERIFIED COLUMNS (2026-05-22):
  # chromosome base_pair_location effect_allele other_allele beta standard_error
  # effect_allele_frequency p_value variant_id pval_heterogeneity
  # NOTE: 03/2024 author correction applied — effect_allele/other_allele are now correct
  # CHR=23 represents X chromosome (numeric encoding) — dropped via autosomal filter
  # variant_id format "chr_bp_ref_alt" — translated to rsID via 1KG lookup
  dt <- fread("data/exposures/Gupta2023_mtDNA_CN_adjusted_GCST90268497.tsv.gz")
  dt <- dt[chromosome %in% 1:22]
  setnames(dt,
           c("variant_id","chromosome","base_pair_location","effect_allele","other_allele",
             "beta","standard_error","p_value","effect_allele_frequency","pval_heterogeneity"),
           c("SNP","CHR","BP","A1","A2","BETA","SE","P","EAF","P_HET"))
  dt <- .translate_gupta_to_rsid(dt)
  dt <- unique(dt, by = "SNP")
  dt[, N := 155998][, INFO := NA_real_][, exposure := "Gupta2023_UKB_WGS_adjusted"]
  out <- dt[, .(SNP, CHR, BP, A1, A2, BETA, SE, P, EAF, N, INFO, exposure)]
  fwrite(out, "data/exposures_harmonized/Gupta2023_adjusted.tsv", sep = "\t")
  cat("Gupta adjusted: ", nrow(out), "SNPs harmonized (autosomal, rsID-mapped)\n")
  invisible(out)
}

harmonize_gupta_raw <- function() {
  dt <- fread("data/exposures/Gupta2023_mtDNA_CN_raw_GCST90268498.tsv.gz")
  dt <- dt[chromosome %in% 1:22]
  setnames(dt,
           c("variant_id","chromosome","base_pair_location","effect_allele","other_allele",
             "beta","standard_error","p_value","effect_allele_frequency"),
           c("SNP","CHR","BP","A1","A2","BETA","SE","P","EAF"))
  dt <- .translate_gupta_to_rsid(dt)
  dt <- unique(dt, by = "SNP")
  dt[, N := 155998][, INFO := NA_real_][, exposure := "Gupta2023_UKB_WGS_raw"]
  out <- dt[, .(SNP, CHR, BP, A1, A2, BETA, SE, P, EAF, N, INFO, exposure)]
  fwrite(out, "data/exposures_harmonized/Gupta2023_raw.tsv", sep = "\t")
  cat("Gupta raw: ", nrow(out), "SNPs harmonized (autosomal, rsID-mapped)\n")
  invisible(out)
}

# ============================================================================
# Run all
# ============================================================================
log_path <- file.path(PROJECT, "logs", sprintf("02_harmonize_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
sink(log_path, split = TRUE)
cat("iter_003 — Exposure harmonization\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Run ONLY AFTER 01_inspect_headers.R confirmed column names
# harmonize_longchamps()
# harmonize_chong()
# harmonize_gupta_adjusted()
# harmonize_gupta_raw()

cat("\n[STUB] Uncomment harmonize_* calls after column-name verification in 01_inspect_headers.R log\n")
sink()
