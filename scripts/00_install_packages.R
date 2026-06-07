#!/usr/bin/env Rscript
# iter_003 — Step 00: Install all R packages required for pipeline
# Run once before any other script
# Date: 2026-05-22

# ============================================================================
# CRAN packages
# ============================================================================
cran_pkgs <- c(
  "data.table",     # FAST I/O ✓
  "dplyr",          # tidy ops ✓
  "ggplot2",        # plotting ✓
  "MendelianRandomization",  # alternative MR estimators
  "metafor",        # random-effects meta-analysis
  "remotes",        # for github installs
  "devtools",       # backup install method
  "BiocManager"     # for VariantAnnotation if needed
)
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("Installing", p, "...\n")
    install.packages(p, repos = "https://cloud.r-project.org/")
  } else {
    cat(p, "already installed.\n")
  }
}

# ============================================================================
# GitHub packages (MR Base ecosystem)
# ============================================================================
gh_pkgs <- list(
  TwoSampleMR = "MRCIEU/TwoSampleMR",
  MRPRESSO    = "rondolab/MR-PRESSO",
  RadialMR    = "WSpiller/RadialMR",
  ieugwasr    = "MRCIEU/ieugwasr"   # for clumping API
)
for (p in names(gh_pkgs)) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("Installing", p, "from GitHub:", gh_pkgs[[p]], "...\n")
    remotes::install_github(gh_pkgs[[p]], upgrade = "never")
  } else {
    cat(p, "already installed.\n")
  }
}

# ============================================================================
# Final verification
# ============================================================================
cat("\n========== FINAL PACKAGE STATUS ==========\n")
all_pkgs <- c(cran_pkgs, names(gh_pkgs))
for (p in all_pkgs) {
  status <- if (requireNamespace(p, quietly = TRUE)) "✓ INSTALLED" else "✗ MISSING"
  cat(sprintf("  %-25s : %s\n", p, status))
}

# Print sessionInfo
cat("\n\n========== SESSION INFO ==========\n")
sessionInfo()
