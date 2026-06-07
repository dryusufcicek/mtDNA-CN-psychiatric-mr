#!/usr/bin/env Rscript
# iter_003 — Step 04: Multi-method meta-MR
# Per outcome: pool 4 IVW estimates from 4 mtDNA-CN exposures using random-effects meta
# Compute I² between datasets — H1 primary endpoint
# Date: 2026-05-22

suppressMessages({
  library(data.table)
  library(dplyr)
  library(metafor)
  library(ggplot2)
})

PROJECT <- Sys.getenv("ITER003_HOME", unset = getwd())
setwd(PROJECT)
dir.create("results/meta_mr", showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Load forward MR summary
# ============================================================================
fwd <- fread("results/forward_mr/SUMMARY_forward_mr.tsv")

# ============================================================================
# Meta-analyze per outcome
# ============================================================================
outcomes <- unique(fwd$outcome)
meta_results <- list()

for (out in outcomes) {
  sub <- fwd %>% filter(outcome == out)
  if (nrow(sub) < 2) {
    cat("Skipping", out, "— only", nrow(sub), "estimates\n")
    next
  }

  # Random-effects meta
  res <- rma(yi = beta, sei = se, slab = exposure, data = sub, method = "REML")

  meta_results[[out]] <- list(
    outcome       = out,
    n_studies     = nrow(sub),
    pooled_beta   = res$beta[1, 1],
    pooled_se     = res$se,
    pooled_p      = res$pval,
    pooled_lo     = res$ci.lb,
    pooled_hi     = res$ci.ub,
    Q             = res$QE,
    Q_p           = res$QEp,
    I2            = res$I2,
    tau2          = res$tau2,
    rma_obj       = res
  )

  # Forest plot
  pdf(sprintf("results/meta_mr/forest_%s.pdf", out), width = 8, height = 5)
  forest(res, header = c("mtDNA-CN GWAS", "OR [95% CI]"), atransf = exp,
         main = sprintf("mtDNA-CN → %s (multi-method meta-MR)", out))
  dev.off()
}

# ============================================================================
# Compile master meta table
# ============================================================================
meta_df <- bind_rows(lapply(meta_results, function(m) {
  data.frame(
    outcome = m$outcome,
    n_studies = m$n_studies,
    pooled_OR = exp(m$pooled_beta),
    pooled_OR_lo = exp(m$pooled_lo),
    pooled_OR_hi = exp(m$pooled_hi),
    pooled_p = m$pooled_p,
    Q = m$Q,
    Q_p = m$Q_p,
    I2 = m$I2,
    tau2 = m$tau2
  )
}))
meta_df <- meta_df %>% mutate(q_BH = p.adjust(pooled_p, method = "BH"))
fwrite(meta_df, "results/meta_mr/META_summary.tsv", sep = "\t")

# ============================================================================
# H1 endpoint: I² > 50% in any disorder?
# ============================================================================
cat("\n========================================\n")
cat("H1 ENDPOINT: Between-dataset I² heterogeneity\n")
cat("========================================\n")
cat(sprintf("%-15s %8s %8s %8s\n", "Outcome", "I²(%)", "Q", "Q_p"))
for (i in seq_len(nrow(meta_df))) {
  cat(sprintf("%-15s %8.1f %8.2f %8.3g\n",
              meta_df$outcome[i], meta_df$I2[i], meta_df$Q[i], meta_df$Q_p[i]))
}

h1_disorders <- meta_df %>% filter(I2 > 50, !outcome %in% c("BMI", "Height"))
cat("\nH1 status:",
    ifelse(nrow(h1_disorders) > 0,
           sprintf("✓ SUPPORTED — I²>50%% in %d psychiatric outcome(s): %s",
                   nrow(h1_disorders), paste(h1_disorders$outcome, collapse = ", ")),
           "✗ FALSIFIED — no psychiatric outcome with I²>50%"),
    "\n")

saveRDS(meta_results, "results/meta_mr/META_full_results.rds")
