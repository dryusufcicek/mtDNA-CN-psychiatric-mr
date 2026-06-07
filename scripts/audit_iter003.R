#!/usr/bin/env Rscript
# DEEP METHOD AUDIT computations for iter_003 mtDNA-CN x psy MR
suppressMessages({ library(data.table); library(dplyr); library(metafor) })
setwd("/arf/scratch/ycicek/iter_003_mtdna_psy_mr")

cat("################ PART A: CORRELATION-AWARE POOLED SE ################\n\n")

fwd <- fread("results/forward_mr/SUMMARY_forward_mr.tsv")
# PAIRWISE_CONCORDANCE.tsv values (from local results/method_comparison, authoritative)
conc <- data.frame(pearson_r = c(0.903, 0.707, 0.553, 0.733, 0.619, 0.687))

# Mean pairwise Pearson r across the 6 panel pairs (correlation of pooled IVW beta across 25 outcomes)
rbar <- mean(conc$pearson_r)
cat(sprintf("Mean pairwise inter-panel Pearson r (across 25 outcomes) = %.4f\n", rbar))
cat(sprintf("Range r = %.3f to %.3f\n", min(conc$pearson_r), max(conc$pearson_r)))

# Effective number of independent panels under exchangeable correlation rho.
# For k correlated estimates with equicorrelation rho, the variance of the simple
# (or IV-weighted, approx equal weights here) mean is inflated by design effect:
#   DEFF = 1 + (k-1)*rho_mean    (for an unweighted mean of k unit-variance items)
# Effective k: k_eff = k / DEFF
k <- 4
for (rho in c(rbar, 0.5, 0.7, 0.9)) {
  deff <- 1 + (k-1)*rho
  keff <- k/deff
  cat(sprintf("  rho=%.3f -> DEFF=%.3f, k_eff=%.2f\n", rho, deff, keff))
}

cat("\n--- Re-derive ASD/ADHD/F3 pooled estimate under correlation-aware variance ---\n")

# Correct fixed-effect IVW pooling under correlated estimates.
# Naive: Var(pooled) = 1 / sum(w_i),  w_i = 1/se_i^2,  pooled = sum(w_i b_i)/sum(w_i)
# Correlated (equicorrelation rho across the k panels, applied to the *estimate* covariance):
#   Cov(b_i,b_j) = rho * se_i * se_j
#   Let w = vector of IVW weights (w_i = 1/se_i^2) normalized to sum 1: a_i = w_i/sum(w)
#   Var(pooled) = a' Sigma a, where Sigma_ij = rho*se_i*se_j (i!=j), se_i^2 (i==j)
corr_pool <- function(b, se, rho) {
  k <- length(b)
  w <- 1/se^2
  a <- w/sum(w)                       # IVW weights, normalized
  Sigma <- outer(se, se) * rho        # off-diag = rho*se_i*se_j
  diag(Sigma) <- se^2                 # diag = se_i^2
  pooled  <- sum(a*b)
  var_pl  <- as.numeric(t(a) %*% Sigma %*% a)
  se_pl   <- sqrt(var_pl)
  z       <- pooled/se_pl
  p       <- 2*pnorm(-abs(z))
  list(beta=pooled, se=se_pl, z=z, p=p, lo=pooled-1.96*se_pl, hi=pooled+1.96*se_pl)
}

for (oc in c("ASD","ADHD","F3")) {
  sub <- fwd[outcome==oc]
  b  <- sub$beta; se <- sub$se
  cat(sprintf("\n=== %s : per-panel beta(se) ===\n", oc))
  for (i in seq_len(nrow(sub))) cat(sprintf("   %-16s b=%+.4f se=%.4f p=%.4g\n",
        sub$exposure[i], b[i], se[i], sub$p[i]))
  # Naive independent (rho=0) — should reproduce metafor FE/REML since tau2=0
  naive <- corr_pool(b, se, 0)
  cat(sprintf(" rho=0   (naive indep) : beta=%+.4f se=%.4f z=%.3f p=%.4g  CI[%.4f,%.4f]\n",
      naive$beta, naive$se, naive$z, naive$p, naive$lo, naive$hi))
  for (rho in c(rbar, 0.55, 0.70, 0.90, 0.99)) {
    cp <- corr_pool(b, se, rho)
    cat(sprintf(" rho=%.3f             : beta=%+.4f se=%.4f z=%.3f p=%.4g  (SE inflation x%.2f)\n",
        rho, cp$beta, cp$se, cp$z, cp$p, cp$se/naive$se))
  }
  # Bound: treat as a single effective study = the single most significant panel
  best <- sub[which.min(sub$p)]
  cat(sprintf(" UPPER BOUND (single best panel = %s): beta=%+.4f se=%.4f p=%.4g\n",
      best$exposure, best$beta, best$se, best$p))
}

cat("\n\n################ PART B: ASD STEIGER RE-CHECK ################\n\n")

r <- readRDS("results/forward_mr/all_forward_results.rds")
asd_keys <- grep("__ASD$", names(r), value=TRUE)

for (k in asd_keys) {
  d <- r[[k]]$dat
  cat(sprintf("\n=== %s (n=%d retained) ===\n", k, nrow(d)))
  cat(sprintf("  samplesize.outcome: unique = %s\n",
              paste(unique(d$samplesize.outcome), collapse=",")))
  cat(sprintf("  samplesize.exposure range: %s - %s\n",
              format(min(d$samplesize.exposure)), format(max(d$samplesize.exposure))))
  cat(sprintf("  rsq.exposure range: %.4g - %.4g (sum=%.4g)\n",
              min(d$rsq.exposure), max(d$rsq.exposure), sum(d$rsq.exposure)))
  cat(sprintf("  rsq.outcome  range: %.4g - %.4g (sum=%.4g)\n",
              min(d$rsq.outcome), max(d$rsq.outcome), sum(d$rsq.outcome)))
  cat(sprintf("  steiger_dir TRUE: %d / %d ; min steiger_pval=%.3g\n",
              sum(d$steiger_dir), nrow(d), min(d$steiger_pval)))
  # How many instruments have rsq.outcome > rsq.exposure (i.e. WRONG direction)?
  wrong <- sum(d$rsq.outcome >= d$rsq.exposure)
  cat(sprintf("  instruments with rsq.outcome >= rsq.exposure (would-be-removed if unfiltered): %d\n", wrong))
  # Steiger is sensitive to the assumed outcome N + prevalence. Re-test with realistic N
  # and with case/control adjustment. ASD Grove constant N used = 46350 (18381+27969).
}

cat("\n--- Sensitivity: does ASD Steiger direction flip if outcome N is MIS-specified? ---\n")
cat("Steiger compares rsq.exposure vs rsq.outcome. rsq scales with N for a fixed Z.\n")
cat("If outcome N were set too LOW, rsq.outcome shrinks -> fewer flips -> MORE instruments kept (anti-conservative).\n")
cat("If outcome N were set too HIGH, rsq.outcome grows -> potential to drop true instruments.\n\n")

# Recompute rsq.outcome under alternative N to test robustness of the kept set.
# get_r_from_lor-style: TwoSampleMR uses get_r_from_pn for continuous; for binary it can
# use get_r_from_lor if prevalence given. Here outcome was treated as continuous (samplesize only).
library(TwoSampleMR)
for (k in asd_keys) {
  d <- r[[k]]$dat
  # Reproduce rsq.outcome via get_r_from_pn (TwoSampleMR default when no ncase/prevalence)
  rsq_out_orig <- d$rsq.outcome
  # Try realistic effective N for a 18381/27969 case-control = 4*Ncase*Ncont/(Nca+Nco)
  Neff <- 4*18381*27969/(18381+27969)
  rsq_out_neff <- get_r_from_pn(d$pval.outcome, rep(Neff, nrow(d)))^2 * sign(1)  # magnitude
  # direction call under each
  rsq_exp <- d$rsq.exposure
  keep_orig <- rsq_exp > rsq_out_orig
  keep_neff <- rsq_exp > rsq_out_neff
  cat(sprintf("%s: outcome-N used=%s | Neff(case-control)=%.0f\n",
      k, unique(d$samplesize.outcome)[1], Neff))
  cat(sprintf("   kept under N-used: %d/%d ; kept under Neff: %d/%d ; DISAGREEMENTS: %d\n",
      sum(keep_orig), nrow(d), sum(keep_neff), nrow(d), sum(keep_orig != keep_neff)))
}

cat("\n\n################ PART C: REPRODUCE METAFOR vs MANUAL (sanity) ################\n\n")
# Confirm rho=0 manual matches metafor REML output in META_summary.tsv
meta <- fread("results/meta_mr/META_summary.tsv")
for (oc in c("ASD","ADHD","F3")) {
  m <- meta[outcome==oc]
  cat(sprintf("%s metafor: pooled_beta=%+.4f pooled_se=%.4f p=%.4g (tau2=%.3g, I2=%.1f)\n",
      oc, m$pooled_beta, m$pooled_se, m$pooled_p, m$tau2, m$I2))
}
cat("\nDONE.\n")
