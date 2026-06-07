#!/usr/bin/env Rscript
suppressMessages({ library(data.table) })
setwd("/arf/scratch/ycicek/iter_003_mtdna_psy_mr")

cat("################ PART D: HIERARCHICAL FDR UNDER CORRELATION CORRECTION ################\n\n")
meta <- fread("results/meta_mr/META_summary.tsv")

# Primary-tier family of 5: ADHD ASD BD MDD SCZ
prim <- c("ADHD","ASD","BD","MDD","SCZ")
# Correlation-aware p (rho=0.70) recomputed for the protective ones; null ones keep naive
# (correlation correction only widens SE -> null stays null)
fwd <- fread("results/forward_mr/SUMMARY_forward_mr.tsv")
corr_p <- function(oc, rho=0.70) {
  sub <- fwd[outcome==oc]; b<-sub$beta; se<-sub$se
  w<-1/se^2; a<-w/sum(w)
  S<-outer(se,se)*rho; diag(S)<-se^2
  pooled<-sum(a*b); v<-as.numeric(t(a)%*%S%*%a); z<-pooled/sqrt(v)
  2*pnorm(-abs(z))
}
cat("--- PRIMARY TIER (k=5) BH-FDR ---\n")
p_naive <- sapply(prim, function(o) meta[outcome==o]$pooled_p)
p_corr  <- sapply(prim, function(o) corr_p(o, 0.70))
q_naive <- p.adjust(p_naive, "BH")
q_corr  <- p.adjust(p_corr, "BH")
for (i in seq_along(prim)) {
  cat(sprintf("  %-5s p_naive=%.4g q_naive=%.4g | p_corr(rho.7)=%.4g q_corr=%.4g  %s->%s\n",
      prim[i], p_naive[i], q_naive[i], p_corr[i], q_corr[i],
      ifelse(q_naive[i]<0.05,"SIG","ns"), ifelse(q_corr[i]<0.05,"SIG","ns")))
}

cat("\n--- FACTOR TIER (k=6) BH-FDR ---\n")
fac <- c("PFactor","F1","F2","F3","F4","F5")
pf_naive <- sapply(fac, function(o) meta[outcome==o]$pooled_p)
qf_naive <- p.adjust(pf_naive, "BH")
# corr only for F3 (others null)
pf_corr <- pf_naive; pf_corr["F3"] <- corr_p("F3",0.70)
qf_corr <- p.adjust(pf_corr, "BH")
for (i in seq_along(fac)) {
  cat(sprintf("  %-8s p_naive=%.4g q_naive=%.4g | p_corr=%.4g q_corr=%.4g  %s->%s\n",
      fac[i], pf_naive[i], qf_naive[i], pf_corr[i], qf_corr[i],
      ifelse(qf_naive[i]<0.05,"SIG","ns"), ifelse(qf_corr[i]<0.05,"SIG","ns")))
}

cat("\n--- POOLED-ACROSS-ALL-FAMILIES sanity (if NOT tiered: single BH over all 19 psychiatric pooled tests) ---\n")
psy <- meta[!outcome %in% c("BMI","Height")]
# Exclude noUKBB duplicates to avoid double counting? Report both ways.
allp <- psy$pooled_p; names(allp)<-psy$outcome
q_all <- p.adjust(allp, "BH")
sig_tiered <- c("ASD","ADHD","F3")
cat(sprintf("  N psychiatric pooled tests (incl noUKBB + factors) = %d\n", length(allp)))
cat("  Single-pool BH q<0.05 hits:\n")
for (nm in names(q_all)[order(q_all)]) {
  if (q_all[nm] < 0.10) cat(sprintf("     %-10s p=%.4g q=%.4g %s\n", nm, allp[nm], q_all[nm],
      ifelse(q_all[nm]<0.05,"SIG","(0.05-0.10)")))
}

cat("\n\n################ PART E: ADHD FRAGILITY + MRlap null check ################\n\n")
ml <- fread("results/mrlap/MRLAP_summary.tsv")
cat("--- MRlap corrected p for SCZ/BD/MDD (claim: remain null) ---\n")
for (oc in c("SCZ","BD","MDD")) {
  sub <- ml[outcome==oc]
  cat(sprintf("  %s: corrected p range = %.3g to %.3g (n=%d pairs); all>0.05? %s\n",
      oc, min(sub$ivw_corrected_p), max(sub$ivw_corrected_p), nrow(sub),
      all(sub$ivw_corrected_p>0.05)))
}
cat("\n--- ADHD: naive pooled p=0.019. Under rho it crosses 0.05. Single-panel max sig: ---\n")
adhd <- fwd[outcome=="ADHD"]
cat(sprintf("  best single-panel ADHD p = %.4g (%s)\n", min(adhd$p), adhd$exposure[which.min(adhd$p)]))
cat(sprintf("  ADHD MRlap corrected p (Gupta_raw, the only nominal): %.4g\n",
    ml[outcome=="ADHD" & exposure=="Gupta_2023_raw"]$ivw_corrected_p))

cat("\n\n################ PART F: MRlap-observed vs forward-MR-IVW DISCREPANCY ################\n\n")
cat("Manuscript headline per-panel ASD betas come from forward MR (LD-clump r2<0.001 + Steiger + F>10).\n")
cat("MRlap 'observed' betas use a DIFFERENT instrument set (distance-prune 500kb, NO LD, NO Steiger).\n\n")
cat(sprintf("%-16s %12s %12s %12s %12s\n","panel","fwd_beta","fwd_nIV","mrlap_obs_b","mrlap_nIV"))
for (p in c("Longchamps_2022","Chong_2022","Gupta_2023_raw","Gupta_2023_adj")) {
  fb <- fwd[outcome=="ASD" & exposure==p]
  mb <- ml[outcome=="ASD" & exposure==p]
  cat(sprintf("%-16s %12.4f %12d %12.4f %12d\n", p, fb$beta, fb$n_ivs, mb$ivw_obs_b, mb$n_iv))
}
cat("\nRatio fwd/mrlap_obs for ASD ~ ", round(mean(fwd[outcome=="ASD"]$beta / ml[outcome=="ASD"]$ivw_obs_b),2), "\n")

cat("\nDONE2.\n")
