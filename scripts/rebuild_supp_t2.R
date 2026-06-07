suppressMessages({library(data.table); library(dplyr); library(tidyr)})
PROJECT <- "/arf/scratch/ycicek/iter_003_mtdna_psy_mr"
setwd(PROJECT)
all_res <- readRDS("results/forward_mr/all_forward_results.rds")

safe_num <- function(x) suppressWarnings(as.numeric(x))

rows <- list()
for (nm in names(all_res)) {
  r <- all_res[[nm]]
  if (is.null(r$main)) next
  egger_int <- if (!is.null(r$pleio) && nrow(r$pleio) > 0) safe_num(r$pleio$egger_intercept[1]) else NA_real_
  egger_intp <- if (!is.null(r$pleio) && nrow(r$pleio) > 0) safe_num(r$pleio$pval[1]) else NA_real_
  q_p_ivw <- if (!is.null(r$heterogen)) {
    ivw_row <- r$heterogen[r$heterogen$method == "Inverse variance weighted", ]
    if (nrow(ivw_row) > 0) safe_num(ivw_row$Q_pval[1]) else NA_real_
  } else NA_real_
  presso_glob_p <- safe_num(if (!is.null(r$presso)) tryCatch(r$presso$`MR-PRESSO results`$`Global Test`$Pvalue, error=function(e) NA) else NA)
  presso_raw_p <- safe_num(if (!is.null(r$presso)) tryCatch(r$presso$`Main MR results`$`P-value`[r$presso$`Main MR results`$`MR Analysis` == "Raw"], error=function(e) NA) else NA)
  presso_corr_b <- safe_num(if (!is.null(r$presso)) tryCatch(r$presso$`Main MR results`$`Causal Estimate`[r$presso$`Main MR results`$`MR Analysis` == "Outlier-corrected"], error=function(e) NA) else NA)
  presso_corr_p <- safe_num(if (!is.null(r$presso)) tryCatch(r$presso$`Main MR results`$`P-value`[r$presso$`Main MR results`$`MR Analysis` == "Outlier-corrected"], error=function(e) NA) else NA)
  if (length(presso_glob_p) == 0 || is.null(presso_glob_p)) presso_glob_p <- NA_real_
  if (length(presso_corr_b) == 0) presso_corr_b <- NA_real_
  if (length(presso_corr_p) == 0) presso_corr_p <- NA_real_
  if (length(presso_raw_p) == 0) presso_raw_p <- NA_real_
  for (i in seq_len(nrow(r$main))) {
    rows[[length(rows)+1]] <- data.frame(
      pair=as.character(r$pair_id), exposure=as.character(r$exposure), outcome=as.character(r$outcome), n_iv=as.integer(r$n_ivs),
      method=as.character(r$main$method[i]), beta=safe_num(r$main$b[i]), se=safe_num(r$main$se[i]), p=safe_num(r$main$pval[i]),
      egger_intercept=egger_int, egger_intercept_p=egger_intp,
      cochran_Q_p_ivw=q_p_ivw,
      presso_global_p=presso_glob_p, presso_raw_p=presso_raw_p,
      presso_outlier_corrected_b=presso_corr_b, presso_outlier_corrected_p=presso_corr_p,
      stringsAsFactors=FALSE
    )
  }
}
df_long <- bind_rows(rows)
cat(sprintf("Long format: %d rows × %d cols\n", nrow(df_long), ncol(df_long)))
fwrite(df_long, "manuscript_draft/supplementary/Supp_Table_2_forward_MR_per_pair.tsv", sep="\t")

# Wide format
df_wide <- df_long %>%
  select(pair, exposure, outcome, n_iv, method, beta, se, p) %>%
  pivot_wider(names_from=method, values_from=c(beta, se, p), names_glue="{method}_{.value}") %>%
  left_join(df_long %>% filter(method == "Inverse variance weighted") %>%
            select(pair, egger_intercept, egger_intercept_p, cochran_Q_p_ivw, presso_global_p, presso_outlier_corrected_b, presso_outlier_corrected_p),
            by="pair")
cat(sprintf("Wide format: %d rows × %d cols\n", nrow(df_wide), ncol(df_wide)))
fwrite(df_wide, "manuscript_draft/supplementary/Supp_Table_2_forward_MR_per_pair_wide.tsv", sep="\t")

# Subset for ASD/ADHD/F3 protective pairs preview
prot <- df_long %>% filter(outcome %in% c("ASD","ADHD","F3")) %>%
  select(exposure, outcome, method, beta, se, p, n_iv) %>%
  arrange(outcome, exposure, method)
cat("\n=== Protective outcomes — 4 methods × 3 outcomes × 4 exposures = 48 rows ===\n")
cat(sprintf("RAPS subset: %d rows\n", sum(prot$method == "Robust adjusted profile score (RAPS)")))
cat("\n=== Markdown header (paste-ready) ===\n")
cat("| Exposure | Outcome | Method | β | SE | p | n_IV |\n")
cat("|---|---|---|---|---|---|---|\n")
for (i in 1:min(20, nrow(prot))) {
  cat(sprintf("| %s | %s | %s | %.3f | %.3f | %.3g | %d |\n",
              prot$exposure[i], prot$outcome[i], prot$method[i],
              prot$beta[i], prot$se[i], prot$p[i], prot$n_iv[i]))
}
