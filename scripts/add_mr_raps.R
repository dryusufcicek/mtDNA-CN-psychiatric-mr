# Efficient MR-RAPS-only patch:
# read all_forward_results.rds → for each pair, run mr_raps on $dat → append to $main
suppressMessages({
  library(TwoSampleMR)
  library(mr.raps)
  library(data.table)
  library(dplyr)
})

PROJECT <- "/arf/scratch/ycicek/iter_003_mtdna_psy_mr"
setwd(PROJECT)

cat(sprintf("[%s] Loading existing all_forward_results.rds\n", format(Sys.time())))
all_res <- readRDS("results/forward_mr/all_forward_results.rds")
cat(sprintf("  N pairs: %d\n", length(all_res)))

# Per-pair add MR-RAPS
new_rows_list <- list()
n_added <- 0
n_skip <- 0
n_fail <- 0
t0 <- Sys.time()

for (nm in names(all_res)) {
  r <- all_res[[nm]]
  if (is.null(r$dat) || is.null(r$main)) { n_skip <- n_skip + 1; next }
  if ("Robust adjusted profile score (RAPS)" %in% r$main$method) {
    n_skip <- n_skip + 1; next
  }
  raps_res <- tryCatch(
    mr(r$dat, method_list = c("mr_raps")),
    error = function(e) NULL
  )
  if (!is.null(raps_res) && nrow(raps_res) > 0) {
    # Match column names of existing main
    common_cols <- intersect(colnames(r$main), colnames(raps_res))
    raps_aligned <- raps_res[, common_cols, drop = FALSE]
    r$main <- rbind(r$main, raps_aligned)
    all_res[[nm]] <- r
    n_added <- n_added + 1
    if (n_added %% 10 == 0) {
      cat(sprintf("  [%s] %d pairs added (%.1f s elapsed)\n",
                  format(Sys.time()), n_added,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
  } else {
    n_fail <- n_fail + 1
    cat(sprintf("  [%s] FAIL: %s\n", format(Sys.time()), nm))
  }
}

cat(sprintf("\n[%s] DONE: %d added, %d skipped, %d failed\n",
            format(Sys.time()), n_added, n_skip, n_fail))

# Save
saveRDS(all_res, "results/forward_mr/all_forward_results.rds")
cat("Saved updated all_forward_results.rds\n")

# Rebuild SUMMARY_forward_mr.tsv to include MR-RAPS row per pair
summary_rows <- list()
for (nm in names(all_res)) {
  r <- all_res[[nm]]
  if (is.null(r$main) || nrow(r$main) == 0) next
  for (i in seq_len(nrow(r$main))) {
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      pair = r$pair_id,
      exposure = r$exposure,
      outcome = r$outcome,
      status = r$status,
      n_ivs = r$n_ivs,
      method = r$main$method[i],
      beta = r$main$b[i],
      se = r$main$se[i],
      p = r$main$pval[i],
      stringsAsFactors = FALSE
    )
  }
}
sum_df <- bind_rows(summary_rows)
cat(sprintf("Long summary: %d rows × %d cols\n", nrow(sum_df), ncol(sum_df)))
fwrite(sum_df, "results/forward_mr/SUMMARY_forward_mr_BY_METHOD.tsv", sep = "\t")
cat("Wrote SUMMARY_forward_mr_BY_METHOD.tsv (long format with method column)\n")
