#!/usr/bin/env Rscript
# iter_003 — Step 01c: Build CHR_BP_A1_A2 -> rsID lookup table from 1KG EUR bim files
# Used to translate Gupta variant_id format -> rsID for PLINK clumping compatibility
# Date: 2026-05-22

suppressMessages({
  library(data.table)
  library(dplyr)
})

PROJECT <- Sys.getenv("ITER003_HOME", unset = getwd())
setwd(PROJECT)

LD_REF_DIR <- "/arf/home/ycicek/EVOSCZ/data/ldsc/sldsc_ref/1000G_EUR_Phase3_plink"
OUT_FILE   <- file.path(PROJECT, "data", "rsid_lookup_1KG_EUR.rds")

cat("Building rsID lookup from 1KG EUR bim files...\n")

# Read all 22 chr bim files
bim_list <- list()
for (chr in 1:22) {
  f <- file.path(LD_REF_DIR, sprintf("1000G.EUR.QC.%d.bim", chr))
  if (!file.exists(f)) {
    cat(sprintf("  chr %d: bim file missing %s\n", chr, f))
    next
  }
  # bim: chr rsID gen_pos bp a1 a2 (PLINK convention: a1=minor/alt, a2=major/ref)
  bim <- fread(f, col.names = c("CHR","rsID","GEN_POS","BP","A1","A2"))
  cat(sprintf("  chr %d: %d variants\n", chr, nrow(bim)))
  bim_list[[chr]] <- bim
}
all_bim <- rbindlist(bim_list)
cat(sprintf("\nTotal: %d variants in 1KG EUR\n", nrow(all_bim)))

# Build lookup keys in BOTH orderings (Gupta may have A1/A2 swapped vs 1KG)
# Key format: "chr_bp_X_Y"
all_bim[, key_a1a2 := sprintf("%s_%s_%s_%s", CHR, BP, A1, A2)]
all_bim[, key_a2a1 := sprintf("%s_%s_%s_%s", CHR, BP, A2, A1)]

# Long format: each row maps a key to rsID
lookup_a1a2 <- all_bim[, .(key = key_a1a2, rsID = rsID)]
lookup_a2a1 <- all_bim[, .(key = key_a2a1, rsID = rsID)]
lookup <- rbindlist(list(lookup_a1a2, lookup_a2a1))
lookup <- unique(lookup, by = "key")  # in case of conflicts, keep first

cat(sprintf("Lookup table: %d unique keys\n", nrow(lookup)))
cat(sprintf("Estimated memory: %.1f MB\n", as.numeric(object.size(lookup)) / 1024^2))

# Save as named character vector for fast lookup
lookup_vec <- setNames(lookup$rsID, lookup$key)
saveRDS(lookup_vec, OUT_FILE)
cat(sprintf("Saved: %s\n", OUT_FILE))
cat(sprintf("Size: %.1f MB on disk\n", file.info(OUT_FILE)$size / 1024^2))
