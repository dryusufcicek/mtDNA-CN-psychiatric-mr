#!/usr/bin/env Rscript
# 08c_threeway_coloc.R — the MISSING leg for mediation: does each candidate
# gene's brain cis-eQTL colocalize with mtDNA-CN itself? Combined with the
# existing eQTL x disorder coloc (Supp 14), a shared causal variant across
# eQTL + mtDNA-CN + disorder is the three-way mediation chain.
suppressMessages({library(data.table); library(coloc); library(dplyr)})
P <- "/arf/scratch/ycicek/iter_003_mtdna_psy_mr"; setwd(P)
log <- function(m) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), m))

GTEX_N <- c(Brain_Amygdala=152,Brain_Anterior_cingulate_cortex_BA24=176,Brain_Caudate_basal_ganglia=246,
 Brain_Cerebellar_Hemisphere=215,Brain_Cerebellum=250,Brain_Cortex=255,Brain_Frontal_Cortex_BA9=209,
 Brain_Hippocampus=197,Brain_Hypothalamus=202,Brain_Nucleus_accumbens_basal_ganglia=246,
 Brain_Putamen_basal_ganglia=205,`Brain_Spinal_cord_cervical_c-1`=159,Brain_Substantia_nigra=140)
EXP <- list(
 Longchamps=list(f="data/exposures/Longchamps2022_mtDNA_CN.ALLm2.bgen.stats.gz", chr="CHR", bp="BP", b="BETA", se="SE", N=465809),
 Chong     =list(f="data/exposures/Chong2022_mtDNA_CN_EUR_GCST90026372.tsv.gz", chr="chromosome", bp="base_pair_location", b="beta", se="standard_error", N=383476),
 Gupta_raw =list(f="data/exposures/Gupta2023_mtDNA_CN_raw_GCST90268498.tsv.gz", chr="chromosome", bp="base_pair_location", b="beta", se="standard_error", N=155998),
 Gupta_adj =list(f="data/exposures/Gupta2023_mtDNA_CN_adjusted_GCST90268497.tsv.gz", chr="chromosome", bp="base_pair_location", b="beta", se="standard_error", N=155998))

# candidate gene x tissue x disorder hits (the 7 MitoCarta genes, PP4>=0.5)
mito <- fread("results/coloc/COLOC_MITO_hits_FIXED.tsv")
hits <- unique(mito[, .(gene_id, gene_name, tissue, outcome, PP4_disorder=PP4)])
log(sprintf("candidate hits: %d (genes: %s)", nrow(hits), paste(unique(hits$gene_name),collapse=",")))

# Pre-load the eQTL dataset per (gene,tissue) used by the hits
eqtl_cache <- list()
for (ti in unique(hits$tissue)) {
  fp <- sprintf("data/coloc_inputs/cis_eqtl_lifted_%s.tsv", ti)
  if (!file.exists(fp)) next
  eq <- fread(fp)
  eq[, key := paste0(chrn, ":", pos37)]
  for (g in unique(hits[tissue==ti]$gene_id)) {
    e <- eq[gene_id==g & !is.na(slope) & slope_se>0 & af>0 & af<1 & pos37>0]
    setorder(e, pval_nominal); e <- e[!duplicated(key)]
    if (nrow(e)>=30) eqtl_cache[[paste(g,ti)]] <- e
  }
}
log(sprintf("eQTL datasets cached: %d", length(eqtl_cache)))

# region per gene-tissue (chr + pos37 range)
res <- hits[, paste0("PP4_mtdna_", names(EXP)) := NA_real_]
res[, n_snps_mtdna := NA_integer_]

for (en in names(EXP)) {
  ex <- EXP[[en]]
  log(sprintf("loading exposure %s", en))
  d <- fread(ex$f, select=c(ex$chr, ex$bp, ex$b, ex$se))
  setnames(d, c(ex$chr, ex$bp, ex$b, ex$se), c("CHR","BP","BETA","SE"))
  d <- d[!is.na(BETA) & !is.na(SE) & SE>0]
  d[, CHR := suppressWarnings(as.integer(gsub("chr","",as.character(CHR))))]
  d[, key := paste0(CHR, ":", BP)]
  for (i in seq_len(nrow(res))) {
    ck <- paste(res$gene_id[i], res$tissue[i])
    e <- eqtl_cache[[ck]]; if (is.null(e)) next
    chrn <- as.integer(e$chrn[1]); lo <- min(e$pos37); hi <- max(e$pos37)
    sub <- d[CHR==chrn & BP>=lo & BP<=hi]
    if (nrow(sub) < 30) next
    sub <- sub[!duplicated(key)]
    common <- intersect(e$key, sub$key); if (length(common) < 30) next
    e2 <- e[match(common, key)]; s2 <- sub[match(common, key)]
    maf <- pmin(e2$af, 1-e2$af)
    d1 <- list(type="quant", beta=e2$slope, varbeta=e2$slope_se^2, snp=common, MAF=maf, N=GTEX_N[[res$tissue[i]]])
    d2 <- list(type="quant", beta=s2$BETA, varbeta=s2$SE^2, snp=common, MAF=maf, N=ex$N)
    cr <- tryCatch(coloc::coloc.abf(d1,d2), error=function(z) NULL)
    if (!is.null(cr)) {
      set(res, i, paste0("PP4_mtdna_",en), as.numeric(cr$summary[["PP.H4.abf"]]))
      if (en==names(EXP)[1] || is.na(res$n_snps_mtdna[i])) set(res, i, "n_snps_mtdna", length(common))
    }
  }
  rm(d); gc()
}

pp <- grep("PP4_mtdna_", names(res), value=TRUE)
res[, max_PP4_mtdna := apply(.SD, 1, function(x) if(all(is.na(x))) NA_real_ else max(x,na.rm=TRUE)), .SDcols=pp]
res[, threeway := !is.na(PP4_disorder) & PP4_disorder>=0.5 & !is.na(max_PP4_mtdna) & max_PP4_mtdna>=0.5]
setorder(res, -threeway, -PP4_disorder)
fwrite(res, "results/coloc/THREEWAY_eqtl_mtdna_disorder.tsv", sep="\t")

cat("\n=== THREE-WAY coloc: eQTL x mtDNA-CN x disorder ===\n")
print(as.data.frame(res[, .(gene_name, tissue, outcome, PP4_disorder=round(PP4_disorder,3),
      max_PP4_mtdna=round(max_PP4_mtdna,3), threeway, n_snps_mtdna)]))
cat(sprintf("\n>>> THREE-WAY hits (eQTL+mtDNA+disorder all PP4>=0.5): %d\n", sum(res$threeway, na.rm=TRUE)))
cat("   by gene:\n"); print(res[threeway==TRUE, .N, by=gene_name])
log("DONE 08c_threeway_coloc.R")
