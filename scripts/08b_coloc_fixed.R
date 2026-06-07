#!/usr/bin/env Rscript
# 08b_coloc_fixed.R — brain cis-eQTL x (ASD/ADHD/F3) colocalization with the
# liftOver fix: match eQTL (now carrying GRCh37 pos37) to outcome GWAS by chr:pos
# (GRCh37). Flags MitoCarta3.0 genes. Replaces the failed rsid-match in 08.
suppressMessages({library(data.table); library(dplyr); library(coloc)})
PROJECT <- "/arf/scratch/ycicek/iter_003_mtdna_psy_mr"; setwd(PROJECT)
log_step <- function(m) cat(sprintf("\n=== %s ===\n", m))
log_info <- function(m) cat(sprintf("[info] %s\n", m))
log_warn <- function(m) cat(sprintf("[warn] %s\n", m))
OUT_DIR <- file.path(PROJECT,"results","coloc"); COLOC_IN <- file.path(PROJECT,"data","coloc_inputs")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

OUTCOME_META <- list(
  ASD  = list(file="data/outcomes/ASD_Grove2019.gz",          n_case=18381, n_control=27969, type="cc", s=18381/(18381+27969)),
  ADHD = list(file="data/outcomes/ADHD_Demontis2023.meta.gz", n_case=38691, n_control=186843, type="cc", s=38691/(38691+186843)),
  F3   = list(file="data/outcomes/CDG3_F3_Neurodev.tsv.gz",   type="quant", N=84760))
BRAIN_TISSUES <- c("Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia",
 "Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9",
 "Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia",
 "Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra")
GTEX_N <- c(Brain_Amygdala=152,Brain_Anterior_cingulate_cortex_BA24=176,Brain_Caudate_basal_ganglia=246,
 Brain_Cerebellar_Hemisphere=215,Brain_Cerebellum=250,Brain_Cortex=255,Brain_Frontal_Cortex_BA9=209,
 Brain_Hippocampus=197,Brain_Hypothalamus=202,Brain_Nucleus_accumbens_basal_ganglia=246,
 Brain_Putamen_basal_ganglia=205,`Brain_Spinal_cord_cervical_c-1`=159,Brain_Substantia_nigra=140)

log_step("read outcomes (key = CHR:BP GRCh37)")
read_outcome <- function(meta){
  d <- fread(meta$file)
  if ("OR" %in% names(d) && !"BETA" %in% names(d)) d[, BETA := log(OR)]
  if (!"BP" %in% names(d)){ pc <- intersect(c("POS","BP","POSITION","base_pair_location"),names(d))[1]; if(!is.na(pc)) setnames(d,pc,"BP")}
  if (!"CHR" %in% names(d)){ cc <- intersect(c("CHROM","Chromosome","chr","chromosome"),names(d))[1]; if(!is.na(cc)) setnames(d,cc,"CHR")}
  d[, CHR := as.integer(gsub("chr","",as.character(CHR)))]
  d[, key := paste0(CHR, ":", BP)]
  d[]
}
outcome_data <- lapply(OUTCOME_META, read_outcome)
for (o in names(outcome_data)) log_info(sprintf("%s: %d SNPs", o, nrow(outcome_data[[o]])))

log_step("coloc.abf over tissue x gene x outcome (chr:pos37 match)")
res_list <- list(); si <- 0L
for (tissue in BRAIN_TISSUES){
  fp <- file.path(COLOC_IN, sprintf("cis_eqtl_lifted_%s.tsv", tissue))
  if (!file.exists(fp)){ log_warn(sprintf("missing %s",tissue)); next }
  eq <- fread(fp); if (nrow(eq)==0) next
  eq[, key := paste0(chrn, ":", pos37)]
  n_eqtl <- GTEX_N[[tissue]]; if (is.null(n_eqtl)||is.na(n_eqtl)) n_eqtl <- 200
  log_info(sprintf("[%s] %d rows, %d genes", tissue, nrow(eq), uniqueN(eq$gene_id)))
  for (gene in unique(eq$gene_id)){
    e <- eq[gene_id==gene]
    e <- e[!is.na(slope)&!is.na(slope_se)&slope_se>0&!is.na(af)&af>0&af<1&pos37>0]
    if (nrow(e)<30) next
    setorder(e, pval_nominal); e <- e[!duplicated(key)]
    d1 <- list(type="quant", beta=e$slope, varbeta=e$slope_se^2, snp=e$key, MAF=pmin(e$af,1-e$af), N=n_eqtl)
    for (oname in names(OUTCOME_META)){
      meta <- OUTCOME_META[[oname]]; om <- outcome_data[[oname]]
      og <- om[key %in% d1$snp]; if (nrow(og)<30) next
      og <- og[!duplicated(key)]
      if (meta$type=="cc"){
        d2 <- list(type="cc", beta=og$BETA, varbeta=og$SE^2, snp=og$key, s=meta$s, N=meta$n_case+meta$n_control)
        mc <- intersect(c("FRQ_U_186843","FRQ_U"), names(og)); if (length(mc)) d2$MAF <- og[[mc[1]]]
      } else {
        d2 <- list(type="quant", beta=og$BETA, varbeta=og$SE^2, snp=og$key, N=meta$N)
        if ("MAF" %in% names(og)) d2$MAF <- og$MAF
      }
      keep <- !is.na(d2$beta) & !is.na(d2$varbeta) & d2$varbeta>0
      d2 <- lapply(d2, function(x) if(length(x)<=1) x else x[keep])
      common <- intersect(d1$snp, d2$snp); if (length(common)<30) next
      i1 <- match(common,d1$snp); i2 <- match(common,d2$snp)
      d1b <- lapply(d1, function(x) if(length(x)<=1) x else x[i1])
      d2b <- lapply(d2, function(x) if(length(x)<=1) x else x[i2])
      cr <- tryCatch(coloc::coloc.abf(d1b,d2b), error=function(ex){log_warn(sprintf("fail %s/%s/%s: %s",tissue,gene,oname,conditionMessage(ex)));NULL})
      if (is.null(cr)) next
      si <- si+1L
      res_list[[si]] <- data.frame(tissue=tissue,gene_id=gene,outcome=oname,n_snps=length(common),
        PP0=cr$summary[["PP.H0.abf"]],PP1=cr$summary[["PP.H1.abf"]],PP2=cr$summary[["PP.H2.abf"]],
        PP3=cr$summary[["PP.H3.abf"]],PP4=cr$summary[["PP.H4.abf"]],stringsAsFactors=FALSE)
    }
  }
}
if (length(res_list)==0){ log_warn("no coloc results"); quit(status=1) }
cdf <- bind_rows(res_list)
fwrite(cdf, file.path(OUT_DIR,"COLOC_per_pair_FIXED.tsv"), sep="\t")
log_info(sprintf("wrote %d coloc tests", nrow(cdf)))

# gene symbol + MitoCarta flag
eg <- tryCatch(fread("/arf/scratch/ycicek/h2_paper2/data/gtex_v10_eqtl_brain/Brain_Cortex.v10.eGenes.txt.gz", select=c("gene_id","gene_name")), error=function(e) data.table(gene_id=character(),gene_name=character()))
cdf <- left_join(cdf, distinct(eg), by="gene_id")
gmt <- tryCatch(readLines("data/magma_inputs/mitocarta3.gmt"), error=function(e) character(0))
mito <- unique(unlist(lapply(strsplit(gmt,"\t"), function(x) if(length(x)>2) x[3:length(x)] else character(0))))
cdf$is_mito <- !is.na(cdf$gene_name) & toupper(cdf$gene_name) %in% toupper(mito)
fwrite(cdf %>% filter(PP4>=0.75) %>% arrange(desc(PP4)), file.path(OUT_DIR,"COLOC_strong_PP4_075_FIXED.tsv"), sep="\t")
fwrite(cdf %>% filter(PP4>=0.5,PP4<0.75) %>% arrange(desc(PP4)), file.path(OUT_DIR,"COLOC_suggestive_FIXED.tsv"), sep="\t")
fwrite(cdf %>% filter(is_mito,PP4>=0.5) %>% arrange(desc(PP4)), file.path(OUT_DIR,"COLOC_MITO_hits_FIXED.tsv"), sep="\t")
cat("\n=== per-outcome summary ===\n")
print(cdf %>% group_by(outcome) %>% summarise(n_tests=n(),max_pp4=round(max(PP4),3),n_strong=sum(PP4>=0.75),n_sugg=sum(PP4>=0.5&PP4<0.75),.groups="drop"))
cat("\n=== MITO-gene coloc PP4>=0.5 ===\n"); print(head(as.data.frame(cdf %>% filter(is_mito,PP4>=0.5) %>% arrange(desc(PP4)) %>% select(outcome,tissue,gene_name,gene_id,PP4,PP3,n_snps)),20))
cat("\n=== top strong overall ===\n"); print(head(as.data.frame(cdf %>% filter(PP4>=0.75) %>% arrange(desc(PP4)) %>% select(outcome,tissue,gene_name,is_mito,PP4,n_snps)),15))
log_step("DONE 08b_coloc_fixed.R")
