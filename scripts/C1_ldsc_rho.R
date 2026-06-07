#!/usr/bin/env Rscript
# C1_ldsc_rho.R (v2, fixed) — sample-overlap correlation among 4 mtDNA-CN panels
# via bivariate LDSC genetic-covariance INTERCEPTS (GenomicSEM). FIX: use d[[col]]
# (not get(col) inside data.table j, which collided with lowercase 'beta' column).
suppressMessages({library(data.table); library(GenomicSEM)})
PROJECT <- "/arf/scratch/ycicek/iter_003_mtdna_psy_mr"; setwd(PROJECT)
dir.create("data/ldsc_rho", showWarnings=FALSE); dir.create("results/ldsc_rho", showWarnings=FALSE)

BIM <- fread("/arf/scratch/ycicek/h2_paper2/data/g1000_eur/g1000_eur.bim", header=FALSE)
setnames(BIM, c("CHR","SNP","CM","BP","A1b","A2b"))
BIM[, keyp := paste(CHR, BP, sep=":")]
bimkey <- BIM[!duplicated(keyp), .(keyp, rsid=SNP)]
cat(sprintf("[info] bim chr:bp->rsid keys: %d\n", nrow(bimkey)))

standardize <- function(path, a1, a2, beta, se, p, N, snpcol=NULL, chrcol=NULL, bpcol=NULL){
  d <- fread(path)
  if (!is.null(snpcol)) {
    SNP <- as.character(d[[snpcol]])
  } else {
    d[, keyp := paste(d[[chrcol]], d[[bpcol]], sep=":")]
    d <- merge(d, bimkey, by="keyp", all.x=TRUE)
    SNP <- as.character(d[["rsid"]])
  }
  out <- data.table(SNP=SNP, A1=toupper(as.character(d[[a1]])), A2=toupper(as.character(d[[a2]])),
                    BETA=as.numeric(d[[beta]]), SE=as.numeric(d[[se]]), P=as.numeric(d[[p]]), N=N)
  out <- out[!is.na(SNP) & startsWith(SNP,"rs") & is.finite(BETA) & is.finite(SE) & is.finite(P) &
             P>0 & P<=1 & A1 %in% c("A","C","G","T") & A2 %in% c("A","C","G","T")]
  out[!duplicated(SNP)]
}

cat("[step] standardize 4 exposures\n")
Ls <- standardize("data/exposures/Longchamps2022_mtDNA_CN.ALLm2.bgen.stats.gz","ALLELE1","ALLELE0","BETA","SE","P_BOLT_LMM",465809, snpcol="SNP")
fwrite(Ls,"data/ldsc_rho/Longchamps.tsv",sep="\t"); cat(sprintf("  Longchamps: %d\n",nrow(Ls))); rm(Ls); gc()
Cs <- standardize("data/exposures/Chong2022_mtDNA_CN_EUR_GCST90026372.tsv.gz","effect_allele","other_allele","beta","standard_error","p_value",383476, snpcol="variant_id")
fwrite(Cs,"data/ldsc_rho/Chong.tsv",sep="\t"); cat(sprintf("  Chong: %d\n",nrow(Cs))); rm(Cs); gc()
GRs <- standardize("data/exposures/Gupta2023_mtDNA_CN_raw_GCST90268498.tsv.gz","effect_allele","other_allele","beta","standard_error","p_value",155998, chrcol="chromosome", bpcol="base_pair_location")
fwrite(GRs,"data/ldsc_rho/Gupta_raw.tsv",sep="\t"); cat(sprintf("  Gupta_raw: %d\n",nrow(GRs))); rm(GRs); gc()
GAs <- standardize("data/exposures/Gupta2023_mtDNA_CN_adjusted_GCST90268497.tsv.gz","effect_allele","other_allele","beta","standard_error","p_value",155998, chrcol="chromosome", bpcol="base_pair_location")
fwrite(GAs,"data/ldsc_rho/Gupta_adj.tsv",sep="\t"); cat(sprintf("  Gupta_adj: %d\n",nrow(GAs))); rm(GAs); gc()

cat("[step] munge\n")
setwd(file.path(PROJECT,"data/ldsc_rho"))
munge(files=c("Longchamps.tsv","Chong.tsv","Gupta_raw.tsv","Gupta_adj.tsv"), hm3="w_hm3.snplist",
      trait.names=c("Longchamps","Chong","Gupta_raw","Gupta_adj"), N=c(465809,383476,155998,155998))

cat("[step] bivariate LDSC\n")
ld <- file.path(PROJECT,"data/ldsc_mrlap/eur_w_ld_chr")
out <- ldsc(traits=c("Longchamps.sumstats.gz","Chong.sumstats.gz","Gupta_raw.sumstats.gz","Gupta_adj.sumstats.gz"),
            sample.prev=rep(NA,4), population.prev=rep(NA,4), ld=ld, wld=ld,
            trait.names=c("Longchamps","Chong","Gupta_raw","Gupta_adj"))
saveRDS(out, file.path(PROJECT,"results/ldsc_rho/LDSCoutput.rds"))
cat("\n=== Intercept matrix I (off-diag = sample-overlap intercept) ===\n"); print(round(out$I,4))
cat("\n=== Genetic correlation ===\n"); print(round(cov2cor(out$S),4))
I <- out$I; od <- I[upper.tri(I)]
cat(sprintf("\n>>> Mean off-diagonal LDSC intercept (sample-overlap rho): %.4f  [%.4f - %.4f]\n", mean(od), min(od), max(od)))
writeLines(c(sprintf("ldsc_rho_mean\t%.4f", mean(od)), sprintf("ldsc_rho_min\t%.4f", min(od)),
             sprintf("ldsc_rho_max\t%.4f", max(od)), "manuscript_rho\t0.70"),
           file.path(PROJECT,"results/ldsc_rho/rho_summary.tsv"))
cat("\n[done] C1\n")
