suppressMessages({library(data.table); library(dplyr)})
setwd("/arf/scratch/ycicek/iter_003_mtdna_psy_mr")

all_res <- readRDS("results/forward_mr/all_forward_results.rds")
conc <- fread("results/method_comparison/PAIRWISE_CONCORDANCE.tsv")
rbar <- mean(conc$pearson_r)   # mean inter-panel correlation

# correlation-aware fixed-effect IVW pool under equicorrelation rho
corr_pool <- function(b, se, rho) {
  k <- length(b); w <- 1/se^2; a <- w/sum(w)
  Sigma <- outer(se, se) * rho; diag(Sigma) <- se^2
  pooled <- sum(a*b)
  var_p  <- as.numeric(t(a) %*% Sigma %*% a)
  se_p <- sqrt(var_p); z <- pooled/se_p; p <- 2*pnorm(-abs(z))
  c(beta=pooled, se=se_p, p=p)
}

primary <- c("SCZ","BD","MDD","ADHD","ASD","F3","ANX","OCD","PTSD","AN","BPD","TS","PPD",
             "F1","F2","F4","F5","PFactor","BMI","Height")
rows <- list()
for (out in primary) {
  pr <- list()
  for (nm in names(all_res)) {
    r <- all_res[[nm]]
    if (toupper(r$outcome)==toupper(out) && !is.null(r$main)) {
      iv <- r$main[r$main$method=="Inverse variance weighted",]
      if (nrow(iv)>0) pr[[r$exposure]] <- c(iv$b, iv$se)
    }
  }
  if (length(pr) < 3) next
  M <- do.call(rbind, pr); b <- M[,1]; se <- M[,2]
  naive <- corr_pool(b, se, 0)         # rho=0 reproduces metafor fixed
  c70 <- corr_pool(b, se, 0.70)
  c55 <- corr_pool(b, se, 0.55)
  c90 <- corr_pool(b, se, 0.90)
  rows[[out]] <- data.frame(
    outcome=out, n_panels=length(b),
    naive_beta=round(naive["beta"],4), naive_p=signif(naive["p"],3),
    corr70_beta=round(c70["beta"],4), corr70_se=round(c70["se"],4), corr70_p=signif(c70["p"],3),
    corr55_p=signif(c55["p"],3), corr90_p=signif(c90["p"],3),
    stringsAsFactors=FALSE)
}
df <- bind_rows(rows)
cat(sprintf("Mean inter-panel r = %.3f; DEFF(rho=0.70,k=4)=%.2f; k_eff=%.2f\n\n",
            rbar, 1+3*0.70, 4/(1+3*0.70)))
print(df, row.names=FALSE)
fwrite(df, "manuscript_draft/supplementary/Supp_Table_12_correlation_aware_meta.tsv", sep="\t")

# Primary-tier BH-FDR on corr70_p (5 primary disorders)
cat("\n--- Primary-tier BH-FDR on correlation-aware p (SCZ,BD,MDD,ADHD,ASD) ---\n")
prim5 <- df[df$outcome %in% c("SCZ","BD","MDD","ADHD","ASD"),]
prim5$q <- p.adjust(prim5$corr70_p, method="BH")
print(prim5[,c("outcome","naive_p","corr70_p","q")], row.names=FALSE)
