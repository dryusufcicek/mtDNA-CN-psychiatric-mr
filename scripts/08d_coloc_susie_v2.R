#!/usr/bin/env Rscript
# 08d_coloc_susie_v2.R — feasible SuSiE colocalization. The GTEx eQTL we have is
# significant-pairs only (too sparse for eQTL fine-mapping), so we fine-map the
# FULL GWAS pair (mtDNA-CN x disorder) with SuSiE (multi-signal robust) at a dense
# +/-250kb window, then test whether each candidate gene's LEAD eQTL variant is in
# the shared credible set (gene attribution; resolves 3p21 NDUFAF3 vs AMT).
suppressMessages({library(data.table); library(coloc); library(susieR)})
P<-"/arf/scratch/ycicek/iter_003_mtdna_psy_mr"; setwd(P)
PLINK<-"/arf/home/ycicek/bin/plink"; BF<-"/arf/scratch/ycicek/h2_paper2/data/g1000_eur/g1000_eur"
log<-function(m) cat(sprintf("[%s] %s\n",format(Sys.time(),"%H:%M:%S"),m))
DIS_N<-c(ADHD=225534,ASD=46350); DIS_S<-c(ADHD=0.1715,ASD=0.3966); DIS_FRQ<-c(ADHD="FRQ_U_186843",ASD=NA)
WIN<-250000
# candidate gene-tissue-disorder (incl BOTH 3p21 genes to compare attribution)
loci<-data.table(
 gene_id=c("ENSG00000160957.15","ENSG00000178057.15","ENSG00000145020.16","ENSG00000175806.15"),
 gene_name=c("RECQL4","NDUFAF3","AMT","MSRA"),
 tissue=c("Brain_Cortex","Brain_Cerebellar_Hemisphere","Brain_Amygdala","Brain_Hippocampus"),
 disorder=c("ADHD","ADHD","ADHD","ASD"))
OUTF<-list(ADHD="data/outcomes/ADHD_Demontis2023.meta.gz",ASD="data/outcomes/ASD_Grove2019.gz")
log("loading Longchamps + outcomes")
LC<-fread("data/exposures/Longchamps2022_mtDNA_CN.ALLm2.bgen.stats.gz",select=c("SNP","CHR","BP","ALLELE1","ALLELE0","BETA","SE","A1FREQ")); LC[,k:=paste0(CHR,":",BP)]
OUT<-list(); for(d in unique(loci$disorder)){
  o<-fread(OUTF[[d]], data.table=FALSE)   # plain data.frame: base-R $<- cannot silently drop a column
  o$BETA <- base::log(as.numeric(o$OR)); o$k <- paste0(o$CHR, ":", o$BP)   # base::log — local log() is the LOGGER (masking caused the 8-iter BETA=NULL bug)
  o$FRQ  <- if("FRQ_U_186843" %in% names(o)) as.numeric(o$FRQ_U_186843) else NA_real_
  o$A1 <- as.character(o$A1); o$A2 <- as.character(o$A2); o$SE <- as.numeric(o$SE)
  setDT(o); OUT[[d]] <- o
  log(sprintf("loaded %s: nrow=%d BETA_present=%s BETA_finite=%d OR_head=%s BETA_head=%s", d, nrow(o), "BETA"%in%names(o), sum(is.finite(o$BETA)), paste(round(head(as.numeric(o$OR),3),3),collapse=","), paste(round(head(o$BETA,3),3),collapse=",")))
}
palin<-function(a,b)(a=="A"&b=="T")|(a=="T"&b=="A")|(a=="C"&b=="G")|(a=="G"&b=="C")

run_one<-function(i){
 g<-loci$gene_id[i]; gn<-loci$gene_name[i]; ti<-loci$tissue[i]; d<-loci$disorder[i]
 log(sprintf("=== %s (%s) / %s ===",gn,ti,d))
 eq<-fread(sprintf("data/coloc_inputs/cis_eqtl_lifted_%s.tsv",ti)); e<-eq[gene_id==g & pos37>0]
 if(nrow(e)==0) return(NULL)
 lead<-e[which.min(pval_nominal)]; chrn<-as.integer(lead$chrn); ctr<-as.integer(lead$pos37); lead_k<-paste0(chrn,":",ctr)
 lo<-ctr-WIN; hi<-ctr+WIN
 reg<-sprintf("/tmp/rg_%s",gn)
 system(sprintf("%s --bfile %s --chr %d --from-bp %d --to-bp %d --make-bed --out %s --memory 6000 >/dev/null 2>&1",PLINK,BF,chrn,lo,hi,reg))
 if(!file.exists(paste0(reg,".bim"))) return(NULL)
 bim<-fread(paste0(reg,".bim")); setnames(bim,c("c","rsid","cm","bp","A1","A2")); bim[,k:=paste0(c,":",bp)]; bim<-bim[!palin(A1,A2)][!duplicated(k)]
 od<-OUT[[d]]; common<-Reduce(intersect,list(bim$k,od$k,LC$k))
 log(sprintf("  GWAS∩1000G common SNPs: %d",length(common))); if(length(common)<80) return(NULL)
 writeLines(bim[k %in% common]$rsid,paste0(reg,".keep"))
 system(sprintf("%s --bfile %s --extract %s --make-bed --out %s2 --memory 6000 >/dev/null 2>&1",PLINK,reg,paste0(reg,".keep"),reg))
 system(sprintf("%s --bfile %s2 --r square --out %sLD --memory 6000 >/dev/null 2>&1",PLINK,reg,reg))
 b3<-fread(paste0(reg,"2.bim")); setnames(b3,c("c","rsid","cm","bp","A1","A2")); b3[,k:=paste0(c,":",bp)]
 LD<-as.matrix(fread(paste0(reg,"LD.ld"))); ord<-b3$k; refA1<-setNames(b3$A1,ord); refA2<-setNames(b3$A2,ord)
 harm<-function(dt,EA,OA,BE,SEc,FRQ){
  m<-dt[k %in% ord][!duplicated(k)]; ix<-match(ord, m$k)
  be<-as.numeric(m[[BE]])[ix]; se<-as.numeric(m[[SEc]])[ix]
  ea<-toupper(as.character(m[[EA]]))[ix]; oa<-toupper(as.character(m[[OA]]))[ix]
  a1<-refA1[ord]; a2<-refA2[ord]
  good<-((ea==a1&oa==a2)|(ea==a2&oa==a1))&is.finite(be)&is.finite(se)&se>0; good[is.na(good)]<-FALSE
  bh<-be; fl<-!is.na(ea)&ea==a2; bh[fl]<- -be[fl]
  frq<-if(!is.na(FRQ)&&FRQ%in%names(m)) as.numeric(m[[FRQ]])[ix] else rep(0.3,length(ord))
  list(beta=bh,se=se,good=good,maf=pmin(frq,1-frq))}
 Hd<-harm(od,"A1","A2","BETA","SE","FRQ"); Hm<-harm(LC,"ALLELE1","ALLELE0","BETA","SE","A1FREQ")
 log(sprintf("  Hd good=%d / Hm good=%d (n_ord=%d)", sum(Hd$good,na.rm=TRUE), sum(Hm$good,na.rm=TRUE), length(ord)))
 keep<-Hd$good&Hm$good; keep[is.na(keep)]<-FALSE; nk<-sum(keep); log(sprintf("  usable SNPs: %d",nk)); if(nk<80) return(NULL)
 idx<-which(keep); LDk<-LD[idx,idx]; snps<-ord[idx]; dimnames(LDk)<-list(snps,snps); pos<-as.integer(sub(".*:","",snps))
 mf<-function(x){x[is.na(x)|x<=0|x>=1]<-0.1;x}
 Dd<-list(beta=Hd$beta[idx],varbeta=Hd$se[idx]^2,snp=snps,position=pos,type="cc",N=DIS_N[[d]],s=DIS_S[[d]],MAF=mf(Hd$maf[idx]),LD=LDk)
 Dm<-list(beta=Hm$beta[idx],varbeta=Hm$se[idx]^2,snp=snps,position=pos,type="quant",N=465809,MAF=mf(Hm$maf[idx]),LD=LDk)
 Sd<-tryCatch(runsusie(Dd,suffix="d"),error=function(z){log(paste("  susie dis fail:",conditionMessage(z)));NULL})
 Sm<-tryCatch(runsusie(Dm,suffix="m"),error=function(z){log(paste("  susie mt fail:",conditionMessage(z)));NULL})
 pp<-NA; shared_in<-NA
 if(!is.null(Sd)&&!is.null(Sm)){cs<-tryCatch(coloc.susie(Sm,Sd),error=function(z)NULL)
   if(!is.null(cs$summary)){pp<-max(cs$summary$PP.H4.abf)}}
 lead_in_mt<-if(!is.null(Sm)) lead_k %in% snps[unlist(Sm$sets$cs)] else NA
 lead_in_dis<-if(!is.null(Sd)) lead_k %in% snps[unlist(Sd$sets$cs)] else NA
 r<-data.table(gene=gn,tissue=ti,disorder=d,n_snps=nk,lead_eqtl=lead_k,
   cs_mtdna=if(!is.null(Sm))length(Sm$sets$cs)else NA, cs_disorder=if(!is.null(Sd))length(Sd$sets$cs)else NA,
   susie_PP4_mtdna_disorder=pp, lead_eqtl_in_mtdna_CS=lead_in_mt, lead_eqtl_in_disorder_CS=lead_in_dis)
 log(sprintf("  CS mt/dis=%s/%s | susie-coloc mtDNAxdis PP4=%s | leadEQTL in mtCS=%s disCS=%s",
   r$cs_mtdna,r$cs_disorder,round(pp,3),lead_in_mt,lead_in_dis))
 r}
out<-list(); for(i in seq_len(nrow(loci))){r<-tryCatch(run_one(i),error=function(z){log(paste("LOCUS FAIL:",conditionMessage(z)));NULL}); if(!is.null(r))out[[length(out)+1]]<-r}
if(length(out)){df<-rbindlist(out,fill=TRUE); fwrite(df,"results/coloc/COLOC_SUSIE_v2.tsv",sep="\t"); cat("\n=== coloc.susie (mtDNA x disorder) + lead-eQTL attribution ===\n"); print(as.data.frame(df))}
log("DONE")
