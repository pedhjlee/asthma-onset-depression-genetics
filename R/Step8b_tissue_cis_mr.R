###############################################################################
## Step 8b — 조직특이 cis-MR (GTEx 폐·뇌·전혈 eQTL → 천식·우울/불안)
## Project: 01_GWAS_Asthma_MDD  (면역축: 조직특이 검정)
##
## Step8a가 저장한 raw/immune/tissue_eqtl/<TISSUE>__<GENE>.tsv 를 읽어:
##   nes→beta, SE=|nes|/z(z=qnorm(p/2)), EA/NEA=variantId(alt/ref), palindrome 제거,
##   cis 도구 clump → cis-MR (결과 binary=OR).
## 핵심: 폐 eQTL→천식(COA/AOA), 뇌 eQTL→우울/불안(MDD/ANX) '조직매칭'.
##   특히 IL33(뇌)→우울, TSLP(폐)→천식 등 혈액서 못본 유전자 검정.
## 실행: source("Step8b_tissue_cis_mr.R")   (Step8a 먼저)
## 주: coloc(B)는 전영역 데이터 필요 → MR 신호 유전자에 한해 후속(GTEx allpairs).
###############################################################################

suppressMessages({ library(data.table); library(TwoSampleMR); library(ggplot2) })

PROJ <- "PATH/TO/PROJECT"
HARM <- file.path(PROJ,"harmonized"); MRD <- file.path(PROJ,"mr")
EQD  <- file.path(PROJ,"raw","immune","tissue_eqtl"); setwd(MRD)

OUTC <- c("COA","AOA","MDD","ANX")
MATCH <- list(Lung=c("COA","AOA"), Brain_Cortex=c("MDD","ANX"),
              Whole_Blood=c("COA","AOA","MDD","ANX"))   # 조직-질환 매칭
CLUMP_R2 <- 0.001; CLUMP_KB <- 1000; MIN_IV <- 1; set.seed(20260714)

## 로컬 clump
LOCAL_BFILE <- "PATH/TO/ldref_ascii/EUR"
{ ad<-"PATH/TO/ldref_ascii"; dir.create(ad,showWarnings=FALSE,recursive=TRUE)
  for(e in c("bed","bim","fam")){d<-file.path(ad,paste0("EUR.",e));s<-file.path(PROJ,"ld_ref",paste0("EUR.",e))
    if(!file.exists(d)&&file.exists(s)) file.copy(s,d)} }
if (!requireNamespace("genetics.binaRies",quietly=TRUE)) remotes::install_github("MRCIEU/genetics.binaRies")
LOCAL_PLINK <- genetics.binaRies::get_plink_binary()
USE_LOCAL <- file.exists(paste0(LOCAL_BFILE,".bed"))
do_clump <- function(dat){
  if (USE_LOCAL){ out<-try(ieugwasr::ld_clump(dplyr::tibble(rsid=dat$SNP,pval=dat$pval.exposure,id=dat$id.exposure),
        clump_kb=CLUMP_KB,clump_r2=CLUMP_R2,bfile=LOCAL_BFILE,plink_bin=LOCAL_PLINK),silent=TRUE)
    if(!inherits(out,"try-error")) return(dat[dat$SNP %in% out$rsid,]) }
  out<-try(clump_data(dat,clump_r2=CLUMP_R2,clump_kb=CLUMP_KB,pop="EUR"),silent=TRUE)
  if(inherits(out,"try-error")) return(dat); out
}

## GTEx tsv → exposure df (nes→beta/se, variantId 파싱, palindrome 제거)
build_exposure <- function(f, gene, tissue){
  d <- fread(f)
  vp <- tstrsplit(d$variantId, "_", fixed=TRUE)
  d[, `:=`(CHR=sub("chr","",vp[[1]]), POS=suppressWarnings(as.integer(vp[[2]])),
           NEA=vp[[3]], EA=vp[[4]])]
  d <- d[EA %in% c("A","C","G","T") & NEA %in% c("A","C","G","T") & EA!=NEA]
  d <- d[is.finite(pValue) & pValue>0 & is.finite(nes)]
  z <- qnorm(d$pValue/2, lower.tail=FALSE)                 # 양수 z
  d[, `:=`(BETA=nes, SE=abs(nes)/z)]
  d <- d[is.finite(SE) & SE>0]
  ## palindrome 제거(EAF 없어 strand 해결 불가)
  d <- d[!((EA=="A"&NEA=="T")|(EA=="T"&NEA=="A")|(EA=="C"&NEA=="G")|(EA=="G"&NEA=="C"))]
  if (!nrow(d)) return(NULL)
  e <- format_data(as.data.frame(d), type="exposure", snp_col="snpId", beta_col="BETA", se_col="SE",
        effect_allele_col="EA", other_allele_col="NEA", pval_col="pValue", chr_col="CHR", pos_col="POS")
  e$exposure <- paste0(gene,"_",tissue); e$id.exposure <- paste0(gene,"_",tissue)
  cl <- do_clump(e)
  message(sprintf("  %-6s x %-12s: %d SNP → clump %d", gene, tissue, nrow(d), nrow(cl)))
  attr(cl,"gene")<-gene; attr(cl,"tissue")<-tissue; cl
}

prep_outcome <- function(tr, snps){
  dt <- as.data.table(readRDS(file.path(HARM,paste0(tr,".rds"))))
  sub <- dt[SNP %in% snps]; if(!nrow(sub)) return(NULL)
  o <- format_data(as.data.frame(sub), type="outcome", snp_col="SNP", beta_col="BETA", se_col="SE",
        effect_allele_col="EA", other_allele_col="NEA", eaf_col="EAF", pval_col="P", samplesize_col="N")
  o$outcome<-tr; o$id.outcome<-tr; o
}
mr_one <- function(e, ou, gene, tissue){
  o <- prep_outcome(ou, e$SNP); if(is.null(o)) return(NULL)
  dat <- harmonise_data(e,o,action=2); dat<-dat[dat$mr_keep,]
  if(nrow(dat) < MIN_IV) return(NULL)
  Fm <- mean((dat$beta.exposure/dat$se.exposure)^2)
  res <- mr(dat, method_list=if(nrow(dat)==1) "mr_wald_ratio" else c("mr_ivw","mr_egger_regression","mr_weighted_median"))
  prim <- if(nrow(dat)==1) res[res$method=="Wald ratio",] else res[res$method=="Inverse variance weighted",]
  data.table(GENE=gene, TISSUE=tissue, OUTCOME=ou, N_IV=nrow(dat), F_MEAN=round(Fm,1),
    METHOD=prim$method, BETA=round(prim$b,4), SE=round(prim$se,4),
    OR=round(exp(prim$b),3), OR_LCI=round(exp(prim$b-1.96*prim$se),3),
    OR_UCI=round(exp(prim$b+1.96*prim$se),3), P=signif(prim$pval,3),
    MATCHED = ou %in% MATCH[[tissue]])
}

## ---- 실행 --------------------------------------------------------------- ##
files <- list.files(EQD, pattern="__.*\\.tsv$", full.names=TRUE)
stopifnot(length(files)>0)
message("=== 조직 cis-MR (", length(files), " gene x tissue 파일) ===")
RES <- rbindlist(lapply(files, function(f){
  bn <- sub("\\.tsv$","",basename(f)); parts <- strsplit(bn,"__")[[1]]
  tissue <- parts[1]; gene <- parts[2]
  e <- try(build_exposure(f, gene, tissue), silent=TRUE)
  if (inherits(e,"try-error")||is.null(e)||!nrow(e)) return(NULL)
  rbindlist(lapply(OUTC, function(ou) mr_one(e, ou, gene, tissue)), fill=TRUE)
}), fill=TRUE)

if (nrow(RES)) {
  ## FDR: 조직매칭 primary 가설 내에서
  RES[MATCHED==TRUE, P_FDR := signif(p.adjust(P,"BH"),3)]
  fwrite(RES, "mr_tissue_cis.csv")
  message("\n=== 조직 cis-MR 요약 (매칭 primary) ===")
  print(RES[MATCHED==TRUE][order(P), .(GENE,TISSUE,OUTCOME,N_IV,F_MEAN,METHOD,OR,OR_LCI,OR_UCI,P,P_FDR)])
  message("\n[참고] 비매칭(교차) 결과는 mr_tissue_cis.csv 의 MATCHED=FALSE 행")
  sig <- RES[MATCHED==TRUE & P<0.05]
  message("\n★ 조직매칭 유의(P<0.05): ",
          if(nrow(sig)) paste(sprintf("%s(%s)→%s OR%.2f", sig$GENE,sig$TISSUE,sig$OUTCOME,sig$OR), collapse=" · ") else "없음")

  ## 그림 (매칭)
  d <- RES[MATCHED==TRUE]; d[, lab:=paste0(GENE," [",TISSUE,"] → ",OUTCOME)]
  d[, sig:=!is.na(P_FDR)&P_FDR<0.05]; setorder(d,TISSUE,GENE,OUTCOME)
  d[, lab:=factor(lab,levels=rev(unique(lab)))]
  p <- ggplot(d, aes(OR,lab)) + geom_vline(xintercept=1,linetype=2,color="grey55") +
    geom_errorbarh(aes(xmin=OR_LCI,xmax=OR_UCI),height=0.2,color="grey35") +
    geom_point(aes(fill=sig),shape=21,size=2.6,color="grey20") +
    scale_fill_manual(values=c(`TRUE`="#B2182B",`FALSE`="white"),guide="none") +
    scale_x_log10() +
    labs(title="Tissue cis-MR: gene expression -> asthma/depression",
         subtitle="Lung->asthma, Brain->depression/anxiety (tissue-matched)",
         x="OR (95% CI) per SD expression", y=NULL) + theme_minimal(base_size=10)
  ggsave("fig_tissue_cis_mr.png", p, width=8, height=5.5, dpi=300, bg="white")
  message("  fig_tissue_cis_mr.png 저장")
}
message("\n================= Step 8b 완료 =================")
message("산출: mr_tissue_cis.csv, fig_tissue_cis_mr.png")
message("해석: 폐발현→천식 & 뇌발현→우울이 유의한 유전자 = 조직특이 공유 면역기전.")
message("다음(coloc): MR 유의 유전자만 GTEx allpairs로 전영역 받아 천식+우울+eQTL 3-way coloc.")
