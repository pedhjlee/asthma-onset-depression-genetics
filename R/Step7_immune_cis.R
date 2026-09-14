###############################################################################
## Step 7 — 면역축 MR
## Project: 01_GWAS_Asthma_MDD  (공유 면역 인과축 검정)
##
## 메인(PART A): 공유 면역유전자 cis-eQTL(eQTLGen 전혈) → 천식(COA/AOA)·우울/불안(MDD/ANX)
##   "우리가 지목한 면역유전자의 발현이 천식·우울에 동시에 인과하나?" (cis-MR)
## 보조(PART B): 호산구(EOS)·CRP polygenic → 같은 4개 질환 (넓은 염증축·양성대조)
##
## 결과형질 4개 모두 binary → OR(SD 발현/노출당) 보고.
## eQTLGen: beta/SE 없음 → AF파일로 Z→beta/SE 변환(Zhu 2016).
##   beta = Z / sqrt(2p(1-p)(N+Z^2)),  SE = 1 / sqrt(2p(1-p)(N+Z^2)),  p=AssessedAllele 빈도
## 인프라: 로컬 plink clump(step4/5와 동일). 순수 R.
## 실행: source("Step7_immune_cis.R")
###############################################################################

suppressMessages({ library(data.table); library(TwoSampleMR); library(ggplot2) })

## --------------------------- 경로 --------------------------------------- ##
PROJ <- "PATH/TO/PROJECT"
RAW  <- file.path(PROJ,"raw"); HARM <- file.path(PROJ,"harmonized")
MRD  <- file.path(PROJ,"mr");  dir.create(MRD, showWarnings=FALSE)
setwd(MRD)

CIS_F <- file.path(RAW,"2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
AF_F  <- file.path(RAW,"2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz")
stopifnot(file.exists(CIS_F), file.exists(AF_F))

## --------------------------- 설정 --------------------------------------- ##
GENES  <- c("IL33","IL1RL1","IL4R","IL13","TSLP","TYK2","IL2RA")   # 공유 Th2/면역
OUTC   <- c("COA","AOA","MDD","ANX")                               # 결과(모두 binary)
CLUMP_R2 <- 0.001; CLUMP_KB <- 1000        # cis 영역 내 독립 도구
MIN_IV <- 1; SEED <- 20260714; set.seed(SEED)
IV_PVAL_POLY <- 5e-8                       # PART B polygenic 도구 임계

## --------------------- 로컬 clumping (ASCII 경로) ----------------------- ##
LOCAL_BFILE <- "PATH/TO/ldref_ascii/EUR"; LOCAL_PLINK <- NULL
{ ad <- "PATH/TO/ldref_ascii"; dir.create(ad, showWarnings=FALSE, recursive=TRUE)
  for (e in c("bed","bim","fam")) { d<-file.path(ad,paste0("EUR.",e)); s<-file.path(PROJ,"ld_ref",paste0("EUR.",e))
    if (!file.exists(d) && file.exists(s)) file.copy(s,d) } }
if (!requireNamespace("genetics.binaRies", quietly=TRUE)) remotes::install_github("MRCIEU/genetics.binaRies")
LOCAL_PLINK <- genetics.binaRies::get_plink_binary()
USE_LOCAL <- file.exists(paste0(LOCAL_BFILE,".bed"))
message("로컬 clumping: ", if (USE_LOCAL) "준비됨" else "★없음 → API 폴백")

do_clump <- function(dat) {
  if (USE_LOCAL) {
    out <- try(ieugwasr::ld_clump(dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure),
                 clump_kb=CLUMP_KB, clump_r2=CLUMP_R2, bfile=LOCAL_BFILE, plink_bin=LOCAL_PLINK), silent=TRUE)
    if (!inherits(out,"try-error")) return(dat[dat$SNP %in% out$rsid,])
  }
  out <- try(clump_data(dat, clump_r2=CLUMP_R2, clump_kb=CLUMP_KB, pop="EUR"), silent=TRUE)
  if (inherits(out,"try-error")) return(dat)   # clump 실패 시 원본(단일 cis면 그대로)
  out
}

## --------------------- eQTLGen 로드 & Z→beta 변환 ----------------------- ##
message("\n=== eQTLGen cis-eQTL 로드 (7유전자) ===")
cis <- fread(CIS_F)                                   # 전체 유의 cis
cis <- cis[GeneSymbol %in% GENES]
message("  대상 유전자 cis-eQTL 행: ", nrow(cis), " | 유전자별: ",
        paste(sort(unique(cis$GeneSymbol)), collapse=", "))
af  <- fread(AF_F, select=c("SNP","AlleleA","AlleleB","AlleleB_all"))
cis <- merge(cis, af, by="SNP")                       # 빈도 병합
rm(af); gc(FALSE)
## AssessedAllele 빈도 p
cis[, eaf := fifelse(AssessedAllele==AlleleB, AlleleB_all,
             fifelse(AssessedAllele==AlleleA, 1-AlleleB_all, NA_real_))]
cis <- cis[is.finite(eaf) & eaf>0 & eaf<1 & is.finite(Zscore) & NrSamples>0]
cis[, den := 2*eaf*(1-eaf)*(NrSamples + Zscore^2)]
cis[, `:=`(BETA = Zscore/sqrt(den), SE = 1/sqrt(den), P = Pvalue)]
cis[, F := (BETA/SE)^2]
message("  변환 후 유전자별 SNP수 / 최대 F:")
print(cis[, .(nSNP=.N, maxF=round(max(F))), by=GeneSymbol][order(-maxF)])

## --------------------- 공통 MR (exposure df → outcome) ------------------ ##
prep_outcome <- function(tr, snps) {
  dt <- as.data.table(readRDS(file.path(HARM, paste0(tr,".rds"))))
  sub <- dt[SNP %in% snps]; if (!nrow(sub)) return(NULL)
  o <- format_data(as.data.frame(sub), type="outcome", snp_col="SNP", beta_col="BETA", se_col="SE",
        effect_allele_col="EA", other_allele_col="NEA", eaf_col="EAF", pval_col="P", samplesize_col="N")
  o$outcome <- tr; o$id.outcome <- tr; o
}
## exposure df(이미 format된) × outcome → 요약 (결과 binary=OR)
mr_one <- function(e, ou, exp_label) {
  o <- prep_outcome(ou, e$SNP); if (is.null(o)) return(NULL)
  dat <- harmonise_data(e, o, action=2); dat <- dat[dat$mr_keep,]
  if (nrow(dat) < MIN_IV) return(NULL)
  Fm <- mean((dat$beta.exposure/dat$se.exposure)^2)
  res <- mr(dat, method_list = if (nrow(dat)==1) "mr_wald_ratio" else
              c("mr_ivw","mr_egger_regression","mr_weighted_median"))
  prim <- if (nrow(dat)==1) res[res$method=="Wald ratio",] else res[res$method=="Inverse variance weighted",]
  st <- try(directionality_test(dat), silent=TRUE)
  data.table(EXPOSURE=exp_label, OUTCOME=ou, N_IV=nrow(dat), F_MEAN=round(Fm,1),
    METHOD=prim$method, BETA=round(prim$b,4), SE=round(prim$se,4),
    OR=round(exp(prim$b),3), OR_LCI=round(exp(prim$b-1.96*prim$se),3),
    OR_UCI=round(exp(prim$b+1.96*prim$se),3), P=signif(prim$pval,3),
    STEIGER_OK=if (!inherits(st,"try-error")) st$correct_causal_direction else NA)
}

## --------------------------- PART A: cis-MR ----------------------------- ##
message("\n### PART A — 면역유전자 발현(cis) → 천식·우울/불안")
buildA <- function(g) {
  sub <- cis[GeneSymbol==g]
  if (!nrow(sub)) { message("  [", g, "] cis-eQTL 없음(전혈 미발현?) → 건너뜀"); return(NULL) }
  e <- format_data(as.data.frame(sub), type="exposure", snp_col="SNP", beta_col="BETA", se_col="SE",
        effect_allele_col="AssessedAllele", other_allele_col="OtherAllele", eaf_col="eaf",
        pval_col="P", samplesize_col="NrSamples", chr_col="SNPChr", pos_col="SNPPos")
  e$exposure <- g; e$id.exposure <- g
  cl <- do_clump(e)
  message("  [", g, "] cis SNP ", nrow(sub), " → clump 후 ", nrow(cl), " 도구")
  cl
}
A <- rbindlist(lapply(GENES, function(g){
  e <- buildA(g); if (is.null(e)) return(NULL)
  rbindlist(lapply(OUTC, function(ou) mr_one(e, ou, g)), fill=TRUE)
}), fill=TRUE)
if (nrow(A)) {
  A[, P_FDR := signif(p.adjust(P,"BH"),3)]
  A[, DOMAIN := fifelse(OUTCOME %in% c("COA","AOA"),"asthma","psychiatric")]
  message("\n  === PART A 요약 (cis-MR) ===")
  print(A[order(EXPOSURE,OUTCOME), .(EXPOSURE,OUTCOME,N_IV,F_MEAN,METHOD,OR,OR_LCI,OR_UCI,P,P_FDR,STEIGER_OK)])
  fwrite(A, "mr_immune_cis.csv")
  ## 공유축: 유전자가 천식(≥1) AND 정신(≥1) 양쪽에 P<0.05
  shared <- A[P<0.05, .(asthma=any(DOMAIN=="asthma"), psych=any(DOMAIN=="psychiatric")), by=EXPOSURE][asthma & psych]
  message("\n  ★ 천식·우울 양쪽에 유의(P<0.05)한 유전자: ",
          if (nrow(shared)) paste(shared$EXPOSURE, collapse=", ") else "없음")
}

## --------------------------- PART B: polygenic 보조 --------------------- ##
message("\n### PART B — 호산구·CRP polygenic → 천식·우울/불안 (보조)")
prep_exposure_poly <- function(tr) {
  dt <- as.data.table(readRDS(file.path(HARM, paste0(tr,".rds"))))
  sig <- dt[P < IV_PVAL_POLY]; if (!nrow(sig)) return(NULL)
  e <- format_data(as.data.frame(sig), type="exposure", snp_col="SNP", beta_col="BETA", se_col="SE",
        effect_allele_col="EA", other_allele_col="NEA", eaf_col="EAF", pval_col="P",
        samplesize_col="N", chr_col="CHR", pos_col="POS")
  e$exposure <- tr; e$id.exposure <- tr
  message("  [", tr, "] p<5e-8 ", nrow(sig), " → clump..."); cl <- do_clump(e)
  message("    clump 후 ", nrow(cl), " 도구"); cl
}
B <- rbindlist(lapply(c("EOS","CRP"), function(tr){
  e <- prep_exposure_poly(tr); if (is.null(e) || nrow(e)<MIN_IV) return(NULL)
  rbindlist(lapply(OUTC, function(ou) mr_one(e, ou, tr)), fill=TRUE)
}), fill=TRUE)
if (nrow(B)) {
  B[, P_FDR := signif(p.adjust(P,"BH"),3)]
  message("\n  === PART B 요약 (polygenic) ===")
  print(B[, .(EXPOSURE,OUTCOME,N_IV,F_MEAN,OR,OR_LCI,OR_UCI,P,P_FDR,STEIGER_OK)])
  fwrite(B, "mr_immune_poly.csv")
}

## --------------------------- 통합 그림 --------------------------------- ##
ALL <- rbindlist(list(A,B), fill=TRUE)
if (nrow(ALL)) {
  fwrite(ALL, "mr_immune_ALL.csv")
  ALL[, lab := paste0(EXPOSURE," → ",OUTCOME)]
  ALL[, grp := fifelse(EXPOSURE %in% GENES, "cis-eQTL (gene expression)", "polygenic (EOS/CRP)")]
  ALL[, sig := !is.na(P_FDR) & P_FDR<0.05]
  setorder(ALL, grp, EXPOSURE, OUTCOME)
  ALL[, lab := factor(lab, levels=rev(unique(lab)))]
  p <- ggplot(ALL, aes(OR, lab)) +
    geom_vline(xintercept=1, linetype=2, color="grey55") +
    geom_errorbarh(aes(xmin=OR_LCI, xmax=OR_UCI), height=0.2, color="grey35") +
    geom_point(aes(fill=sig), shape=21, size=2.6, color="grey20") +
    scale_fill_manual(values=c(`TRUE`="#B2182B",`FALSE`="white"), guide="none") +
    scale_x_log10() + facet_grid(grp~., scales="free_y", space="free_y") +
    labs(title="면역축 MR: 발현/노출 → 천식·우울/불안", x="OR (95% CI)", y=NULL) +
    theme_minimal(base_size=10)
  ggsave("fig_immune_mr.png", p, width=8, height=9, dpi=300, bg="white")
  message("\n  fig_immune_mr.png 저장")
}

message("\n================= Step 7 완료 =================")
message("산출: mr_immune_cis.csv, mr_immune_poly.csv, mr_immune_ALL.csv, fig_immune_mr.png")
message("해석: cis-MR에서 천식·우울 '양쪽'에 유의한 면역유전자 = 공유 인과축 직접 증거.")
message("주의: eQTLGen=전혈. IL13/IL33/TSLP 등 조직특이 유전자는 cis 없을 수 있음(로그 확인).")
writeLines(capture.output(sessionInfo()), file.path(MRD, paste0("_step7_sessionInfo_", Sys.Date(), ".txt")))
