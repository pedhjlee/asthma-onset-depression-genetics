###############################################################################
## Step 5 — 폐기능 축 MR: 방향성(direction) + 판별(discriminant)
## Project: 01_GWAS_Asthma_MDD
##
## 목적
##   "호흡곤란 증상축"에서 폐기능의 역할을 인과적으로 규명. Step4(천식↔정신)와
##   별개로, 폐기능(FEV1·FVC 메인, FEV1/FVC ratio=LUNG 폐쇄대조)에 대해:
##
##   PART A  천식 → 폐기능      : 천식 유전소인이 폐기능을 낮추나(질병→결과, 예상방향)
##   PART B  폐기능 → 천식      : 기저 낮은 폐기능이 천식 위험인가(체질적 기도협착 방향)
##   PART C  우울/불안 → 폐기능 : (판별) 우울이 폐기능에 인과효과? 볼륨(FVC)이냐 폐쇄냐
##   PART D  폐기능 → 우울/불안 : 낮은 폐기능이 우울을 유발하나(BMC Med 280k 문헌방향)
##
## 설계 근거(역인과 방어)
##   관찰연구상 "천식→폐기능저하"는 질병상태 관점에서 성립(암→체중감소 유형).
##   유전 관점에선 "기저 기도구조→천식"도 가능 → A/B 양방향으로 어느 쪽이 우세한지 확정.
##   우울-폐기능은 볼륨(FVC↓, 제한성/PRISm)과 더 연관(문헌) → C에서 FVC vs FEV1 대비,
##   LUNG(ratio)까지 넣어 "폐쇄가 아니라 볼륨"인지 직접 판별.
##
## ⚠️ 보고 스케일
##   폐기능=연속형 → 효과는 beta(결과 SD 변화)로 보고(OR 아님).
##   천식/정신질환=binary(결과일 때) → OR=exp(beta)로 보고.
##   run_mr()가 결과형질 TYPE 보고 자동 분기.
##
## 인프라: Step4와 동일(로컬 plink clumping, IVW/Egger/WM/mode, F, Steiger, LOO, PRESSO).
##   이 스크립트는 self-contained — Step4 없이 단독 실행 가능.
## 실행: source("step5_mr_lung.R")  (harmonized/FEV1.rds·FVC.rds·LUNG.rds 필요)
###############################################################################

suppressMessages({ library(data.table); library(TwoSampleMR); library(MRPRESSO); library(ggplot2) })

## --------------------------- 경로 --------------------------------------- ##
PROJ  <- "PATH/TO/PROJECT"
HARM  <- file.path(PROJ, "harmonized")
MRDIR <- file.path(PROJ, "mr")
dir.create(MRDIR, showWarnings = FALSE, recursive = TRUE)
setwd(MRDIR)

## --------------------------- 설정 --------------------------------------- ##
IV_PVAL   <- 5e-8; IV_PVAL_REL <- 5e-6
CLUMP_R2  <- 0.001; CLUMP_KB <- 10000
MIN_IV    <- 3; F_MIN <- 10; SEED <- 20260714
set.seed(SEED)

ASTHMA <- c("COA", "AOA")
MENTAL <- c("MDD", "ANX")
LUNG   <- c("FEV1", "FVC", "LUNG")          # FEV1·FVC 메인, LUNG=FEV1/FVC ratio(폐쇄대조)

## 결과형질 스케일 분기용
TRAIT_TYPE <- c(COA="binary", AOA="binary", MDD="binary", ANX="binary",
                FEV1="continuous", FVC="continuous", LUNG="continuous")

## --------------------- 로컬 clumping 준비 (ASCII 경로) ------------------- ##
USE_LOCAL   <- TRUE
LOCAL_BFILE <- "PATH/TO/ldref_ascii/EUR"
LOCAL_PLINK <- NULL
{  ## 한글·공백 경로 → ASCII(C:)로 EUR.{bed,bim,fam} 복사(없을 때만)
  ascii_dir <- "PATH/TO/ldref_ascii"; dir.create(ascii_dir, showWarnings = FALSE, recursive = TRUE)
  src <- file.path(PROJ, "ld_ref")
  for (ext in c("bed","bim","fam")) {
    d <- file.path(ascii_dir, paste0("EUR.", ext)); s <- file.path(src, paste0("EUR.", ext))
    if (!file.exists(d) && file.exists(s)) { message("  LD패널 ASCII 복사: EUR.", ext); file.copy(s, d, overwrite = TRUE) }
  }
}
if (USE_LOCAL) {
  if (is.null(LOCAL_PLINK)) {
    if (!requireNamespace("genetics.binaRies", quietly = TRUE)) remotes::install_github("MRCIEU/genetics.binaRies")
    LOCAL_PLINK <- genetics.binaRies::get_plink_binary()
  }
  ok <- file.exists(paste0(LOCAL_BFILE, ".bed"))
  message("  로컬 clumping: ", if (ok) "준비됨" else "★참조패널 없음 → API 폴백")
  if (!ok) USE_LOCAL <- FALSE
}

## --------------------------- 공통 함수 ---------------------------------- ##
do_clump <- function(dat) {
  if (USE_LOCAL) {
    out <- try(ieugwasr::ld_clump(
      dplyr::tibble(rsid = dat$SNP, pval = dat$pval.exposure, id = dat$id.exposure),
      clump_kb = CLUMP_KB, clump_r2 = CLUMP_R2, bfile = LOCAL_BFILE, plink_bin = LOCAL_PLINK), silent = TRUE)
    if (!inherits(out, "try-error")) return(dat[dat$SNP %in% out$rsid, ])
    message("    로컬 clumping 실패 → API 시도")
  }
  out <- try(clump_data(dat, clump_r2 = CLUMP_R2, clump_kb = CLUMP_KB, pop = "EUR"), silent = TRUE)
  if (inherits(out, "try-error")) { message("    ★clumping 실패"); return(NULL) }
  out
}
prep_exposure <- function(tr, pval = IV_PVAL) {
  dt <- as.data.table(readRDS(file.path(HARM, paste0(tr, ".rds"))))
  sig <- dt[P < pval]; if (!nrow(sig)) { message("    ", tr, ": p<", pval, " SNP 없음"); return(NULL) }
  message("    ", tr, ": p<", pval, " SNP = ", nrow(sig))
  e <- format_data(as.data.frame(sig), type = "exposure", snp_col="SNP", beta_col="BETA", se_col="SE",
                   effect_allele_col="EA", other_allele_col="NEA", eaf_col="EAF", pval_col="P",
                   samplesize_col="N", chr_col="CHR", pos_col="POS")
  e$exposure <- tr; e$id.exposure <- tr
  cl <- do_clump(e); if (is.null(cl) || !nrow(cl)) return(NULL)
  message("    clumping 후: ", nrow(cl), " SNP"); cl
}
prep_outcome <- function(tr, snps) {
  dt <- as.data.table(readRDS(file.path(HARM, paste0(tr, ".rds"))))
  sub <- dt[SNP %in% snps]; if (!nrow(sub)) return(NULL)
  o <- format_data(as.data.frame(sub), type = "outcome", snp_col="SNP", beta_col="BETA", se_col="SE",
                   effect_allele_col="EA", other_allele_col="NEA", eaf_col="EAF", pval_col="P", samplesize_col="N")
  o$outcome <- tr; o$id.outcome <- tr; o
}

## run_mr — 결과형질 스케일 자동 분기(연속=beta, binary=OR)
run_mr <- function(ex, ou, pval = IV_PVAL, tag = "main") {
  out_bin <- identical(TRAIT_TYPE[[ou]], "binary")
  message("\n  ── ", ex, " → ", ou, " [", tag, "] (결과 ", TRAIT_TYPE[[ou]], ") ──")
  e <- prep_exposure(ex, pval); if (is.null(e) || nrow(e) < MIN_IV) { message("    도구변수 부족 → 건너뜀"); return(NULL) }
  o <- prep_outcome(ou, e$SNP); if (is.null(o)) { message("    결과데이터 없음"); return(NULL) }
  dat <- harmonise_data(e, o, action = 2); dat <- dat[dat$mr_keep, ]
  message("    최종 도구변수: ", nrow(dat)); if (nrow(dat) < MIN_IV) { message("    <", MIN_IV, " → 건너뜀"); return(NULL) }
  dat$F_stat <- (dat$beta.exposure/dat$se.exposure)^2; Fmean <- mean(dat$F_stat)
  dat$R2 <- 2*dat$eaf.exposure*(1-dat$eaf.exposure)*dat$beta.exposure^2; R2t <- sum(dat$R2, na.rm=TRUE)
  message(sprintf("    평균 F=%.1f%s | R2=%.4f", Fmean, if (Fmean<F_MIN) " ★약한도구" else "", R2t))
  res <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
  print(res[, c("method","nsnp","b","se","pval")])
  concordant <- length(unique(sign(res$b))) == 1
  message("    방법 간 방향 일치: ", if (concordant) "예" else "★아니오")
  plei <- try(mr_pleiotropy_test(dat), silent=TRUE); if (inherits(plei,"try-error")) plei <- NULL else
    message(sprintf("    Egger절편 p=%.3g%s", plei$pval, if (plei$pval<0.05) " ★다면발현의심" else ""))
  het <- try(mr_heterogeneity(dat), silent=TRUE); if (inherits(het,"try-error")) het <- NULL
  st <- try(directionality_test(dat), silent=TRUE); if (inherits(st,"try-error")) st <- NULL else
    message("    Steiger 방향올바름: ", st$correct_causal_direction, sprintf(" (p=%.3g)", st$steiger_pval))
  loo <- try(mr_leaveoneout(dat), silent=TRUE); if (inherits(loo,"try-error")) loo <- NULL
  presso <- NULL
  if (nrow(dat) >= 10) {
    presso <- try(mr_presso(BetaOutcome="beta.outcome", BetaExposure="beta.exposure",
                  SdOutcome="se.outcome", SdExposure="se.exposure", OUTLIERtest=TRUE, DISTORTIONtest=TRUE,
                  data=as.data.frame(dat), NbDistribution=2000, SignifThreshold=0.05), silent=TRUE)
    if (inherits(presso,"try-error")) presso <- NULL else
      message("    MR-PRESSO Global p=", format(presso$`MR-PRESSO results`$`Global Test`$Pvalue))
  } else message("    MR-PRESSO 생략(<10 도구변수)")
  saveRDS(list(dat=dat,res=res,plei=plei,het=het,steiger=st,loo=loo,presso=presso,Fmean=Fmean,R2=R2t),
          paste0("MR_", tag, "_", ex, "_", ou, ".rds"))
  ivw <- res[res$method=="Inverse variance weighted", ]
  b <- ivw$b; se <- ivw$se
  data.table(TAG=tag, EXPOSURE=ex, OUTCOME=ou, OUT_TYPE=TRAIT_TYPE[[ou]],
    N_IV=nrow(dat), F_MEAN=round(Fmean,1), R2=round(R2t,5),
    ## beta 스케일(항상): 연속형 결과의 주지표 = 결과 SD 변화
    IVW_B=round(b,4), IVW_SE=round(se,4),
    IVW_B_LCI=round(b-1.96*se,4), IVW_B_UCI=round(b+1.96*se,4), IVW_P=signif(ivw$pval,3),
    ## OR 스케일: 결과가 binary일 때만 의미(연속형이면 NA)
    IVW_OR=if (out_bin) round(exp(b),3) else NA_real_,
    OR_LCI=if (out_bin) round(exp(b-1.96*se),3) else NA_real_,
    OR_UCI=if (out_bin) round(exp(b+1.96*se),3) else NA_real_,
    EGGER_B=round(res$b[res$method=="MR Egger"],4), WM_B=round(res$b[res$method=="Weighted median"],4),
    MODE_B=round(res$b[res$method=="Weighted mode"],4), CONCORDANT=concordant,
    EGGER_INT_P=if (!is.null(plei)) signif(plei$pval,3) else NA_real_,
    Q_P=if (!is.null(het) && nrow(het[het$method=="Inverse variance weighted",])) signif(het$Q_pval[het$method=="Inverse variance weighted"],3) else NA_real_,
    STEIGER_OK=if (!is.null(st)) st$correct_causal_direction else NA,
    PRESSO_P=if (!is.null(presso)) presso$`MR-PRESSO results`$`Global Test`$Pvalue else NA_real_)
}

## 한 파트 실행 + 요약 출력 헬퍼
run_family <- function(ex_set, ou_set, tag, title) {
  message("\n### ", title)
  s <- rbindlist(lapply(ex_set, function(ex)
         rbindlist(lapply(ou_set, function(ou) run_mr(ex, ou, tag=tag)), fill=TRUE)), fill=TRUE)
  if (nrow(s)) {
    s[, IVW_P_FDR := signif(p.adjust(IVW_P,"BH"),3)]
    message("\n  === ", title, " 요약 ===")
    print(s[, .(EXPOSURE,OUTCOME,OUT_TYPE,N_IV,F_MEAN,IVW_B,IVW_B_LCI,IVW_B_UCI,IVW_OR,IVW_P,IVW_P_FDR,CONCORDANT,STEIGER_OK)])
    fwrite(s, paste0("mr_summary_", tag, ".csv"))
  } else message("  (유효 결과 없음)")
  s
}

## --------------------------- 실행 --------------------------------------- ##
A <- run_family(ASTHMA, LUNG, "lung_asthma2lung", "PART A — 천식 → 폐기능 (질병→결과 방향)")
B <- run_family(LUNG, ASTHMA, "lung_lung2asthma", "PART B — 폐기능 → 천식 (체질적 기도 방향)")
C <- run_family(MENTAL, LUNG, "lung_mental2lung", "PART C — 우울/불안 → 폐기능 (판별: 볼륨 vs 폐쇄)")
D <- run_family(LUNG, MENTAL, "lung_lung2mental", "PART D — 폐기능 → 우울/불안 (낮은 폐기능→우울?)")

## --------------------------- 민감도(완화 5e-6): 폐기능 노출 IV 적을 때 ---- ##
## SpiroMeta 단독(N~75k)이라 폐기능 genome-wide IV가 적을 수 있음 → 완화 임계값 보조.
message("\n### 민감도 — 폐기능 노출(B·D)에서 IV<10 이면 5e-6 재실행")
lowIV <- rbindlist(list(
  if (!is.null(B) && nrow(B)) B[N_IV < 10, .(EX=EXPOSURE, OU=OUTCOME, tg="lung_lung2asthma_relax")] else NULL,
  if (!is.null(D) && nrow(D)) D[N_IV < 10, .(EX=EXPOSURE, OU=OUTCOME, tg="lung_lung2mental_relax")] else NULL), fill=TRUE)
if (nrow(lowIV)) {
  relax <- rbindlist(lapply(seq_len(nrow(lowIV)), function(i)
    run_mr(lowIV$EX[i], lowIV$OU[i], pval=IV_PVAL_REL, tag=lowIV$tg[i])), fill=TRUE)
  if (nrow(relax)) { fwrite(relax, "mr_summary_lung_relaxed.csv")
    print(relax[, .(EXPOSURE,OUTCOME,N_IV,F_MEAN,IVW_B,IVW_OR,IVW_P)]) }
} else message("  해당 없음(모두 IV 10+ 또는 결과 없음)")

## --------------------------- 통합 요약 & 그림 --------------------------- ##
ALL <- rbindlist(list(A,B,C,D), fill=TRUE)
if (nrow(ALL)) {
  fwrite(ALL, "mr_summary_lung_ALL.csv")

  ## PART C(우울/불안→폐기능) 판별 그림: 형질별 beta(SD) forest
  fpC <- ALL[TAG=="lung_mental2lung"]
  if (nrow(fpC)) {
    fpC[, label := paste0(EXPOSURE," → ",OUTCOME)]
    fpC[, label := factor(label, levels=rev(label))]
    p <- ggplot(fpC, aes(IVW_B, label)) + geom_vline(xintercept=0, linetype=2, colour="grey50") +
      geom_errorbarh(aes(xmin=IVW_B_LCI, xmax=IVW_B_UCI), height=0.2, colour="#b0413e") +
      geom_point(size=2.5, colour="#b0413e") +
      labs(x="폐기능 변화 (SD, IVW) per 우울/불안 유전소인", y=NULL,
           title="우울/불안 → 폐기능 (판별: FVC=볼륨 vs FEV1 vs LUNG=폐쇄)") +
      theme_minimal(base_size=11)
    ggsave("fig_forest_mental2lung.png", p, width=8, height=4.5, dpi=150)
    message("  fig_forest_mental2lung.png 저장")
  }
}

## --------------------------- 해석 가이드 ------------------------------- ##
message("\n### 해석 포인트")
message("  A vs B: 천식-폐기능에서 어느 방향이 유의한가.")
message("          A(천식→폐기능↓) 우세 = 질병이 폐기능 저하시킴 / B(폐기능→천식) 우세 = 기저 기도협착이 위험.")
message("  C     : 우울/불안이 폐기능에 인과효과 있나. FVC(볼륨) beta가 LUNG(ratio,폐쇄)보다 크면 '제한성' 축.")
message("          → 우울-천식 공유축이 '단순 기도폐쇄'가 아니라 '전신/볼륨 감소'임을 시사(면역기전과 정합).")
message("  D     : 낮은 폐기능→우울(BMC Med 280k 방향) 재현되나. C·D 양방향 패턴으로 인과구조 정리.")
message("  ⚠️ 폐기능=연속형 → beta(SD)가 주지표. OR 칼럼은 결과가 binary(천식/정신)일 때만 봄.")
message("\n저장: ", MRDIR, "  (mr_summary_lung_*.csv, mr_summary_lung_ALL.csv, fig_forest_mental2lung.png)")
writeLines(capture.output(sessionInfo()), file.path(MRDIR, paste0("_step5_sessionInfo_", Sys.Date(), ".txt")))
message("== Step 5 폐기능 MR 완료 ==")
