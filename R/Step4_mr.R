###############################################################################
## Step 4 — 양방향 멘델 무작위화 (Bidirectional MR)
## Project: 01_GWAS_Asthma_MDD
##
## 인과방향: 천식(COA/AOA) ↔ 우울·불안(MDD/ANX)
##   정방향: 천식 → 정신질환 (유전적 소인)
##   역방향: 정신질환 → 천식 (음성/보완 대조)
## 방법: IVW(주) + Egger + weighted median + mode, Egger절편(다면발현),
##       이질성 Q, Steiger(방향), leave-one-out, MR-PRESSO, 검정력.
##
## 08 계승: 로컬 plink clumping — 한글·공백 경로를 plink.exe가 못 읽으므로
##   1000G EUR을 ASCII 경로(PATH/TO/ldref_ascii)로 복사해서 사용.
## ⚠️ LDSC에서 AOA×MDD/ANX rg 유의(0.25)였으니 MR도 신호 기대 가능(08과 다름).
###############################################################################

suppressMessages({ library(data.table); library(TwoSampleMR); library(MRPRESSO); library(ggplot2) })

## --------------------------- 경로 --------------------------------------- ##
PROJ    <- "PATH/TO/PROJECT"
HARM    <- file.path(PROJ, "harmonized")
MRDIR   <- file.path(PROJ, "mr")
dir.create(MRDIR, showWarnings = FALSE, recursive = TRUE)
setwd(MRDIR)

## --------------------------- 설정 --------------------------------------- ##
IV_PVAL     <- 5e-8; IV_PVAL_REL <- 5e-6
CLUMP_R2    <- 0.001; CLUMP_KB <- 10000
MIN_IV      <- 3; F_MIN <- 10; SEED <- 20260714
set.seed(SEED)

ASTHMA <- c("COA", "AOA")
MENTAL <- c("MDD", "ANX")          # 1차 내재화. 필요시 BIP,SCZ,LUNG 추가

## 표본크기(검정력·보고용). ⚠️ 케이스/대조는 근사 — 최종은 PI 확정.
N_INFO <- data.table(
  TRAIT  = c("COA","AOA","MDD","ANX"),
  N_CASE = c(13962, 26582, 310128, 116304),
  N_CTRL = c(300671, 300671, 1035355, 707209),
  PREV   = c(0.10, 0.05, 0.15, 0.16))

## --------------------- 로컬 clumping 준비 (ASCII 경로) ------------------- ##
USE_LOCAL   <- TRUE
LOCAL_BFILE <- "PATH/TO/ldref_ascii/EUR"     # 확장자 없이
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

run_mr <- function(ex, ou, pval = IV_PVAL, tag = "main") {
  message("\n  ── ", ex, " → ", ou, " [", tag, "] ──")
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
                  data=as.data.frame(dat), NbDistribution=1000, SignifThreshold=0.05), silent=TRUE)
    if (inherits(presso,"try-error")) presso <- NULL else
      message("    MR-PRESSO Global p=", format(presso$`MR-PRESSO results`$`Global Test`$Pvalue))
  } else message("    MR-PRESSO 생략(<10 도구변수)")
  saveRDS(list(dat=dat,res=res,plei=plei,het=het,steiger=st,loo=loo,presso=presso,Fmean=Fmean,R2=R2t),
          paste0("MR_", tag, "_", ex, "_", ou, ".rds"))
  ivw <- res[res$method=="Inverse variance weighted", ]
  data.table(TAG=tag, EXPOSURE=ex, OUTCOME=ou, N_IV=nrow(dat), F_MEAN=round(Fmean,1), R2=round(R2t,5),
    IVW_B=round(ivw$b,4), IVW_SE=round(ivw$se,4), IVW_OR=round(exp(ivw$b),3),
    IVW_LCI=round(exp(ivw$b-1.96*ivw$se),3), IVW_UCI=round(exp(ivw$b+1.96*ivw$se),3), IVW_P=signif(ivw$pval,3),
    EGGER_B=round(res$b[res$method=="MR Egger"],4), WM_B=round(res$b[res$method=="Weighted median"],4),
    MODE_B=round(res$b[res$method=="Weighted mode"],4), CONCORDANT=concordant,
    EGGER_INT_P=if (!is.null(plei)) signif(plei$pval,3) else NA_real_,
    Q_P=if (!is.null(het) && nrow(het[het$method=="Inverse variance weighted",])) signif(het$Q_pval[het$method=="Inverse variance weighted"],3) else NA_real_,
    STEIGER_OK=if (!is.null(st)) st$correct_causal_direction else NA,
    PRESSO_P=if (!is.null(presso)) presso$`MR-PRESSO results`$`Global Test`$Pvalue else NA_real_)
}

## --------------------------- PART A 정방향 (천식→정신) ------------------- ##
message("\n### PART A — 정방향: 천식 → 우울/불안")
fwd <- rbindlist(lapply(ASTHMA, function(ex) rbindlist(lapply(MENTAL, function(ou) run_mr(ex, ou)), fill=TRUE)), fill=TRUE)
if (nrow(fwd)) { fwd[, IVW_P_FDR := signif(p.adjust(IVW_P,"BH"),3)]
  message("\n  === 정방향 요약 ==="); print(fwd[, .(EXPOSURE,OUTCOME,N_IV,F_MEAN,IVW_OR,IVW_LCI,IVW_UCI,IVW_P,IVW_P_FDR,CONCORDANT)])
  fwrite(fwd, "mr_summary_forward.csv") }

## --------------------------- PART B 역방향 (정신→천식) ------------------- ##
message("\n### PART B — 역방향: 우울/불안 → 천식")
rev <- rbindlist(lapply(MENTAL, function(ex) rbindlist(lapply(ASTHMA, function(ou) run_mr(ex, ou, tag="reverse")), fill=TRUE)), fill=TRUE)
if (nrow(rev)) { rev[, IVW_P_FDR := signif(p.adjust(IVW_P,"BH"),3)]
  message("\n  === 역방향 요약 ==="); print(rev[, .(EXPOSURE,OUTCOME,N_IV,F_MEAN,IVW_OR,IVW_LCI,IVW_UCI,IVW_P,IVW_P_FDR,CONCORDANT)])
  fwrite(rev, "mr_summary_reverse.csv") }

## --------------------------- PART C 민감도(완화 5e-6) ------------------- ##
message("\n### PART C — 민감도: 완화 임계값(5e-6, 주분석 아님)")
need <- if (nrow(fwd)) fwd[N_IV < 10, .(EXPOSURE,OUTCOME)] else data.table()
if (nrow(need)) {
  relax <- rbindlist(lapply(seq_len(nrow(need)), function(i) run_mr(need$EXPOSURE[i], need$OUTCOME[i], pval=IV_PVAL_REL, tag="relaxed")), fill=TRUE)
  if (nrow(relax)) { fwrite(relax, "mr_summary_relaxed.csv"); print(relax[, .(EXPOSURE,OUTCOME,N_IV,F_MEAN,IVW_OR,IVW_P)]) }
} else message("  모든 조합 도구변수 10+ → 생략")

## --------------------------- PART D 검정력 ------------------------------ ##
message("\n### PART D — 검정력 분석 (⚠️공식·가정 PI 확정)")
min_or_80 <- function(R2,N,K,alpha=0.05,power=0.80){ if (is.na(R2)||R2<=0) return(NA_real_)
  crit<-qchisq(1-alpha,1); ncp<-uniroot(function(x) pchisq(crit,1,x,lower.tail=FALSE)-power, c(0,200))$root
  exp(sqrt(ncp*(1-R2)/(R2*N*K*(1-K)))) }
if (nrow(fwd)) {
  pw <- merge(fwd[, .(EXPOSURE,OUTCOME,N_IV,R2,IVW_OR,IVW_LCI,IVW_UCI,IVW_P)], N_INFO, by.x="OUTCOME", by.y="TRAIT", all.x=TRUE)
  pw[, N_TOTAL:=N_CASE+N_CTRL][, K:=N_CASE/N_TOTAL][, MIN_OR_80:=mapply(min_or_80,R2,N_TOTAL,K)]
  print(pw[, .(EXPOSURE,OUTCOME,N_IV,R2=round(R2,5),MIN_OR_80=round(MIN_OR_80,3),OBS_OR=IVW_OR,P=IVW_P)])
  fwrite(pw, "mr_power_analysis.csv")
}

## --------------------------- PART E 그림 ------------------------------- ##
message("\n### PART E — forest plot")
if (nrow(fwd)) {
  fp <- copy(fwd); fp[, label:=paste0(EXPOSURE," → ",OUTCOME)][, label:=factor(label, levels=rev(label))]
  p <- ggplot(fp, aes(IVW_OR, label)) + geom_vline(xintercept=1, linetype=2, colour="grey50") +
    geom_errorbarh(aes(xmin=IVW_LCI, xmax=IVW_UCI), height=0.2, colour="#2f5f9e") +
    geom_point(size=2.5, colour="#2f5f9e") + scale_x_log10() +
    labs(x="OR (95% CI), IVW", y=NULL, title="천식 → 우울/불안 인과효과 (IVW)") + theme_minimal(base_size=11)
  ggsave("fig_forest_forward.png", p, width=8, height=4.5, dpi=150)
  for (i in seq_len(nrow(fwd))) { f<-paste0("MR_main_",fwd$EXPOSURE[i],"_",fwd$OUTCOME[i],".rds"); if (!file.exists(f)) next
    o<-readRDS(f); tryCatch({ png(paste0("fig_scatter_",fwd$EXPOSURE[i],"_",fwd$OUTCOME[i],".png"),800,700,res=120)
      print(mr_scatter_plot(o$res,o$dat)[[1]]); dev.off() }, error=function(e) NULL) }
  message("  fig_forest_forward.png + scatter 저장")
}

## --------------------------- 종합 -------------------------------------- ##
message("\n### 종합")
if (nrow(fwd)) {
  message(sprintf("  정방향 FDR유의 %d개 | 약한도구(F<10) %d개 | 방향불일치 %d개",
    sum(fwd$IVW_P_FDR<0.05,na.rm=TRUE), sum(fwd$F_MEAN<F_MIN,na.rm=TRUE), sum(!fwd$CONCORDANT,na.rm=TRUE)))
  message("  점검: 정/역 방향 어느 쪽이 유의한가(인과방향) · Steiger · Egger절편(다면발현) · 역방향 음성대조.")
}
message("\n저장: ", MRDIR)
writeLines(capture.output(sessionInfo()), file.path(MRDIR, paste0("_step4_sessionInfo_", Sys.Date(), ".txt")))
