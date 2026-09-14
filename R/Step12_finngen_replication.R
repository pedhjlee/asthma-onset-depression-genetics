## Step12_finngen_replication.R
## 우울 -> 천식 MR의 FinnGen 재현 (2026-09-11)
## 입력: raw/finngen_R13_J10_ASTHMA_EXMORE.gz, raw/finngen_R13_ASTHMA_CHILD_EXMORE.gz,
##       raw/PGC_UKB_depression_genome-wide.txt, raw/daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.gz
## 출력: harmonized/{FG_ASTHMA,FG_COA,MDD_HOWARD,MDD_WRAY}.rds,
##       mr/REP_*.rds, Result/tables_csv/TableS11_FinnGen_replication.csv (잠정)
##
## ★ 2026-09-12 갱신: MR-PRESSO 를 seed 20260911·10,000회로 통일 (원고 최종값과 동일 설정).
##   이 스크립트의 Steiger(directionality_test)와 Table S11 은 잠정값임.
##   최종 Steiger(이분형·유효표본수 통일)와 최종 S8/S11/S12 는 반드시
##   Step13_final_sensitivity.R 을 이어서 실행해 만들 것. 이 파일만 다시 돌리면 S11 이 잠정값으로 덮임.

library(data.table); library(TwoSampleMR); library(MRPRESSO)
PRESSO_SEED <- 20260911; PRESSO_N <- 10000
PROJ  <- "PATH/TO/PROJECT"
RAW   <- file.path(PROJ,"raw"); HARM <- file.path(PROJ,"harmonized")
OUT   <- file.path(PROJ,"mr");  TABO <- file.path(PROJ,"Result","tables_csv")
dir.create(HARM, showWarnings=FALSE)

## 참조패널: plink 는 공백/한글 경로를 못 읽으므로 C 드라이브로 복사해 사용
BFILE <- "PATH/TO/ldref_ascii/EUR"
PLINK <- genetics.binaRies::get_plink_binary()
if (!file.exists(paste0(BFILE,".bed"))) {
  dir.create("PATH/TO/ldref_ascii", showWarnings=FALSE)
  for (e in c("bed","bim","fam"))
    file.copy(file.path(PROJ,"ld_ref",paste0("EUR.",e)), paste0(BFILE,".",e))
}

ACGT <- c("A","C","G","T")
clean <- function(d, tag) {
  n0 <- nrow(d)
  d <- d[!is.na(SNP) & grepl("^rs", SNP)]
  d <- d[A1 %in% ACGT & A2 %in% ACGT & A1 != A2]
  d <- d[!is.na(BETA) & !is.na(SE) & SE > 0 & is.finite(BETA)]
  d <- d[is.na(EAF) | (EAF > 0.01 & EAF < 0.99)]
  d <- d[!SNP %in% d[duplicated(SNP), SNP]]
  cat(sprintf("  %s: %s -> %s\n", tag, format(n0,big.mark=","), format(nrow(d),big.mark=",")))
  d[]
}

## ---------------- 1) harmonize ----------------
## FinnGen: 효과는 alt 대립유전자 기준, GRCh38 이지만 rsID 로 매칭
fg <- list(FG_ASTHMA = c("finngen_R13_J10_ASTHMA_EXMORE.gz",   61196, 250433),
           FG_COA    = c("finngen_R13_ASTHMA_CHILD_EXMORE.gz",  8428, 250433))
for (nm in names(fg)) {
  p <- file.path(RAW, fg[[nm]][1]); ca <- as.numeric(fg[[nm]][2]); co <- as.numeric(fg[[nm]][3])
  cat("\n[",nm,"]\n")
  d <- fread(p, select=c("#chrom","pos","ref","alt","rsids","beta","sebeta","pval","af_alt"))
  setnames(d, c("CHR","BP","A2","A1","SNP","BETA","SE","P","EAF"))
  d[, SNP := tstrsplit(SNP, ",", fixed=TRUE)[[1]]]
  d[, `:=`(A1=toupper(A1), A2=toupper(A2), N=4/(1/ca+1/co))]
  saveRDS(clean(d[, .(SNP,CHR,BP,A1,A2,EAF,BETA,SE,P,N)], nm), file.path(HARM, paste0(nm,".rds")))
}

## Howard 2019: A1 = 효과 대립유전자 (소문자), N 열 없음 -> 유효표본수 직접 계산
cat("\n[ MDD_HOWARD ]\n")
d <- fread(file.path(RAW,"PGC_UKB_depression_genome-wide.txt"))
setnames(d, c("MarkerName","A1","A2","Freq","LogOR","StdErrLogOR","P"),
            c("SNP","A1","A2","EAF","BETA","SE","P"))
d[, `:=`(A1=toupper(A1), A2=toupper(A2), CHR=NA_integer_, BP=NA_integer_,
         N=4/(1/170756 + 1/329443))]
saveRDS(clean(d[, .(SNP,CHR,BP,A1,A2,EAF,BETA,SE,P,N)], "MDD_HOWARD"),
        file.path(HARM,"MDD_HOWARD.rds"))

## Wray 2018 (daner): 가교 분석용. 도구변수 부족으로 최종 분석에는 쓰지 못함
pw <- file.path(RAW,"daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.gz")
if (file.exists(pw)) {
  cat("\n[ MDD_WRAY ]\n")
  d <- fread(pw)
  fa <- grep("^FRQ_A_", names(d), value=TRUE)[1]; fu <- grep("^FRQ_U_", names(d), value=TRUE)[1]
  ca <- as.numeric(sub("FRQ_A_","",fa)); co <- as.numeric(sub("FRQ_U_","",fu))
  d[, `:=`(BETA=log(OR), EAF=(get(fa)*ca+get(fu)*co)/(ca+co), N=4/(1/ca+1/co))]
  saveRDS(clean(d[, .(SNP,CHR,BP,A1=toupper(A1),A2=toupper(A2),EAF,BETA,SE,P,N)], "MDD_WRAY"),
          file.path(HARM,"MDD_WRAY.rds"))
}

## ---------------- 2) 재현 MR ----------------
rd <- function(x) readRDS(file.path(HARM, paste0(x,".rds")))
run_mr <- function(expo, outc, tag, kb=10000) {
  cat("\n=====", tag, "=====\n")
  e <- rd(expo)[P < 5e-8]; cat("  P<5e-8:", nrow(e), "\n")
  if (nrow(e) < 3) { cat("  중단\n"); return(invisible(NULL)) }
  cl <- try(ieugwasr::ld_clump(data.frame(rsid=e$SNP, pval=e$P),
        clump_kb=kb, clump_r2=0.001, bfile=BFILE, plink_bin=PLINK), silent=TRUE)
  if (inherits(cl,"try-error") || nrow(cl) < 3) { cat("  도구변수 부족 -> 중단\n"); return(invisible(NULL)) }
  cat("  clump 후:", nrow(cl), "\n")
  ex <- format_data(as.data.frame(e[SNP %in% cl$rsid]), type="exposure", snp_col="SNP",
        beta_col="BETA", se_col="SE", eaf_col="EAF", pval_col="P",
        effect_allele_col="A1", other_allele_col="A2", samplesize_col="N")
  ou <- format_data(as.data.frame(rd(outc)[SNP %in% cl$rsid]), type="outcome", snp_col="SNP",
        beta_col="BETA", se_col="SE", eaf_col="EAF", pval_col="P",
        effect_allele_col="A1", other_allele_col="A2", samplesize_col="N")
  ## action=3: 회문형 SNP 제거 (핀란드 인구는 빈도 기반 방향 추정이 위험)
  dat <- harmonise_data(ex, ou, action=3); dat <- dat[dat$mr_keep, ]
  cat("  harmonise 후:", nrow(dat), "\n")
  if (nrow(dat) < 3) { cat("  중단\n"); return(invisible(NULL)) }
  res  <- mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
  plei <- mr_pleiotropy_test(dat); het <- mr_heterogeneity(dat)
  loo  <- mr_leaveoneout(dat); stg <- directionality_test(dat)
  set.seed(PRESSO_SEED)
  pr <- if (nrow(dat) >= 10) try(mr_presso(BetaOutcome="beta.outcome", BetaExposure="beta.exposure",
        SdOutcome="se.outcome", SdExposure="se.exposure", OUTLIERtest=TRUE, DISTORTIONtest=TRUE,
        data=dat, NbDistribution=PRESSO_N, SignifThreshold=0.05), silent=TRUE) else NULL
  print(as.data.table(res)[, .(method, nsnp, OR=round(exp(b),3),
        LCI=round(exp(b-1.96*se),3), UCI=round(exp(b+1.96*se),3), P=signif(pval,3))])
  saveRDS(list(dat=dat,res=res,plei=plei,het=het,loo=loo,steiger=stg,presso=pr,
               presso_seed=PRESSO_SEED, presso_n=PRESSO_N),
          file.path(OUT, paste0("REP_",tag,".rds")))
  invisible(NULL)
}
run_mr("MDD_HOWARD","FG_ASTHMA","HOWARD_FGASTHMA")
run_mr("MDD_HOWARD","FG_COA",   "HOWARD_FGCOA")
## 가교 시도: 10 Mb / 1 Mb 모두 clump 후 2개 -> 분석 불가 (원고 Methods 에 기술)
run_mr("MDD_WRAY","FG_ASTHMA","WRAY_FGASTHMA")
run_mr("MDD_WRAY","FG_COA",   "WRAY_FGCOA")

## ---------------- 3) 전체 vs 소아 차이 검정 ----------------
## 대조군 공유 + 소아 사례가 전체에 포함 -> 표본오차 상관 0.5 가정
a <- as.data.table(readRDS(file.path(OUT,"REP_HOWARD_FGASTHMA.rds"))$dat)
b_ <- as.data.table(readRDS(file.path(OUT,"REP_HOWARD_FGCOA.rds"))$dat)
m <- merge(a[, .(SNP, bx=beta.exposure, b1=beta.outcome, s1=se.outcome)],
           b_[, .(SNP, b2=beta.outcome, s2=se.outcome)], by="SNP")
r <- 0.5
m[, `:=`(d=b1-b2, sd=sqrt(s1^2 + s2^2 - 2*r*s1*s2))]
fit <- lm(d ~ bx - 1, weights=1/sd^2, data=m)
bb <- coef(summary(fit))[1,1]; se <- coef(summary(fit))[1,2]/min(1, summary(fit)$sigma)
cat(sprintf("\n전체 vs 소아 (nSNP %d): OR 비 %.3f (%.3f-%.3f), P_difference = %.3g\n",
            nrow(m), exp(bb), exp(bb-1.96*se), exp(bb+1.96*se), 2*pnorm(-abs(bb/se))))
DIFF <- data.table(Outcome="Difference (all vs childhood)", `N IV`=nrow(m),
  `IVW OR (95% CI)`=sprintf("%.2f (%.2f-%.2f)", exp(bb), exp(bb-1.96*se), exp(bb+1.96*se)),
  `IVW P`=signif(2*pnorm(-abs(bb/se)),3))

## ---------------- 4) Table S11 (잠정: Steiger 열은 Step13 에서 최종화) ----------------
NM <- c(HOWARD_FGASTHMA="FinnGen asthma (all)",
        HOWARD_FGCOA   ="FinnGen childhood-onset asthma (age<16)")
ORc <- function(b, se) sprintf("%.2f (%.2f-%.2f)", exp(b), exp(b-1.96*se), exp(b+1.96*se))
tabS11 <- rbindlist(lapply(names(NM), function(tg) {
  o <- readRDS(file.path(OUT, paste0("REP_",tg,".rds")))
  r <- as.data.table(o$res); g <- function(mm,w) r[method==mm][[w]]
  pr <- o$presso
  corr <- if (!is.null(pr) && !inherits(pr,"try-error")) pr$`Main MR results`[2,] else NULL
  data.table(Outcome=NM[tg], `N IV`=r$nsnp[1],
    `F (mean)`=round(mean(o$dat$beta.exposure^2/o$dat$se.exposure^2),1),
    `IVW OR (95% CI)`=ORc(g("Inverse variance weighted","b"), g("Inverse variance weighted","se")),
    `IVW P`=signif(g("Inverse variance weighted","pval"),3),
    `Weighted median`=ORc(g("Weighted median","b"), g("Weighted median","se")),
    `Weighted mode`  =ORc(g("Weighted mode","b"),   g("Weighted mode","se")),
    `MR-Egger`       =ORc(g("MR Egger","b"),        g("MR Egger","se")),
    `Egger intercept P`=signif(o$plei$pval,3),
    `Cochran Q P`=signif(o$het$Q_pval[o$het$method=="Inverse variance weighted"],3),
    `MR-PRESSO global P`=if (!is.null(pr) && !inherits(pr,"try-error"))
        as.character(pr$`MR-PRESSO results`$`Global Test`$Pvalue) else NA_character_,
    `Outlier-corrected OR`=if (!is.null(corr) && !is.na(corr$`Causal Estimate`))
        ORc(corr$`Causal Estimate`, corr$Sd) else "no outlier",
    `Steiger P`=signif(o$steiger$steiger_pval,3))
}))
tabS11 <- rbind(tabS11, DIFF, fill=TRUE)
fwrite(tabS11, file.path(TABO,"TableS11_FinnGen_replication.csv"))
print(tabS11[, .(Outcome, `N IV`, `IVW OR (95% CI)`, `IVW P`)])
cat("\nDONE — 이어서 Step13_final_sensitivity.R 을 실행해 S11/S12 를 최종화할 것\n")
