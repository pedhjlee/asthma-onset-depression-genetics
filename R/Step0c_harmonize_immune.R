###############################################################################
## Step 0c — 면역형질 harmonize (GWAS-VCF → 공통 스키마)
## Project: 01_GWAS_Asthma_MDD  (면역축 확장)
##
## 대상(Tier 1): 호산구(EOS, Astle 2016 ebi-a-GCST004606),
##               CRP(Ligthart 2018 ieu-b-35).  둘 다 IEU OpenGWAS GWAS-VCF.
## GWAS-VCF 포맷: FORMAT=ES:SE:LP:AF:SS:ID (순서는 파일의 FORMAT 문자열로 동적 파싱)
##   ES = ALT대립유전자 기준 효과 → EA=ALT, NEA=REF, BETA=ES
##   LP = -log10(P) → P=10^(-LP);  AF=ALT빈도=EAF;  SS=N;  좌표 GRCh37
##
## 기존 Step0 finalize()/QC 재사용(STEP0_NORUN). 연속형.
## 산출: harmonized/EOS.rds, harmonized/CRP.rds
## 실행: source("Step0c_harmonize_immune.R")
###############################################################################

## ---- 기존 Step0 함수 로드 (실행부 건너뜀) --------------------------------
STEP0_NORUN <- TRUE
SRC <- "PATH/TO/PROJECT/R script/Step0_harmonize.R"
if (!file.exists(SRC)) SRC <- "Step0_harmonize.R"
source(SRC)
stopifnot(exists("finalize"), exists("RAW"), exists("HARM"))
suppressPackageStartupMessages(library(data.table))

## ---- GWAS-VCF 리더 -------------------------------------------------------
## FORMAT 문자열을 첫 데이터행에서 읽어 필드 위치를 동적으로 잡음(소스별 순서차 대응).
read_gwasvcf <- function(fname, trait, n_const = NA_real_) {
  f <- file.path(RAW, fname)
  stopifnot(file.exists(f))
  message("  읽는 중: ", fname)
  d <- fread(f, skip = "#CHROM", header = TRUE, sep = "\t", showProgress = FALSE)
  setnames(d, 1, "CHROM")                          # "#CHROM" → CHROM
  stopifnot(all(c("POS","ID","REF","ALT","FORMAT") %in% names(d)))
  smpcol <- names(d)[ncol(d)]                       # 마지막 컬럼 = 형질값(FORMAT 적용대상)
  fmt <- strsplit(as.character(d[["FORMAT"]][1]), ":", fixed = TRUE)[[1]]
  message("    FORMAT: ", paste(fmt, collapse=":"), " | 값컬럼: ", smpcol)
  parts <- tstrsplit(d[[smpcol]], ":", fixed = TRUE)
  names(parts) <- fmt
  getf <- function(k) if (k %in% fmt) suppressWarnings(as.numeric(parts[[k]])) else rep(NA_real_, nrow(d))
  ES <- getf("ES"); SE <- getf("SE"); LP <- getf("LP"); AF <- getf("AF"); SS <- getf("SS")
  ## SS(표본크기) 필드가 없는 소스(예: Astle 호산구)는 상수 N 사용
  if (all(is.na(SS)) && is.finite(n_const)) {
    SS <- rep(n_const, nrow(d)); message("    SS 필드 없음 → N=", format(n_const, big.mark=","), " 상수 적용")
  }
  P  <- 10^(-LP); P[!is.finite(P)] <- NA_real_; P[P == 0] <- 1e-300   # LP 큰 값 언더플로 방지
  x <- data.table(SNP = as.character(d[["ID"]]),
                  CHR = d[["CHROM"]], POS = d[["POS"]],
                  EA  = as.character(d[["ALT"]]), NEA = as.character(d[["REF"]]),
                  EAF = AF, BETA = ES, SE = SE, P = P, N = SS)
  rm(d, parts); gc(FALSE)
  finalize(x, trait, "continuous", info = NULL)
}

## ---- 형질별 --------------------------------------------------------------
harm_EOS <- function() read_gwasvcf("ebi-a-GCST004606.vcf.gz", "EOS", n_const = 173480)  # 호산구(Astle2016 N고정)
harm_CRP <- function() read_gwasvcf("ieu-b-35.vcf.gz",         "CRP")   # C-반응성 단백

## ---- 실행 & 저장 ---------------------------------------------------------
message("== Step 0c: 면역형질 harmonize ==")
QC_LOG <- list()
for (tn in c("EOS","CRP")) {
  message("-- ", tn)
  h <- if (tn=="EOS") harm_EOS() else harm_CRP()
  saveRDS(h, file.path(HARM, paste0(tn, ".rds")))
}

## ---- QC 리포트 -----------------------------------------------------------
qc <- rbindlist(QC_LOG); setnames(qc, c("trait","step","n_remaining"))
fwrite(qc, file.path(HARM, "_step0c_qc_report_immune.csv"))

summ <- rbindlist(lapply(c("EOS","CRP"), function(tn){
  d <- as.data.table(readRDS(file.path(HARM, paste0(tn,".rds"))))
  data.table(TRAIT=d$TRAIT[1], TYPE=d$TYPE[1], nSNP=nrow(d),
             medN=round(median(d$N)), EAF_min=round(min(d$EAF),3),
             lambdaGC=round(median(qchisq(1-d$P,1))/qchisq(0.5,1),3))
}))
message("\n== 면역형질 harmonized 요약 =="); print(summ)
message("저장: ", HARM, "  (EOS.rds, CRP.rds)")
message("점검: nSNP 규모 / lambdaGC 1.0~1.3 / EAF·BETA 방향(ES=ALT기준)")
message("== Step 0c 완료 ==")
