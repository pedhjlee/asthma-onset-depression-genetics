###############################################################################
## Step 0 — Harmonize GWAS summary statistics to a common schema
## Project: 01_GWAS_Asthma_MDD  (호흡곤란 축 × 천식-정신질환 공유유전)
##
## 목표: 새 6개 형질(MDD, BIP, SCZ, ANX, LUNG)을 기존 천식(COA/AOA) harmonized와
##       동일한 12컬럼 스키마로 통일.
##   SNP CHR POS EA NEA EAF BETA SE P N TRAIT TYPE
##   - EA = effect allele, BETA = EA 기준 효과 (binary: log-odds / continuous: SD)
##   - N  = effective N (binary), 실제 N (continuous)
##   - 좌표: GRCh37 (모든 파일 확인됨)
##
## 실행: R에서 source("Step0_harmonize.R") 또는 Rscript Step0_harmonize.R
###############################################################################

## ----------------------------- 0. 설정 ------------------------------------ ##
PROJ    <- "PATH/TO/PROJECT"
RAW     <- file.path(PROJ, "raw")
HARM    <- file.path(PROJ, "harmonized")
dir.create(HARM, showWarnings = FALSE, recursive = TRUE)

## QC 임계값 (필요시 조정)
INFO_MIN <- 0.90     # imputation INFO 최소 (INFO 컬럼 있는 형질만 적용)
MAF_MIN  <- 0.01     # minor allele frequency 최소
ANX_N_FRAC <- 0.70   # ANX: per-SNP N >= 0.70 * max(N) 만 유지 (프로토콜)
RS_ONLY  <- TRUE     # rsID(^rs) SNP만 유지 (reference 매칭·중복 방지)
UPDATE_ASTHMA_NEFF <- TRUE  # 천식 COA/AOA의 N을 effective N으로 갱신

## 천식 effective N (Ferreira 2019, PMID 30929738)
##   COA 13,962 cases / 300,671 controls -> Neff = 4*ca*co/(ca+co)
##   AOA 26,582 cases / 300,671 controls
NEFF_COA <- 4 * 13962 * 300671 / (13962 + 300671)   # ≈ 53,371
NEFF_AOA <- 4 * 26582 * 300671 / (26582 + 300671)   # ≈ 97,692

## ------------------------- 1. 패키지 -------------------------------------- ##
.repo <- "https://cloud.r-project.org"
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos = .repo)
if (!requireNamespace("R.utils",    quietly = TRUE)) install.packages("R.utils",    repos = .repo)  # fread gz 지원
suppressPackageStartupMessages(library(data.table))
setDTthreads(0)  # 모든 코어

## ------------------------- 2. 공용 함수 ----------------------------------- ##
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

## case/control 빈도 -> 결합 EAF (샘플수 가중)
combined_eaf <- function(fca, fco, nca, nco) {
  (nca * fca + nco * fco) / (nca + nco)
}

## 표준 QC + 스키마 정렬. 입력 dt는 이미 캐논 컬럼(SNP,CHR,POS,EA,NEA,EAF,BETA,SE,P,N) 보유
finalize <- function(dt, trait, type, info = NULL, n_frac = NA_real_) {
  n0 <- nrow(dt)
  step <- function(msg, n) QC_LOG[[length(QC_LOG)+1]] <<- data.table(trait = trait, step = msg, n = n)
  step("raw rows", n0)

  ## 대립유전자 대문자화, ACGT 단일염기(SNP)만
  dt[, `:=`(EA = toupper(EA), NEA = toupper(NEA))]
  dt <- dt[EA %in% c("A","C","G","T") & NEA %in% c("A","C","G","T") & EA != NEA]
  step("ACGT 단일염기 SNP만", nrow(dt))

  ## 상염색체 1-22
  dt[, CHR := suppressWarnings(as.integer(CHR))]
  dt <- dt[CHR %in% 1:22]
  dt[, POS := suppressWarnings(as.integer(POS))]
  dt <- dt[!is.na(POS)]
  step("상염색체 1-22 & POS 정수", nrow(dt))

  ## 값 유효성
  dt <- dt[is.finite(BETA) & is.finite(SE) & SE > 0 &
             is.finite(P) & P > 0 & P <= 1 & is.finite(EAF) & is.finite(N) & N > 0]
  step("BETA/SE/P/EAF/N 유효값", nrow(dt))

  ## INFO 필터 (있을 때만)
  if (!is.null(info)) { dt <- dt[is.finite(info_tmp) & info_tmp >= INFO_MIN]; step(paste0("INFO>=",INFO_MIN), nrow(dt)) }

  ## MAF
  dt <- dt[EAF >= MAF_MIN & EAF <= (1 - MAF_MIN)]
  step(paste0("MAF>=",MAF_MIN), nrow(dt))

  ## per-SNP N 필터 (ANX 등)
  if (!is.na(n_frac)) { thr <- n_frac * max(dt$N, na.rm = TRUE)
    dt <- dt[N >= thr]; step(sprintf("N>=%.2f*max(N)=%.0f", n_frac, thr), nrow(dt)) }

  ## rsID만
  if (RS_ONLY) { dt <- dt[grepl("^rs", SNP)]; step("rsID(^rs)만", nrow(dt)) }

  ## 중복 rsID 제거 (양쪽 다 제거)
  dup <- dt$SNP[duplicated(dt$SNP)]
  if (length(dup)) dt <- dt[!SNP %in% dup]
  step("중복 rsID 제거", nrow(dt))

  ## 스키마 정렬
  dt[, `:=`(TRAIT = trait, TYPE = type)]
  dt <- dt[, .(SNP = as.character(SNP), CHR = as.integer(CHR), POS = as.integer(POS),
               EA = as.character(EA), NEA = as.character(NEA),
               EAF = as.numeric(EAF), BETA = as.numeric(BETA), SE = as.numeric(SE),
               P = as.numeric(P), N = as.numeric(N),
               TRAIT = as.character(TRAIT), TYPE = as.character(TYPE))]

  ## QC 요약 지표
  lambda <- median(qchisq(1 - dt$P, df = 1), na.rm = TRUE) / qchisq(0.5, 1)
  top    <- dt[which.min(P)]
  message(sprintf("  [%s] %d -> %d SNP | λGC=%.3f | top: %s chr%d:%d P=%.2e",
                  trait, n0, nrow(dt), lambda, top$SNP, top$CHR, top$POS, top$P))
  attr(dt, "lambda") <- lambda
  dt
}

## reader: gz 파일을 skip 이후 header=TRUE 로 읽기
read_raw <- function(file, skip) fread(file.path(RAW, file), skip = skip, header = TRUE,
                                       sep = "\t", fill = TRUE,   # 컬럼 모자란 줄(indel 등)에서 멈추지 않고 NA로 채워 끝까지
                                       showProgress = FALSE, data.table = TRUE)

QC_LOG <- list()

## ------------------------- 3. 형질별 harmonize ---------------------------- ##

## ---- MDD (PGC-MDD2025, noUKBB) -------------------------------------------
## 헤더(54 meta): #CHROM POS ID EA NEA BETA SE PVAL FCAS FCON IMPINFO NEFF NCAS NCON HETI HETDF HETPVAL
harm_MDD <- function() {
  d <- read_raw("pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv.gz", skip = "#CHROM")
  setnames(d, 1, "CHR")          # 1번 컬럼 = 염색체 (#CHROM)
  x <- d[, .(SNP = ID, CHR, POS, EA, NEA,
             EAF  = combined_eaf(FCAS, FCON, NCAS, NCON),
             BETA = BETA, SE = SE, P = PVAL,
             N    = NEFF,            # NEFF = effective N (직접)
             info_tmp = IMPINFO)]
  finalize(x, "MDD", "binary", info = x$info_tmp)
}

## ---- BIP (PGC3 BIP, noUKBB, daner) ---------------------------------------
## 헤더: CHR SNP BP A1 A2 FRQ_A_40463 FRQ_U_313436 INFO OR SE P ngt Direction ... Nca Nco Neff_half
harm_BIP <- function() {
  d <- read_raw("daner_bip_pgc3_nm_noukbiobank.gz", skip = 0)
  fa <- grep("^FRQ_A", names(d), value = TRUE)[1]; fu <- grep("^FRQ_U", names(d), value = TRUE)[1]
  x <- d[, .(SNP, CHR, POS = BP, EA = A1, NEA = A2,
             EAF  = combined_eaf(get(fa), get(fu), Nca, Nco),
             BETA = log(OR), SE = SE, P = P,
             N    = 2 * Neff_half,   # daner: effective N = 2*Neff_half
             info_tmp = INFO)]
  finalize(x, "BIP", "binary", info = x$info_tmp)
}

## ---- SCZ (PGC3 SCZ wave3, public) ----------------------------------------
## 헤더(73 meta): CHROM ID POS A1 A2 FCAS FCON IMPINFO BETA SE PVAL NCAS NCON NEFF
harm_SCZ <- function() {
  d <- read_raw("PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz", skip = "CHROM\tID")
  x <- d[, .(SNP = ID, CHR = CHROM, POS, EA = A1, NEA = A2,
             EAF  = combined_eaf(FCAS, FCON, NCAS, NCON),
             BETA = BETA, SE = SE, P = PVAL,
             N    = NEFF,            # NEFF = effective N (직접)
             info_tmp = IMPINFO)]
  finalize(x, "SCZ", "binary", info = x$info_tmp)
}

## ---- ANX (2026 fullANX v12, woUTAH, daner; UKB 포함) ----------------------
## 헤더: CHR SNP BP A1 A2 FRQ_A_122083 FRQ_U_729602 INFO OR SE P ... Nca Nco Neff_half
harm_ANX <- function() {
  d <- read_raw("ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz", skip = 0)
  fa <- grep("^FRQ_A", names(d), value = TRUE)[1]; fu <- grep("^FRQ_U", names(d), value = TRUE)[1]
  x <- d[, .(SNP, CHR, POS = BP, EA = A1, NEA = A2,
             EAF  = combined_eaf(get(fa), get(fu), Nca, Nco),
             BETA = log(OR), SE = SE, P = P,
             N    = 2 * Neff_half,
             info_tmp = INFO)]
  finalize(x, "ANX", "binary", info = x$info_tmp, n_frac = ANX_N_FRAC)
}

## ---- LUNG (Shrine 2019 SpiroMeta FEV1/FVC; 연속형 대조) --------------------
## 헤더: #SNP Chromosome Position_b37 Coded Non_coded N Neff Coded_freq beta SE P
harm_LUNG <- function() {
  d <- read_raw("Shrine_30804560_SpiroMeta_FEV1_to_FVC_RATIO.txt.gz", skip = 0)
  setnames(d, 1, "SNP")          # 1번 컬럼 = rsID (#SNP)
  x <- d[, .(SNP, CHR = Chromosome, POS = Position_b37, EA = Coded, NEA = Non_coded,
             EAF  = Coded_freq, BETA = beta, SE = SE, P = P,
             N    = N)]           # 연속형: 실제 N
  finalize(x, "LUNG", "continuous", info = NULL)
}

## ------------------------- 4. 실행 & 저장 --------------------------------- ##
## STEP0_NORUN 가 정의돼 있으면 함수만 로드하고 실행부는 건너뜀
## (일부 형질만 다시 돌릴 때: STEP0_NORUN<-TRUE; source(...); saveRDS(harm_BIP(), file.path(HARM,"BIP.rds")))
if (!exists("STEP0_NORUN")) {
message("== Step 0 harmonize 시작 ==")
traits <- list(MDD = harm_MDD, BIP = harm_BIP, SCZ = harm_SCZ, ANX = harm_ANX, LUNG = harm_LUNG)
lambdas <- list()
for (tn in names(traits)) {
  message("-- ", tn)
  h <- traits[[tn]]()
  saveRDS(h, file.path(HARM, paste0(tn, ".rds")))
  lambdas[[tn]] <- attr(h, "lambda")
}

## ---- 천식 N을 effective N으로 갱신 (선택) --------------------------------
if (UPDATE_ASTHMA_NEFF) {
  for (a in c("COA", "AOA")) {
    f <- file.path(HARM, paste0(a, ".rds"))
    if (file.exists(f)) {
      z <- as.data.table(readRDS(f))
      z[, N := if (a == "COA") NEFF_COA else NEFF_AOA]
      saveRDS(z, f)
      message(sprintf("  [%s] N -> effective %.0f", a, z$N[1]))
    }
  }
}

## ------------------------- 5. QC 리포트 ----------------------------------- ##
qc <- rbindlist(QC_LOG)
setnames(qc, c("trait", "step", "n_remaining"))
fwrite(qc, file.path(HARM, "_step0_qc_report.csv"))

## 최종 요약 테이블 (8개 형질)
summ <- rbindlist(lapply(list.files(HARM, pattern = "\\.rds$", full.names = TRUE), function(f) {
  d <- as.data.table(readRDS(f))
  data.table(TRAIT = d$TRAIT[1], TYPE = d$TYPE[1], nSNP = nrow(d),
             medN = round(median(d$N)), EAF_min = round(min(d$EAF), 3),
             lambdaGC = round(median(qchisq(1 - d$P, 1)) / qchisq(0.5, 1), 3))
}))
setorder(summ, TYPE, TRAIT)
fwrite(summ, file.path(HARM, "_step0_summary.csv"))
message("\n== 최종 harmonized 요약 =="); print(summ)
message("\n저장 위치: ", HARM)
message("리포트: _step0_qc_report.csv, _step0_summary.csv")
sessionInfo_txt <- capture.output(sessionInfo())
writeLines(sessionInfo_txt, file.path(HARM, "_step0_sessionInfo.txt"))
message("== Step 0 완료 ==")
}  # end if(!exists("STEP0_NORUN"))
