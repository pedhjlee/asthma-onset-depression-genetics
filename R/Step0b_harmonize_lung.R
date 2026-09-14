###############################################################################
## Step 0b — FEV1 · FVC harmonize  (폐기능 메인 축 추가)
## Project: 01_GWAS_Asthma_MDD
##
## 배경 / 설계 근거
##   폐기능 축을 FEV1/FVC ratio(폐쇄지표) 대신 FEV1·FVC(볼륨지표) 2개로 확장.
##   - 관찰연구 근거: 우울은 폐쇄(ratio↓)보다 볼륨감소(FVC↓, 제한성/PRISm 패턴)와
##     더 일관되게 연관 (UK Biobank PRISm: 우울 HR 1.12 / BMC Med 280k: FVC·FEV1↓ → 우울,
##     염증 부분매개). 천식은 폐쇄(FEV1↓)와 연관.
##   - 따라서 FEV1 = 천식·기류제한 축, FVC = 우울·제한성/전신 축 으로 판별(discriminant) 대비.
##   - FEV1/FVC ratio(LUNG.rds)는 예비 보관(리뷰어 대응용).
##
## 방식
##   기존 Step0_harmonize.R 의 finalize()/read_raw()/QC 파라미터를 그대로 재사용
##   (STEP0_NORUN <- TRUE 로 source → 실행부는 건너뛰고 함수·경로만 로드).
##   SpiroMeta FEV1·FVC 는 ratio 와 헤더·포맷 100% 동일하므로 harm_LUNG 을 복제.
##   헤더: #SNP Chromosome Position_b37 Coded Non_coded N Neff Coded_freq beta SE P
##   좌표 GRCh37, beta 단위 = SD (rank inverse-normal transform), 연속형(실제 N 사용).
##
## 산출: harmonized/FEV1.rds, harmonized/FVC.rds
##       + _step0b_qc_report_lung.csv, _step0_summary.csv(10형질로 갱신)
##
## 실행: R에서   source("Step0b_harmonize_lung.R")
##       (Step0_harmonize.R 와 같은 "R script" 폴더에 두고 실행)
###############################################################################

## ------------------------- 0. 기존 Step0 함수 로드 ------------------------ ##
## 실행부는 돌리지 않고 finalize()/read_raw()/RAW/HARM/QC파라미터만 가져온다.
STEP0_NORUN <- TRUE
SRC <- "PATH/TO/PROJECT/R script/Step0_harmonize.R"
if (!file.exists(SRC)) SRC <- "Step0_harmonize.R"   # 같은 폴더에서 실행하는 경우 fallback
source(SRC)
stopifnot(exists("finalize"), exists("read_raw"), exists("RAW"), exists("HARM"))

## ------------------------- 1. FEV1 · FVC harmonize ------------------------ ##
## harm_LUNG 과 동일 로직 — 파일명·형질명만 교체.
harm_FEV1 <- function() {
  d <- read_raw("Shrine_30804560_SpiroMeta_FEV1.txt.gz", skip = 0)
  setnames(d, 1, "SNP")                       # 1번 컬럼 = rsID (#SNP)
  x <- d[, .(SNP, CHR = Chromosome, POS = Position_b37,
             EA = Coded, NEA = Non_coded, EAF = Coded_freq,
             BETA = beta, SE = SE, P = P, N = N)]   # 연속형: 실제 N
  finalize(x, "FEV1", "continuous", info = NULL)
}

harm_FVC <- function() {
  d <- read_raw("Shrine_30804560_SpiroMeta_FVC.txt.gz", skip = 0)
  setnames(d, 1, "SNP")
  x <- d[, .(SNP, CHR = Chromosome, POS = Position_b37,
             EA = Coded, NEA = Non_coded, EAF = Coded_freq,
             BETA = beta, SE = SE, P = P, N = N)]
  finalize(x, "FVC", "continuous", info = NULL)
}

## ------------------------- 2. 실행 & 저장 --------------------------------- ##
message("== Step 0b: FEV1 · FVC harmonize 시작 ==")
QC_LOG <- list()                                # 이번 실행분 QC 로그 초기화
for (tn in c("FEV1", "FVC")) {
  message("-- ", tn)
  h <- if (tn == "FEV1") harm_FEV1() else harm_FVC()
  saveRDS(h, file.path(HARM, paste0(tn, ".rds")))
}

## ------------------------- 3. QC 리포트 ----------------------------------- ##
qc <- rbindlist(QC_LOG); setnames(qc, c("trait", "step", "n_remaining"))
fwrite(qc, file.path(HARM, "_step0b_qc_report_lung.csv"))

## 갱신된 전체 요약 (이제 FEV1·FVC 포함 최대 10형질)
summ <- rbindlist(lapply(list.files(HARM, pattern = "\\.rds$", full.names = TRUE), function(f) {
  d <- as.data.table(readRDS(f))
  data.table(TRAIT = d$TRAIT[1], TYPE = d$TYPE[1], nSNP = nrow(d),
             medN = round(median(d$N)), EAF_min = round(min(d$EAF), 3),
             lambdaGC = round(median(qchisq(1 - d$P, 1)) / qchisq(0.5, 1), 3))
}))
setorder(summ, TYPE, TRAIT)
fwrite(summ, file.path(HARM, "_step0_summary.csv"))
message("\n== 갱신 harmonized 요약 (FEV1·FVC 포함) =="); print(summ)
message("\n저장: ", HARM, "  (FEV1.rds, FVC.rds)")
message("점검: FEV1·FVC 의 nSNP 가 LUNG(ratio) 과 비슷한지 / lambdaGC 가 1.0~1.2 수준인지")
message("== Step 0b 완료 ==")
