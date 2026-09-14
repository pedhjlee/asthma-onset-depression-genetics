###############################################################################
## Step 1 — LDSC 전역 유전상관(rg) + 표본중복(intercept) 행렬
## Project: 01_GWAS_Asthma_MDD  (호흡곤란 축 × 천식-정신질환 공유유전)
##
## 입력: harmonized/*.rds (Step 0 + 0b 결과, 최대 10형질: 정신4·천식2·폐기능3·... )
## 방법: GenomicSEM munge + ldsc (08 파이프라인 계승)
## 출력(ldsc/):
##   - ldsc_result.rds        : LDSCoutput 전체(S, V, I, N, m)
##   - h2_results.csv         : 형질별 h2(관측/liability) + intercept
##   - rg_matrix.csv          : 7x7 rg 점추정 행렬
##   - rg_results.csv         : 전 페어 rg + SE + P + FDR (관심페어 표시)
##   - intercept_matrix.csv   : LDSC intercept 행렬 I (대각=단변량, 비대각=교차)
##   - sample_overlap.csv     : LAVA용 표본중복 행렬(=I) + cov2cor(I)
##
## ⚠️ effective N 사용 → binary는 sample.prev=0.5 (liability 변환용).
##    rg는 스케일 무관하므로 유병률 값에 영향받지 않음.
###############################################################################

suppressPackageStartupMessages({ library(GenomicSEM); library(data.table) })

## --------------------------- 0. 경로 & 파라미터 --------------------------- ##
PROJ  <- "PATH/TO/PROJECT"
HARM  <- file.path(PROJ, "harmonized")
LDSC  <- file.path(PROJ, "ldsc")
dir.create(LDSC, showWarnings = FALSE, recursive = TRUE)

## LD reference (08에서 09로 복사; 이미 있으면 건너뜀) — 프로젝트 자체 보관용
LD08  <- "PATH/TO/eur_w_ld_chr"
LDDIR <- file.path(PROJ, "eur_w_ld_chr")
HM3   <- file.path(LDDIR, "w_hm3.snplist")
if (!file.exists(HM3)) {
  message("LD reference 09로 복사 중 (eur_w_ld_chr)...")
  file.copy(LD08, PROJ, recursive = TRUE)
}
stopifnot(file.exists(HM3))

REMUNGE <- FALSE   # TRUE면 기존 .sumstats.gz 있어도 다시 munge

## 형질 순서 (binary 먼저, 연속형 폐기능 3종 마지막)
## 폐기능: FEV1·FVC = 메인 축(볼륨), LUNG = FEV1/FVC ratio(폐쇄) 예비대조
TRAITS <- c("COA", "AOA", "MDD", "BIP", "SCZ", "ANX", "FEV1", "FVC", "LUNG")
TYPE   <- c(COA="binary", AOA="binary", MDD="binary", BIP="binary",
            SCZ="binary", ANX="binary",
            FEV1="continuous", FVC="continuous", LUNG="continuous")

## liability-scale 모집단 유병률 (⚠️ 문헌값 — 사용자 확정 필요; rg엔 무관)
##   effective N 사용 → sample.prev = 0.5 (binary), 연속형은 NA
POP_PREV <- c(COA=0.10, AOA=0.05, MDD=0.15, BIP=0.01, SCZ=0.01, ANX=0.16,
              FEV1=NA, FVC=NA, LUNG=NA)
samp.prev <- ifelse(TYPE[TRAITS] == "binary", 0.5, NA)
pop.prev  <- POP_PREV[TRAITS]

setwd(LDSC)  # munge는 wd에 .sumstats.gz 기록

## ----------------- 1. .rds → LDSC 입력 txt (SNP,A1,A2,effect,P,N) --------- ##
message("\n=== [1/4] .rds → txt ===")
for (tr in TRAITS) {
  d <- as.data.table(readRDS(file.path(HARM, paste0(tr, ".rds"))))
  stopifnot(all(c("SNP","EA","NEA","BETA","P","N") %in% names(d)))
  fwrite(d[, .(SNP = SNP, A1 = EA, A2 = NEA, effect = BETA, P = P, N = N)],
         file.path(LDSC, paste0(tr, ".txt")), sep = "\t", quote = FALSE)
  message("  ", tr, ": ", format(nrow(d), big.mark = ","), " SNP (N중앙값 ",
          format(round(median(d$N)), big.mark = ","), ")")
}

## ------------------------ 2. munge (HapMap3 표준화) ----------------------- ##
message("\n=== [2/4] munge ===")
need <- if (REMUNGE) TRAITS else TRAITS[!file.exists(file.path(LDSC, paste0(TRAITS, ".sumstats.gz")))]
if (length(need)) {
  munge(files       = file.path(LDSC, paste0(need, ".txt")),
        hm3         = HM3,
        trait.names = need,
        info.filter = 0.9,     # INFO 컬럼 없으면 자동 무시
        maf.filter  = 0.01)    # MAF 컬럼 없으면 자동 무시(Step0에서 이미 필터함)
} else message("  모든 .sumstats.gz 존재 → munge 건너뜀 (REMUNGE=TRUE로 강제)")

## ------------------------ 3. ldsc (rg + intercept) ----------------------- ##
message("\n=== [3/4] ldsc ===")
LDSCoutput <- ldsc(
  traits          = file.path(LDSC, paste0(TRAITS, ".sumstats.gz")),
  sample.prev     = samp.prev,
  population.prev = pop.prev,
  ld              = LDDIR,
  wld             = LDDIR,
  trait.names     = TRAITS
)
saveRDS(LDSCoutput, file.path(LDSC, "ldsc_result.rds"))

S <- LDSCoutput$S; V <- LDSCoutput$V; I <- LDSCoutput$I; k <- nrow(S)
rownames(S) <- colnames(S) <- TRAITS
if (!is.null(I)) rownames(I) <- colnames(I) <- TRAITS

## vech(하삼각, column-major) 인덱스 — GenomicSEM V 순서
vech_idx <- matrix(NA_integer_, k, k); cc <- 1L
for (j in 1:k) for (i in j:k) { vech_idx[i,j] <- cc; vech_idx[j,i] <- cc; cc <- cc + 1L }

## ------------------------ 4. 결과 정리 ------------------------------------ ##
message("\n=== [4/4] 결과 정리 ===")

## (a) h2 + intercept
h2    <- diag(S)
h2_se <- sqrt(sapply(1:k, function(i) V[vech_idx[i,i], vech_idx[i,i]]))
h2_tbl <- data.table(
  TRAIT = TRAITS, SCALE = ifelse(TYPE[TRAITS]=="binary","liability","obs"),
  h2 = round(h2,4), h2_SE = round(h2_se,4), h2_Z = round(h2/h2_se,2),
  intercept = if (!is.null(I)) round(diag(I),4) else NA_real_
)
message("\n[h2 + 단변량 intercept]"); print(h2_tbl)
fwrite(h2_tbl, file.path(LDSC, "h2_results.csv"))

## (b) rg 행렬 (점추정)
RG <- cov2cor(S); rownames(RG) <- colnames(RG) <- TRAITS
message("\n[rg 행렬]"); print(round(RG,3))
fwrite(data.table(TRAIT=TRAITS, round(as.data.table(RG),3)), file.path(LDSC, "rg_matrix.csv"))

## (c) 전 페어 rg + delta-method SE + P + FDR
rg_se <- function(i, j) {
  a <- as.numeric(S[i,i]); b <- as.numeric(S[j,j]); rg <- as.numeric(RG[i,j])
  ia <- vech_idx[i,i]; ib <- vech_idx[j,j]; ic <- vech_idx[i,j]
  Sig <- matrix(c(V[ic,ic],V[ic,ia],V[ic,ib],
                  V[ia,ic],V[ia,ia],V[ia,ib],
                  V[ib,ic],V[ib,ia],V[ib,ib]), 3, 3)
  grad <- c(1/sqrt(a*b), -rg/(2*a), -rg/(2*b))
  as.numeric(sqrt(t(grad) %*% Sig %*% grad))
}
cmb <- t(combn(k, 2))
rg_all <- rbindlist(lapply(1:nrow(cmb), function(r){
  i <- cmb[r,1]; j <- cmb[r,2]
  rg <- as.numeric(RG[i,j]); se <- rg_se(i,j); z <- rg/se
  gcov_int <- if (!is.null(I)) round(I[i,j],4) else NA_real_
  data.table(T1=TRAITS[i], T2=TRAITS[j], rg=round(rg,4), SE=round(se,4),
             Z=round(z,3), P=signif(2*pnorm(-abs(z)),3), gcov_intercept=gcov_int)
}))
rg_all[, P_FDR := signif(p.adjust(P, "BH"), 3)]

## 관심 페어 표시: 천식(COA/AOA) × 정신·폐기능, 그리고 양성대조 COA×AOA
asthma <- c("COA","AOA"); psy <- c("MDD","ANX","BIP","SCZ"); ctrl <- c("FEV1","FVC","LUNG")
rg_all[, focus := fifelse(
  (T1 %in% asthma & T2 %in% c(psy,ctrl)) | (T2 %in% asthma & T1 %in% c(psy,ctrl)), "asthma_x_trait",
  fifelse((T1=="COA"&T2=="AOA")|(T1=="AOA"&T2=="COA"), "pos_control_COAxAOA", ""))]
setorder(rg_all, -focus, P)
message("\n[전 페어 rg]"); print(rg_all)
fwrite(rg_all, file.path(LDSC, "rg_results.csv"))

## (d) 표본중복 행렬 (LAVA 입력)
if (!is.null(I)) {
  fwrite(data.table(TRAIT=TRAITS, round(as.data.table(I),5)), file.path(LDSC, "intercept_matrix.csv"))
  I_cor <- cov2cor(I); rownames(I_cor) <- colnames(I_cor) <- TRAITS
  fwrite(data.table(TRAIT=TRAITS, round(as.data.table(I_cor),4)), file.path(LDSC, "sample_overlap.csv"))
  message("\n[표본중복(교차 intercept) — 0에서 멀수록 overlap 큼]")
  print(round(I,4))
  message("→ 천식×불안(ANX) 교차 intercept가 유의하게 0이 아니면 overlap 확인됨(예상). LAVA에 I 사용.")
}

message("\n================= Step 1 완료 =================")
message("저장: ", LDSC)
message("핵심 산출: rg_results.csv(focus 페어), sample_overlap.csv(LAVA용 intercept 행렬)")
message("점검: (1) COA×AOA rg 높음(양성대조) (2) h2 Z>4 (3) 천식×ANX gcov_intercept≠0")
writeLines(capture.output(sessionInfo()), file.path(LDSC, "_step1_sessionInfo.txt"))
