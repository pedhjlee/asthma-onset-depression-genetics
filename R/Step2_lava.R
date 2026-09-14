###############################################################################
## Step 2 — LAVA 국소 유전상관 (local genetic correlation)
## Project: 01_GWAS_Asthma_MDD
##
## 목적: 전역 rg가 어느 유전자좌에서 오는지(AOA×MDD/ANX), 그리고 전역 null이
##       국소 상쇄인지(COA×MDD/ANX)를 2495 LD 블록에서 국소적으로 검정.
##
## 입력: harmonized/*.rds, ldsc/ldsc_result.rds(→표본중복), ld_ref/EUR.*(1000G),
##       lava/blocks_...locfile
## 출력(lava/): 형질별 sumstats, input.info.txt, sample_overlap.txt,
##       lava_univ.csv(국소 h2), lava_bivar.csv(국소 rg, focus·FDR)
###############################################################################

## --------------------------- 0. 패키지 ----------------------------------- ##
.repo <- "https://cloud.r-project.org"
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos = .repo)
if (!requireNamespace("LAVA", quietly = TRUE)) {
  if (!requireNamespace("remotes",  quietly = TRUE)) install.packages("remotes",  repos = .repo)
  if (!requireNamespace("pkgbuild", quietly = TRUE)) install.packages("pkgbuild", repos = .repo)
  if (!pkgbuild::has_build_tools(debug = FALSE))
    stop("LAVA 설치엔 C++ 컴파일러(Rtools)가 필요합니다.\n",
         "  1) https://cran.r-project.org/bin/windows/Rtools/ 에서 R 버전에 맞는 Rtools 설치(기본옵션)\n",
         "  2) R(또는 RStudio) 완전히 껐다 다시 켜기\n",
         "  3) pkgbuild::has_build_tools()  가 TRUE인지 확인 후 이 스크립트 재실행")
  remotes::install_github("josefin-werme/LAVA")
}
suppressPackageStartupMessages({ library(LAVA); library(data.table) })

## --------------------------- 1. 경로 ------------------------------------- ##
PROJ    <- "PATH/TO/PROJECT"
HARM    <- file.path(PROJ, "harmonized")
LAVADIR <- file.path(PROJ, "lava")
REFDIR  <- file.path(PROJ, "ld_ref")
REF     <- file.path(REFDIR, "EUR")                 # plink prefix (EUR.bed/bim/fam)
LOCFILE <- file.path(LAVADIR, "blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
LDSCRDS <- file.path(PROJ, "ldsc", "ldsc_result.rds")
dir.create(LAVADIR, showWarnings = FALSE, recursive = TRUE)

## 1000G EUR plink 참조패널
##  ⚠️ 성능 핵심: LAVA는 유전자좌마다 .bed를 랜덤 액세스로 읽는데, G:(Google Drive)에
##     두면 I/O 지연으로 극도로 느림(수 시간). → 로컬 C: 임시디스크로 복사해서 사용.
LD08 <- "PATH/TO/ld_ref"
LOCAL_REF_DIR <- file.path(Sys.getenv("TEMP"), "lava_ref")   # 예: C:/Users/.../AppData/Local/Temp/lava_ref (ASCII, 로컬)
LOCAL_REF     <- file.path(LOCAL_REF_DIR, "EUR")
dir.create(LOCAL_REF_DIR, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(paste0(LOCAL_REF, ".bed"))) {
  src <- if (file.exists(paste0(REF, ".bed"))) REFDIR else LD08   # 09에 이미 있으면 거기서, 없으면 08에서
  message("1000G EUR plink → 로컬 C: 임시로 복사 중 (~1.3GB, 일회성)... 몇 분 걸림")
  for (ext in c("bed","bim","fam"))
    file.copy(file.path(src, paste0("EUR.", ext)), paste0(LOCAL_REF, ".", ext), overwrite = TRUE)
  message("  복사 완료: ", LOCAL_REF_DIR)
}
stopifnot(file.exists(paste0(LOCAL_REF, ".bed")), file.exists(LOCFILE), file.exists(LDSCRDS))

TRAITS <- c("COA","AOA","MDD","BIP","SCZ","ANX","LUNG")
ASTHMA <- c("COA","AOA")          # 관심: 천식이 포함된 페어만 bivariate 보고
UNIV_THRESH <- 0.05               # 국소 bivariate 전 단변량 사전필터(양쪽 다 통과 시)
MAXSNP      <- 25000              # 참조 SNP 이보다 많은 초대형 유전자좌는 제외(계산 폭증 방지)
ADAP_THRESH <- c(1e-100, 1e-100)  # 적응적 순열 escalation 사실상 끔(스크린 속도↑; top hit은 나중에 정밀 재계산)

## ------------------ 2. LAVA 입력 sumstats (HapMap3 제한; SNP A1 A2 N Z) --- ##
## ⚠️ 속도 핵심: 전체 1000G SNP(유전자좌당 ~3000)로 LAVA 고유분해 → R 기본 BLAS에서
##   유전자좌당 ~11초(전체 7시간+). HapMap3(유전자좌당 ~500)로 제한 → ~100배↑(전체 10~15분).
##   LDSC munge와 동일 SNP셋이라 방법론적으로 일관, 국소 rg엔 영향 미미(흔한변이 태깅).
message("\n=== [1/5] LAVA sumstats 생성 (HapMap3 제한) ===")
hm3 <- fread(file.path(PROJ, "eur_w_ld_chr", "w_hm3.snplist"))$SNP
message("  HapMap3 SNP ", format(length(hm3), big.mark = ","), "개로 제한")
for (tr in TRAITS) {
  ss <- file.path(LAVADIR, paste0(tr, ".sumstats.gz"))
  d  <- as.data.table(readRDS(file.path(HARM, paste0(tr, ".rds"))))
  d  <- d[SNP %in% hm3]                          # HapMap3 교집합
  fwrite(d[, .(SNP = SNP, A1 = EA, A2 = NEA, N = N, Z = BETA / SE)], ss, sep = "\t")
  message("  ", tr, " → ", format(nrow(d), big.mark = ","), " SNP")
}

## ------------------ 3. info 파일 + 표본중복 행렬 -------------------------- ##
message("\n=== [2/5] info 파일 + sample overlap ===")
## 전 형질 N-기반(quantitative)으로 처리 → 국소 rg는 스케일 무관하게 유효
## ⚠️ filename은 공백 없는 '상대경로'(bare)로 — LAVA read.table가 공백기준 파싱이라
##    공백이 포함된 절대경로를 넣으면 컬럼이 쪼개짐. wd=lava로 잡음.
info <- data.table(phenotype = TRAITS, cases = NA, controls = NA,
                   filename = paste0(TRAITS, ".sumstats.gz"))
INFO_FILE <- file.path(LAVADIR, "input.info.txt")
write.table(info, INFO_FILE, row.names = FALSE, quote = FALSE, sep = "\t")

## 표본중복 = cross-trait LDSC intercept 정규화(cov2cor) — LAVA 표준
L <- readRDS(LDSCRDS); I <- L$I
stopifnot(!is.null(I)); dimnames(I) <- list(TRAITS, TRAITS)
ovl <- cov2cor(I)
OVL_FILE <- file.path(LAVADIR, "sample_overlap.txt")
write.table(round(ovl, 5), OVL_FILE, quote = FALSE)

## ------------------ 4. process.input + 유전자좌 ----------------------------- ##
message("\n=== [3/5] process.input ===")
## wd=lava 로 고정 → info의 bare 파일명 + 상대 plink 경로가 resolve됨.
## ref는 상대 ASCII("../ld_ref/EUR")로 넘겨 C++ 리더가 ★/공백 문자열을 직접 파싱하지 않게 함.
setwd(LAVADIR)
input <- process.input(input.info.file    = INFO_FILE,
                       sample.overlap.file = OVL_FILE,
                       ref.prefix          = LOCAL_REF,   # 로컬 C: 절대경로(ASCII) — 빠름
                       phenos              = TRAITS)
loci <- read.loci(LOCFILE)
message("  유전자좌 ", nrow(loci), "개, 형질 ", length(TRAITS), "개")

## --- 유전자좌별 참조 SNP 수 계산 → 초대형(계산 폭증) 유전자좌 식별/제외 ---
message("  유전자좌별 참조 SNP 수 계산 중...")
bim <- fread(paste0(LOCAL_REF, ".bim"), header = FALSE, select = c(1, 4))
setnames(bim, c("bCHR", "bPOS"))
loci$nsnp_ref <- 0L
for (ch in unique(loci$CHR)) {
  p  <- sort(bim$bPOS[bim$bCHR == ch])
  ix <- which(loci$CHR == ch)
  loci$nsnp_ref[ix] <- findInterval(loci$STOP[ix], p) - findInterval(loci$START[ix] - 1L, p)
}
message(sprintf("  참조 SNP/유전자좌: 중앙 %d · 90%% %d · 최대 %d | >10k %d개 · >%d(제외) %d개",
                as.integer(median(loci$nsnp_ref)), as.integer(quantile(loci$nsnp_ref, .9)),
                max(loci$nsnp_ref), sum(loci$nsnp_ref > 10000L), MAXSNP, sum(loci$nsnp_ref > MAXSNP)))
big <- loci[loci$nsnp_ref > MAXSNP, , drop = FALSE]
if (nrow(big)) { fwrite(big, file.path(LAVADIR, "_skipped_large_loci.csv"))
  message("  초대형 유전자좌 ", nrow(big), "개 제외 → _skipped_large_loci.csv") }

## --- targeted LAVA: 천식(COA/AOA)에 신호 있는 유전자좌만 ---
## LAVA 국소 rg는 두 형질 모두 국소 h2 필요 → 천식 국소신호 있는 곳만 정보적.
## (전수 스캔은 이 머신에서 비현실적: 유전자좌당 ~6초 + 간헐 시스템 정체)
ASTHMA_P <- 1e-4                              # 천식 시사 임계(완화하려면 1e-3, 엄격 5e-8)
sig <- unique(rbindlist(lapply(ASTHMA, function(tr) {
  d <- as.data.table(readRDS(file.path(HARM, paste0(tr, ".rds"))))
  d[P < ASTHMA_P, .(CHR, POS)]
})))
loci$has_asthma <- FALSE
for (ch in unique(loci$CHR)) {
  ps <- sort(sig$POS[sig$CHR == ch]); if (!length(ps)) next
  ix <- which(loci$CHR == ch)
  loci$has_asthma[ix] <- (findInterval(loci$STOP[ix], ps) - findInterval(loci$START[ix] - 1L, ps)) > 0
}
run_idx <- which(loci$nsnp_ref <= MAXSNP & loci$has_asthma)
message(sprintf("  천식 신호(p<%.0e) 유전자좌 %d개만 분석 (전체 %d 중)", ASTHMA_P, length(run_idx), nrow(loci)))

## ------------------ 5. 유전자좌 루프 (중간저장/이어하기) ----------------- ##
## ⚠️ 이 머신은 간헐적으로 실행이 끊김(절전/스왑 추정). 50개마다 체크포인트 저장 →
##   끊겨도 다시 source하면 그 지점부터 이어감. (ASTHMA_P 등 바꿔 새로 돌릴 땐
##   _step2_checkpoint.rds 를 먼저 지울 것)
CKPT <- file.path(LAVADIR, "_step2_checkpoint.rds")
univ_list <- list(); bivar_list <- list(); done_locs <- integer(0); n_skip <- 0L
if (file.exists(CKPT)) {
  ck <- readRDS(CKPT)
  if (length(ck$univ) + length(ck$bivar) > 0) {       # 결과 있는 체크포인트만 이어감
    univ_list <- ck$univ; bivar_list <- ck$bivar; done_locs <- ck$done
    message("  체크포인트 발견 → 기존 ", length(done_locs), "개 이어서 진행")
  } else {                                             # 빈(결과0) 체크포인트는 폐기
    message("  빈 체크포인트 무시(결과 0) → 처음부터")
    file.remove(CKPT)
  }
}
todo <- run_idx[!(loci$LOC[run_idx] %in% done_locs)]
message("\n=== [4/5] 유전자좌 루프: 남은 ", length(todo), "/", length(run_idx), "개 (25단위 진행) ===")
save_ckpt <- function() saveRDS(list(done = done_locs, univ = univ_list, bivar = bivar_list), CKPT)
t0 <- proc.time()[3]; done <- 0L; t_pl <- 0; t_an <- 0
for (i in todo) {
  ok <- tryCatch({
    locus <- NULL; res <- NULL; tp <- 0; ta <- 0
    invisible(capture.output(suppressWarnings(suppressMessages({
      tp <- system.time(locus <- process.locus(loci[i, ], input))[3]
      if (!is.null(locus))   # target 안 씀(단일형질만 허용) → 전 페어 계산 후 focus로 추림
        ta <- system.time(res <- run.univ.bivar(locus, univ.thresh = UNIV_THRESH,
                                                adap.thresh = ADAP_THRESH))[3]
    }))))
    t_pl <- t_pl + tp; t_an <- t_an + ta
    if (!is.null(locus)) {
      lb <- data.table(LOC = loci$LOC[i], CHR = loci$CHR[i],
                       START = loci$START[i], STOP = loci$STOP[i],
                       N_SNP = tryCatch(locus$n.snps, error = function(e) NA_integer_))
      if (!is.null(res$univ)  && nrow(res$univ))
        univ_list[[length(univ_list)+1]]  <- cbind(lb, as.data.table(res$univ))
      if (!is.null(res$bivar) && nrow(res$bivar))
        bivar_list[[length(bivar_list)+1]] <- cbind(lb, as.data.table(res$bivar))
    }
    TRUE
  }, error = function(e) FALSE)
  if (!ok) n_skip <- n_skip + 1L
  done_locs <- c(done_locs, loci$LOC[i])       # 처리기록(성공/실패 무관 → 재시도 안함)
  done <- done + 1L
  if (done %% 25 == 0) {
    el <- proc.time()[3] - t0
    message(sprintf("  ...%d/%d (loc%d) %.0f초, %.2f초/loc [locus %.2f + 분석 %.2f] skip%d | univ%d bivar%d",
                    done, length(todo), loci$LOC[i], el, el/done, t_pl/done, t_an/done, n_skip,
                    length(univ_list), length(bivar_list)))
    flush.console()
  }
  if (done %% 50 == 0) { save_ckpt(); gc(FALSE) }   # 50개마다 중간저장 + 메모리정리
}
save_ckpt()
message("  완료. 처리불가 유전자좌(누적): ", n_skip)

## ------------------ 6. 정리·저장 ---------------------------------------- ##
message("\n=== [5/5] 결과 정리 ===")
univ  <- rbindlist(univ_list,  fill = TRUE)
bivar <- rbindlist(bivar_list, fill = TRUE)
if (nrow(univ))  fwrite(univ,  file.path(LAVADIR, "lava_univ.csv"))
stopifnot(nrow(bivar) > 0)   # 천식 포함 국소 rg가 하나도 없으면 설정 재점검

## target=ASTHMA 라 모든 bivar 행이 천식 포함(focus=TRUE)
bivar[, focus := (phen1 %in% ASTHMA) | (phen2 %in% ASTHMA)]
bivar[, pair := paste(phen1, phen2, sep = "_")]
## 다중검정: focus(천식 포함) bivariate 검정 전체에 BH-FDR + Bonferroni
bivar[focus == TRUE, p.FDR := p.adjust(p, "BH")]
n_focus_tests <- bivar[focus == TRUE, .N]
bivar[focus == TRUE, p.Bonf := pmin(1, p * n_focus_tests)]
setorder(bivar, -focus, p)
fwrite(bivar, file.path(LAVADIR, "lava_bivar.csv"))

message("\n[요약] 국소 rg 검정 수(천식 포함): ", n_focus_tests,
        "  | Bonferroni 임계 p<", signif(0.05 / max(n_focus_tests,1), 3))
message("[천식 포함 유의 국소 rg (FDR<0.05)]")
sig <- bivar[focus == TRUE & p.FDR < 0.05][order(p)]
if (nrow(sig)) print(sig[, .(LOC, CHR, START, STOP, phen1, phen2,
                             rho = round(rho,3), rho.lower = round(rho.lower,3),
                             rho.upper = round(rho.upper,3), r2 = round(r2,3),
                             p = signif(p,3), p.FDR = signif(p.FDR,3))]) else
  message("  (FDR<0.05 없음 — nominal p<0.05는 lava_bivar.csv 참조)")

message("\n================= Step 2 완료 =================")
message("저장: ", LAVADIR)
message("  lava_univ.csv (국소 h2), lava_bivar.csv (국소 rg; focus=천식포함, p.FDR/p.Bonf)")
message("점검: AOA×MDD·AOA×ANX 국소 rg가 특정 locus에 집중되는지; COA×정신질환이 국소 +/− 상쇄인지")
writeLines(capture.output(sessionInfo()), file.path(LAVADIR, "_step2_sessionInfo.txt"))
if (file.exists(CKPT)) file.remove(CKPT)   # 완주 → 체크포인트 삭제(다음 실행은 새로)
