###############################################################################
## Step 3 — PLACO+ 다면발현(pleiotropy) 스캔
## Project: 01_GWAS_Asthma_MDD
##
## 천식(COA/AOA) × 정신질환/폐기능 각 페어에서 두 형질에 동시 연관된
## 다면발현 SNP을 PLACO+로 스캔(표본중복/형질상관 CorZ 보정).
##
## 이 스크립트 = 발견(discovery): 사전스크리닝(SCREEN_P) 켜서 빠르게.
##   → 최종 논문 rigor는 순열검증(아래 RUN_PERM=TRUE, SCREEN_P<-1) 별도 실행.
##
## 08 교훈 반영: 병렬(PSOCK), 워커에 PLACO source, 사전검증(NA체크),
##   페어별 체크포인트(재실행 시 이어감), lambda 캘리브레이션.
## ⚠️ 실행 중 노트북 절전 끄기(중간에 끊기면 안 끝난 페어만 재실행).
###############################################################################

suppressMessages({ library(data.table); library(parallel) })

## --------------------------- 경로 --------------------------------------- ##
PROJ  <- "PATH/TO/PROJECT"
HARM  <- file.path(PROJ, "harmonized")
PLACO <- file.path(PROJ, "placo")
dir.create(PLACO, showWarnings = FALSE, recursive = TRUE)
PLACO_SRC <- file.path(PLACO, "PLACO_v0.2.0.R")
stopifnot(file.exists(PLACO_SRC))
source(PLACO_SRC)
stopifnot(all(c("var.placo","cor.pearson","placo","placo.plus") %in% ls()))

## --------------------------- 설정 --------------------------------------- ##
EXPOS  <- c("COA", "AOA")                                      # 천식
OUTC   <- c("MDD", "ANX", "BIP", "SCZ", "FEV1", "FVC", "LUNG") # 정신질환 + 폐기능(FEV1·FVC 메인, LUNG=ratio 대조)
## 체크포인트(SUMM_*.rds)라 기존 페어는 건너뛰고 천식×FEV1·FVC 4페어만 새로 실행됨
SCREEN_P <- 1e-3        # 발견용 사전스크리닝. 최종검증은 1 (전체SNP)
Z2_MAX   <- 80          # 저자 권장
VAR_P    <- 1e-4        # var.placo / cor.pearson p.threshold
GW       <- 5e-8
N_CORE   <- max(1, detectCores() - 1)
CHUNK    <- 5000

## --------------------------- 공통 함수 ---------------------------------- ##
load_pair <- function(ex, ou) {
  t1 <- as.data.table(readRDS(file.path(HARM, paste0(ex, ".rds"))))
  t2 <- as.data.table(readRDS(file.path(HARM, paste0(ou, ".rds"))))
  t1[, Z := BETA/SE]; t2[, Z := BETA/SE]
  t1 <- t1[is.finite(Z) & !is.na(P)]; t2 <- t2[is.finite(Z) & !is.na(P)]
  cm <- merge(t1[, .(SNP, CHR, POS, EA1=EA, NEA1=NEA, Z1=Z, P1=P)],
              t2[, .(SNP, EA2=EA, NEA2=NEA, Z2=Z, P2=P)], by="SNP")
  rm(t1, t2); gc(FALSE)
  cm[, mt := fifelse(EA1==EA2 & NEA1==NEA2, "same",
             fifelse(EA1==NEA2 & NEA1==EA2, "flip", "mismatch"))]
  cm <- cm[mt != "mismatch"]; cm[mt=="flip", Z2 := -Z2]
  cm <- cm[Z1^2 <= Z2_MAX & Z2^2 <= Z2_MAX]
  cm[, mt := NULL]; cm[]
}

count_loci <- function(sig, window_kb = 500) {
  if (!nrow(sig)) return(0L)
  sig <- sig[order(P_PLACO)]; taken <- rep(FALSE, nrow(sig)); n <- 0L
  for (i in seq_len(nrow(sig))) { if (taken[i]) next
    taken[i] <- TRUE; n <- n + 1L
    taken[!taken & sig$CHR==sig$CHR[i] & abs(sig$POS-sig$POS[i]) < window_kb*1000] <- TRUE }
  n
}

placo_parallel <- function(Zmat, VarZ, CorZ, use_plus, cl, chunk = CHUNK) {
  n <- nrow(Zmat); starts <- seq(1, n, by = chunk)
  chunks <- lapply(starts, function(s) Zmat[s:min(s+chunk-1, n), , drop=FALSE])
  res <- parLapply(cl, chunks, function(Zc, V, C, plus) {
    vapply(seq_len(nrow(Zc)), function(i) {
      z <- Zc[i, ]
      if (plus) placo.plus(Z=z, VarZ=V, CorZ=C)$p.placo.plus else placo(Z=z, VarZ=V)$p.placo
    }, numeric(1))
  }, V=VarZ, C=CorZ, plus=use_plus)
  unlist(res, use.names = FALSE)
}

## --------------------------- 클러스터 ----------------------------------- ##
message("=== 클러스터 시작 (", N_CORE, " 워커) ===")
cl <- makeCluster(N_CORE)
clusterExport(cl, "PLACO_SRC", envir = environment())
invisible(clusterEvalQ(cl, { suppressMessages(source(PLACO_SRC)); NULL }))
if (!all(unlist(clusterEvalQ(cl, all(c("placo","placo.plus") %in% ls()))))) {
  stopCluster(cl); stop("워커 PLACO 로드 실패") }
message("  워커 PLACO 로드 확인")

## --------------------------- 페어별 실행 -------------------------------- ##
run_pair <- function(ex, ou) {
  sum_rds <- file.path(PLACO, paste0("SUMM_", ex, "_", ou, ".rds"))
  if (file.exists(sum_rds)) { message("  [건너뜀] ", ex, " × ", ou); return(readRDS(sum_rds)) }
  message("\n  ── ", ex, " × ", ou, " ──"); t0 <- proc.time()[3]
  cm <- load_pair(ex, ou)
  Zall <- as.matrix(cm[, .(Z1, Z2)]); Pall <- as.matrix(cm[, .(P1, P2)])
  VarZ <- var.placo(Zall, Pall, p.threshold = VAR_P)
  CorZ <- cor.pearson(Zall, Pall, p.threshold = VAR_P, returnMatrix = FALSE)
  use_plus <- abs(CorZ) > 0.02
  message(sprintf("    공통SNP %s | VarZ %.3f,%.3f | CorZ %+.4f | %s",
                  format(nrow(cm), big.mark=","), VarZ[1], VarZ[2], CorZ,
                  if (use_plus) "PLACO+" else "PLACO"))
  cand <- which(cm$P1 < SCREEN_P | cm$P2 < SCREEN_P)
  message("    스크리닝(P<", SCREEN_P, "): ", format(length(cand), big.mark=","), " SNP 계산 (", N_CORE, "코어)")
  pv <- placo_parallel(Zall[cand, , drop=FALSE], VarZ, CorZ, use_plus, cl)
  res <- data.table(SNP=cm$SNP[cand], CHR=cm$CHR[cand], POS=cm$POS[cand],
                    Z1=cm$Z1[cand], P1=cm$P1[cand], Z2=cm$Z2[cand], P2=cm$P2[cand],
                    T_PLACO=cm$Z1[cand]*cm$Z2[cand], P_PLACO=pv,
                    METHOD=if (use_plus) "PLACO+" else "PLACO")
  res[, FDR := p.adjust(P_PLACO, "BH")]
  sig <- res[P_PLACO < GW]; n_loci <- count_loci(sig)
  ## 참고 lambda(스크리닝된 집합이라 편향 — 최종검증은 전체SNP)
  pf <- res$P_PLACO[is.finite(res$P_PLACO) & res$P_PLACO > 0]
  lam <- median(qchisq(pf, 1, lower.tail=FALSE)) / qchisq(0.5, 1)
  mins <- (proc.time()[3]-t0)/60
  message(sprintf("    ★ 유의(5e-8) %d | loci %d | FDR<0.05 %d | NA %d  (%.1f분)",
                  nrow(sig), n_loci, sum(res$FDR<0.05, na.rm=TRUE), sum(is.na(pv)), mins))
  fwrite(res[order(P_PLACO)][1:min(10000, .N)], file.path(PLACO, paste0("PLACO_", ex, "_", ou, "_top.csv")))
  if (nrow(sig)) fwrite(sig[order(P_PLACO)], file.path(PLACO, paste0("SIG_", ex, "_", ou, ".csv")))
  saveRDS(res, file.path(PLACO, paste0("PLACO_", ex, "_", ou, ".rds")))
  s <- data.table(EXPOSURE=ex, OUTCOME=ou, METHOD=if (use_plus) "PLACO+" else "PLACO",
                  CORZ=round(CorZ,4), N_TESTED=nrow(res), N_NA=sum(is.na(pv)),
                  LAMBDA_scr=round(lam,3), N_SIG_GW=nrow(sig), N_LOCI=n_loci,
                  N_SIG_FDR=sum(res$FDR<0.05, na.rm=TRUE),
                  MIN_P=signif(min(res$P_PLACO, na.rm=TRUE),3), MIN=round(mins,1))
  saveRDS(s, sum_rds); rm(cm, res); gc(FALSE); s
}

## --------------------------- 사전검증(2만 SNP NA체크) ------------------- ##
message("\n=== 사전검증 (2만 SNP) ===")
vc <- load_pair(EXPOS[1], OUTC[1])
VarZv <- var.placo(as.matrix(vc[, .(Z1,Z2)]), as.matrix(vc[, .(P1,P2)]), p.threshold=VAR_P)
pv_test <- placo_parallel(as.matrix(vc[1:20000, .(Z1,Z2)]), VarZv, 0, FALSE, cl)
if (sum(is.na(pv_test)) > 0) { stopCluster(cl); stop("사전검증 실패: NA ", sum(is.na(pv_test)), "개") }
message("  통과 (NA 0). 본실행 시작."); rm(vc, pv_test); gc(FALSE)

## --------------------------- 실행 --------------------------------------- ##
summ <- rbindlist(lapply(EXPOS, function(ex)
          rbindlist(lapply(OUTC, function(ou) run_pair(ex, ou)))), fill = TRUE)
stopCluster(cl)

fwrite(summ, file.path(PLACO, "placo_summary.csv"))
message("\n================= Step 3 PLACO 완료 =================")
print(summ[, .(EXPOSURE, OUTCOME, METHOD, CORZ, N_TESTED, N_SIG_GW, N_LOCI, MIN_P)])
message("\n저장: ", PLACO, "  (PLACO_*.rds, SIG_*.csv, placo_summary.csv)")
message("점검: 천식×MDD/ANX 다면발현 loci가 LAVA 유의 locus와 겹치는지 / COA vs AOA")
message("다음(최종 rigor): SCREEN_P<-1(전체SNP) + 순열검증(exposure Z 셔플) → 인공물 배제")
writeLines(capture.output(sessionInfo()), file.path(PLACO, paste0("_step3_sessionInfo_", Sys.Date(), ".txt")))
