###############################################################################
## Step 11b — MVMR 민감도: 약한 도구 보정 (Q-minimization, MVMR::qhet_mvmr)
## Project: 01_GWAS_Asthma_MDD
##
## 배경: Step11 주분석에서 SmkInit 조건부 F=3.5 (WithoutUKB, 9 SNP). 약한 도구는
##       IVW-MVMR 추정치를 편향시킬 수 있음 → Q-minimization 추정(Sanderson 2021)으로
##       MDD 직접효과가 유지되는지 확인. 노출 간 표현형 상관(pcor)은 개인자료 없이
##       추정 불가 → 문헌값 + 0/보수값으로 시나리오 민감도.
## 실행 전략 (Windows = 단일코어): 점추정은 전 시나리오, 부트스트랩 CI는 주분석
##       (AOA × pcor='lit') 한 개만 N_BOOT회.
## 입력: mr/MVMR_main_AOA.rds, MVMR_main_COA.rds  (Step11 산출)
## 출력: mr/mvmr_qhet_summary.csv
###############################################################################
suppressMessages({ library(data.table); library(MVMR) })
PROJ  <- "PATH/TO/PROJECT"
MRDIR <- file.path(PROJ, "mr"); setwd(MRDIR)
SEED  <- 20260714; set.seed(SEED)
N_BOOT  <- 200                        # BCa 부트스트랩 (주분석 1개만)
CI_FOR  <- list(OUTCOME = "AOA", PCOR = "lit")

## 표현형 상관 시나리오. 문헌 근사: MDD–BMI ≈ 0.10, MDD–SmkInit ≈ 0.20, BMI–SmkInit ≈ 0.05
PCOR <- list(
  none = c(MDD_BMI = 0,    MDD_SMK = 0,    BMI_SMK = 0),
  lit  = c(MDD_BMI = 0.10, MDD_SMK = 0.20, BMI_SMK = 0.05),
  high = c(MDD_BMI = 0.20, MDD_SMK = 0.30, BMI_SMK = 0.10))
make_pcor <- function(v, ord) {
  m <- diag(3); dimnames(m) <- list(ord, ord)
  m["MDD","BMI"] <- m["BMI","MDD"] <- v["MDD_BMI"]
  m["MDD","SMK"] <- m["SMK","MDD"] <- v["MDD_SMK"]
  m["BMI","SMK"] <- m["SMK","BMI"] <- v["BMI_SMK"]; m
}
## qhet CI는 paste(lcb, ucb, sep="-") 문자열 → 구분자 '-'와 음수 부호가 섞이므로 순서대로 파싱
parse_ci <- function(s) {
  s <- trimws(s); sgn1 <- 1
  if (startsWith(s, "-")) { sgn1 <- -1; s <- substring(s, 2) }
  p <- strsplit(s, "-", fixed = TRUE)[[1]]          # "a-b" → c(a,b); "a--b" → c(a,"",b)
  lo <- sgn1 * as.numeric(p[1])
  hi <- if (length(p) >= 3 && p[2] == "") -as.numeric(p[3]) else as.numeric(p[2])
  c(lo, hi)
}

run_qhet <- function(outc) {
  x   <- readRDS(paste0("MVMR_main_", outc, ".rds"))$mvdat
  ord <- colnames(x$exposure_beta)
  rin <- format_mvmr(BXGs = x$exposure_beta, BY = x$outcome_beta, seBXGs = x$exposure_se,
                     seBY = x$outcome_se, RSID = rownames(x$exposure_beta))
  rbindlist(lapply(names(PCOR), function(sc) {
    message(sprintf("\n  ── %s | pcor=%s ──", outc, sc))
    pc  <- make_pcor(PCOR[[sc]], ord)
    gc  <- phenocov_mvmr(pcor = pc, seBXGs = x$exposure_se)
    cF  <- suppressWarnings(strength_mvmr(rin, gencov = gc))
    Q   <- suppressWarnings(pleiotropy_mvmr(rin, gencov = gc))
    pt  <- suppressWarnings(qhet_mvmr(rin, pcor = pc, CI = FALSE, iterations = 1, ncores = 1))
    est <- as.numeric(pt[, 1]); lo <- hi <- rep(NA_real_, length(est))
    message("    Q-min OR: ", paste(sprintf("%s=%.3f", ord, exp(est)), collapse = ", "),
            " | 조건부 F: ", paste(sprintf("%s=%.1f", ord, as.numeric(cF)), collapse = ", "))
    if (outc == CI_FOR$OUTCOME && sc == CI_FOR$PCOR) {
      message("    부트스트랩 CI ", N_BOOT, "회 시작 (", format(Sys.time(), "%H:%M"), ") — 단일코어, 수 분~수십 분")
      qh <- suppressWarnings(qhet_mvmr(rin, pcor = pc, CI = TRUE, iterations = N_BOOT, ncores = 1))
      ci <- lapply(as.character(qh[, 2]), parse_ci); lo <- sapply(ci, `[`, 1); hi <- sapply(ci, `[`, 2)
      message("    CI 완료 (", format(Sys.time(), "%H:%M"), ")")
    }
    data.table(OUTCOME = outc, PCOR = sc, EXPOSURE = ord, METHOD = "Q-minimization",
               N_SNP = nrow(x$exposure_beta), N_BOOT = ifelse(is.na(lo), NA_integer_, N_BOOT),
               B = round(est, 4), OR = round(exp(est), 3), LCI = round(exp(lo), 3), UCI = round(exp(hi), 3),
               COND_F = round(as.numeric(cF), 1), Q = round(Q$Qstat, 1), Q_P = signif(Q$Qpval, 3))
  }))
}

res <- rbind(run_qhet("AOA"), run_qhet("COA"))
fwrite(res, "mvmr_qhet_summary.csv")

message("\n### 판정 (MDD → AOA, Q-minimization)")
print(res[OUTCOME == "AOA" & EXPOSURE == "MDD", .(PCOR, OR, LCI, UCI, COND_F)])
message("### COA 대조 (MDD)")
print(res[OUTCOME == "COA" & EXPOSURE == "MDD", .(PCOR, OR, COND_F)])
message("저장: mvmr_qhet_summary.csv")
writeLines(capture.output(sessionInfo()), file.path(MRDIR, paste0("_step11b_sessionInfo_", Sys.Date(), ".txt")))
