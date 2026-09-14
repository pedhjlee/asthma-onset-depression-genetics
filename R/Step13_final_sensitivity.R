###############################################################################
## Step 13 — 최종 민감도 분석 (재현 가능 버전)  2026-09-12
## Project: 01_GWAS_Asthma_MDD
##
## 목적: 2026-09-11 콘솔에서만 실행됐던 최종 계산을 스크립트로 고정.
##   (A) MR-PRESSO  — seed 20260911, NbDistribution 10,000 (주분석 8 + 재현 2)
##   (B) Steiger    — 이분형(liability) 척도, 유효표본수(Neff)로 통일 [A안]
##                    r 계산과 P 계산 모두 dat$samplesize(=Neff)를 기준으로 하고,
##                    사례:대조 비율만 원 GWAS 값을 유지.  9개 유병률 조합 × 10 분석.
##   (C) Table S8   — 주분석 8개 민감도 표 (이상치 보정 CI는 t 분포 기준)
##   (D) Table S11  — FinnGen 재현 표 (Steiger 열을 (B)로 갱신)
##   (E) Table S12  — Steiger 이분형 90행
##
## 실행: RStudio 콘솔에서
##   source("PATH/TO/PROJECT/R script/Step13_final_sensitivity.R")
## 실행 후: source(Step10_tables.R) 로 워크북 재생성.
##
## 입력: mr/MR_main_*.rds, mr/MR_reverse_*.rds, mr/REP_HOWARD_*.rds  (Step4·Step12 산출)
## 출력: 위 RDS에 presso/presso_n/presso_seed/steiger_binary 갱신 저장,
##       Result/tables_csv/TableS8_MR_sensitivity.csv, TableS11_FinnGen_replication.csv,
##       TableS12_steiger_binary.csv
###############################################################################

suppressMessages({ library(data.table); library(TwoSampleMR); library(MRPRESSO) })

## --------------------------- 경로·설정 ----------------------------------- ##
PROJ  <- "PATH/TO/PROJECT"
MRDIR <- file.path(PROJ, "mr")
TABO  <- file.path(PROJ, "Result", "tables_csv"); dir.create(TABO, showWarnings = FALSE, recursive = TRUE)

PRESSO_SEED <- 20260911
PRESSO_N    <- 10000
## TRUE 로 바꾸면 10개 모두 다시 돌림 (총 약 3시간). FALSE 면 RDS에 저장된
## seed 20260911·10,000회 결과가 있을 때 그대로 쓰고, 없는 것만 다시 돌림.
RERUN_PRESSO <- FALSE

## 분석 목록 (RDS 파일명, 노출, 결과)
AN <- data.table(
  FILE = c("MR_main_COA_MDD","MR_main_COA_ANX","MR_main_AOA_MDD","MR_main_AOA_ANX",
           "MR_reverse_MDD_COA","MR_reverse_MDD_AOA","MR_reverse_ANX_COA","MR_reverse_ANX_AOA",
           "REP_HOWARD_FGASTHMA","REP_HOWARD_FGCOA"),
  EXPO = c("COA","COA","AOA","AOA","MDD","MDD","ANX","ANX","MDD_HOWARD","MDD_HOWARD"),
  OUTC = c("MDD","ANX","MDD","ANX","COA","AOA","COA","AOA","FG_ASTHMA","FG_COA"),
  SET  = c(rep("main",8), rep("rep",2)))

## 사례/대조 수 (사례 비율 계산용) + 유병률 격자 (기준값은 2번째)
## 원 GWAS 보고값. Neff 는 dat$samplesize 에서 가져오므로 여기 N 은 비율에만 쓰임.
CC <- data.table(
  TRAIT  = c("COA","AOA","MDD","ANX","FG_ASTHMA","FG_COA","MDD_HOWARD"),
  N_CASE = c(13962, 26582, 357636, 122083, 61196, 8428, 170756),
  N_CTRL = c(300671, 300671, 1281936, 729602, 250433, 250433, 329443))
PREV_GRID <- list(
  COA = c(0.05, 0.10, 0.15), AOA = c(0.02, 0.05, 0.10),
  MDD = c(0.10, 0.15, 0.20), ANX = c(0.10, 0.16, 0.20),
  FG_ASTHMA = c(0.05, 0.10, 0.15), FG_COA = c(0.02, 0.05, 0.10),
  MDD_HOWARD = c(0.10, 0.15, 0.20))
BASE_PREV <- sapply(PREV_GRID, `[`, 2)

fmtP <- function(p) ifelse(is.na(p), NA_character_,
                    ifelse(p < 1e-300, "<1e-300", formatC(signif(p, 3), format = "g")))
## MR-PRESSO global P 는 NbDistribution 미만이면 "<1e-04" 문자열로 반환됨 → 그대로 유지
fmtPresso <- function(p) { if (is.null(p)) NA_character_ else if (is.numeric(p)) fmtP(p) else as.character(p) }
## 이상치가 없으면 MR-PRESSO 가 Distortion Test 를 반환하지 않음 → NULL 이면 NA (행 누락 방지)
nz <- function(x) if (is.null(x) || !length(x)) NA_real_ else x
ORc  <- function(b, se, q = 1.96) sprintf("%.2f (%.2f–%.2f)", exp(b), exp(b - q*se), exp(b + q*se))

## =========================== (A) MR-PRESSO =============================== ##
message("\n### (A) MR-PRESSO  seed ", PRESSO_SEED, ", NbDistribution ", PRESSO_N)
for (i in seq_len(nrow(AN))) {
  f <- file.path(MRDIR, paste0(AN$FILE[i], ".rds")); o <- readRDS(f)
  have <- !is.null(o$presso) && !inherits(o$presso, "try-error") &&
          isTRUE(as.numeric(o$presso_seed) == PRESSO_SEED) && isTRUE(as.numeric(o$presso_n) == PRESSO_N)
  if (have && !RERUN_PRESSO) { message("  ", AN$FILE[i], ": 저장된 최종 PRESSO 사용 (건너뜀)"); next }
  nIV <- nrow(o$dat); message("  ", AN$FILE[i], ": PRESSO 실행 (", nIV, " IV) ... ", format(Sys.time(), "%H:%M"))
  set.seed(PRESSO_SEED)
  pr <- try(mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                      SdOutcome = "se.outcome", SdExposure = "se.exposure",
                      OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = as.data.frame(o$dat),
                      NbDistribution = PRESSO_N, SignifThreshold = 0.05), silent = TRUE)
  if (inherits(pr, "try-error")) { message("    ★실패: ", pr); next }
  o$presso <- pr; o$presso_seed <- PRESSO_SEED; o$presso_n <- PRESSO_N
  saveRDS(o, f)
  message("    Global P = ", pr$`MR-PRESSO results`$`Global Test`$Pvalue, "  ", format(Sys.time(), "%H:%M"))
}

## =========================== (B) Steiger, Neff 통일 ====================== ##
## r_i = get_r_from_lor(beta_i, eaf_i, ncase_i, nctrl_i, prev)
##   ncase_i = samplesize_i × (N_CASE/(N_CASE+N_CTRL)),  nctrl_i = samplesize_i − ncase_i
##   → 총합 = samplesize_i (=Neff),  사례 비율 = 원 GWAS 비율
## P : mr_steiger(p_exp, p_out, n_exp = mean(samplesize.exposure), n_out = mean(samplesize.outcome), r_exp, r_out)
message("\n### (B) Steiger — liability 척도, 유효표본수 통일")
steiger_binary <- function(dat, expo, outc, prev_e, prev_o) {
  pe <- CC[TRAIT == expo]; po <- CC[TRAIT == outc]
  fe <- pe$N_CASE / (pe$N_CASE + pe$N_CTRL); fo <- po$N_CASE / (po$N_CASE + po$N_CTRL)
  ne <- dat$samplesize.exposure; no <- dat$samplesize.outcome
  r_e <- get_r_from_lor(dat$beta.exposure, dat$eaf.exposure, ne*fe, ne*(1-fe), prev_e)
  r_o <- get_r_from_lor(dat$beta.outcome,  dat$eaf.outcome,  no*fo, no*(1-fo), prev_o)
  R2e <- sum(r_e^2, na.rm = TRUE); R2o <- sum(r_o^2, na.rm = TRUE)
  st  <- mr_steiger(p_exp = dat$pval.exposure, p_out = dat$pval.outcome,
                    n_exp = mean(ne), n_out = mean(no), r_exp = sqrt(R2e), r_out = sqrt(R2o))
  list(R2_exp = R2e, R2_out = R2o, correct = st$correct_causal_direction, p = st$steiger_test,
       n_exp = mean(ne), n_out = mean(no))
}
S12 <- list()
for (i in seq_len(nrow(AN))) {
  f <- file.path(MRDIR, paste0(AN$FILE[i], ".rds")); o <- readRDS(f)
  ex <- AN$EXPO[i]; ou <- AN$OUTC[i]
  if (is.null(o$dat$samplesize.exposure) || is.null(o$dat$samplesize.outcome))
    stop(AN$FILE[i], ": dat 에 samplesize 열이 없음")
  grid <- CJ(pe = PREV_GRID[[ex]], po = PREV_GRID[[ou]])
  rows <- rbindlist(lapply(seq_len(nrow(grid)), function(k) {
    s <- steiger_binary(o$dat, ex, ou, grid$pe[k], grid$po[k])
    data.table(Analysis = AN$FILE[i], Exposure = ex, Outcome = ou,
               `Prev exposure` = grid$pe[k], `Prev outcome` = grid$po[k],
               `Neff exposure` = round(s$n_exp), `Neff outcome` = round(s$n_out),
               `R2 exposure` = signif(s$R2_exp, 3), `R2 outcome` = signif(s$R2_out, 3),
               `Correct direction` = s$correct, `Steiger P` = fmtP(s$p), P_raw = s$p,
               Scenario = ifelse(grid$pe[k] == BASE_PREV[ex] & grid$po[k] == BASE_PREV[ou], "base", "sensitivity"))
  }))
  S12[[i]] <- rows
  base <- rows[Scenario == "base"]
  o$steiger_binary <- list(table = rows, base = base, method = "Neff-unified (Step13, 2026-09-12)")
  saveRDS(o, f)
  message(sprintf("  %-22s base: R2e=%.4f R2o=%.5f correct=%s P=%s | 범위 P %s ~ %s",
    AN$FILE[i], base$`R2 exposure`, base$`R2 outcome`, base$`Correct direction`, base$`Steiger P`,
    fmtP(min(rows$P_raw)), fmtP(max(rows$P_raw))))
}
S12 <- rbindlist(S12)
fwrite(S12[, !"P_raw"], file.path(TABO, "TableS12_steiger_binary.csv"))
message("  → TableS12_steiger_binary.csv (", nrow(S12), " 행)  방향 올바름: ", sum(S12$`Correct direction`), "/", nrow(S12))

## =========================== (C) Table S8 (주분석 8) ===================== ##
## 이상치 보정 CI: MR-PRESSO 는 절편 없는 가중 lm 이므로 df = (남은 IV 수) − 1, t 분위수 사용.
message("\n### (C) Table S8")
loo_range <- function(loo) { l <- loo[loo$SNP != "All", ]; sprintf("%.2f–%.2f", exp(min(l$b)), exp(max(l$b))) }
s8 <- rbindlist(lapply(which(AN$SET == "main"), function(i) {
  o <- readRDS(file.path(MRDIR, paste0(AN$FILE[i], ".rds")))
  r <- as.data.table(o$res); g <- function(mm, w) r[method == mm][[w]]
  ivw <- "Inverse variance weighted"; nIV <- nrow(o$dat)
  pr <- o$presso; main <- pr$`Main MR results`; corr <- main[2, ]
  out_idx <- pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  out_snps <- if (is.numeric(out_idx)) o$dat$SNP[out_idx] else character(0)
  n_out <- length(out_snps); df <- nIV - n_out - 1
  t_q <- qt(0.975, df)
  p_t <- 2 * pt(-abs(corr$`Causal Estimate` / corr$Sd), df)
  if (is.finite(p_t) && !is.na(corr$`P-value`) && abs(p_t - corr$`P-value`) > 1e-6)
    message(sprintf("    ⚠ %s: t 재계산 P=%.4g vs PRESSO 저장 P=%.4g (df=%d) 불일치 — df 가정 점검",
                    AN$FILE[i], p_t, corr$`P-value`, df))
  het <- o$het[o$het$method == ivw, ]
  data.table(
    Direction = ifelse(AN$SET[i] == "main" & AN$EXPO[i] %in% c("COA","AOA"), "Asthma → psychiatric", "Psychiatric → asthma"),
    Exposure = AN$EXPO[i], Outcome = AN$OUTC[i], `N IV` = nIV,
    `IVW OR (95% CI)` = ORc(g(ivw,"b"), g(ivw,"se")), `IVW P` = fmtP(g(ivw,"pval")),
    `Weighted median OR (95% CI)` = ORc(g("Weighted median","b"), g("Weighted median","se")),
    `Weighted median P` = fmtP(g("Weighted median","pval")),
    `Weighted mode OR (95% CI)` = ORc(g("Weighted mode","b"), g("Weighted mode","se")),
    `MR-Egger OR (95% CI)` = ORc(g("MR Egger","b"), g("MR Egger","se")),
    `Egger intercept P` = fmtP(o$plei$pval),
    `Cochran Q (df), P` = sprintf("%.1f (%d), %s", het$Q, het$Q_df, formatC(het$Q_pval, format = "e", digits = 2)),
    `Leave-one-out OR range` = loo_range(o$loo),
    `MR-PRESSO global P` = fmtPresso(pr$`MR-PRESSO results`$`Global Test`$Pvalue),
    `MR-PRESSO outliers` = if (n_out) paste0(n_out, ": ", paste(sort(out_snps), collapse = ", ")) else "0",
    `Outlier-corrected OR (95% CI)` = if (!is.na(corr$`Causal Estimate`)) ORc(corr$`Causal Estimate`, corr$Sd, t_q) else "—",
    `Outlier-corrected P` = fmtP(corr$`P-value`),
    `Distortion P` = fmtP(nz(pr$`MR-PRESSO results`$`Distortion Test`$Pvalue)),
    `Steiger P` = o$steiger_binary$base$`Steiger P`)
}))
fwrite(s8, file.path(TABO, "TableS8_MR_sensitivity.csv"))
print(s8[, .(Exposure, Outcome, `N IV`, `IVW OR (95% CI)`, `IVW P`, `Outlier-corrected OR (95% CI)`, `Outlier-corrected P`, `Steiger P`)])

## =========================== (D) Table S11 (FinnGen) ===================== ##
message("\n### (D) Table S11")
NM <- c(REP_HOWARD_FGASTHMA = "FinnGen asthma (all)", REP_HOWARD_FGCOA = "FinnGen childhood-onset asthma (age<16)")
s11 <- rbindlist(lapply(which(AN$SET == "rep"), function(i) {
  o <- readRDS(file.path(MRDIR, paste0(AN$FILE[i], ".rds")))
  r <- as.data.table(o$res); g <- function(mm, w) r[method == mm][[w]]; ivw <- "Inverse variance weighted"
  pr <- o$presso; corr <- pr$`Main MR results`[2, ]
  out_idx <- pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  n_out <- if (is.numeric(out_idx)) length(out_idx) else 0L; df <- nrow(o$dat) - n_out - 1
  data.table(Outcome = NM[AN$FILE[i]], `N IV` = nrow(o$dat),
    `F (mean)` = round(mean(o$dat$beta.exposure^2 / o$dat$se.exposure^2), 1),
    `IVW OR (95% CI)` = ORc(g(ivw,"b"), g(ivw,"se")), `IVW P` = fmtP(g(ivw,"pval")),
    `Weighted median` = ORc(g("Weighted median","b"), g("Weighted median","se")),
    `Weighted mode` = ORc(g("Weighted mode","b"), g("Weighted mode","se")),
    `MR-Egger` = ORc(g("MR Egger","b"), g("MR Egger","se")),
    `Egger intercept P` = fmtP(o$plei$pval),
    `Cochran Q P` = fmtP(o$het$Q_pval[o$het$method == ivw]),
    `MR-PRESSO global P` = fmtPresso(pr$`MR-PRESSO results`$`Global Test`$Pvalue),
    `MR-PRESSO outliers` = n_out,
    `Outlier-corrected OR` = if (!is.na(corr$`Causal Estimate`)) ORc(corr$`Causal Estimate`, corr$Sd, qt(0.975, df)) else "no outlier",
    `Outlier-corrected P` = if (!is.na(corr$`Causal Estimate`)) fmtP(corr$`P-value`) else NA_character_,
    `Distortion P` = fmtP(nz(pr$`MR-PRESSO results`$`Distortion Test`$Pvalue)),
    `Steiger P` = o$steiger_binary$base$`Steiger P`)
}))
stopifnot(nrow(s11) == 2)   # 소아천식 행 누락 시 즉시 중단
## 전체 vs 소아 차이 (Step12 §3 과 동일: 오차상관 0.5 가정, 탐색적)
a  <- as.data.table(readRDS(file.path(MRDIR, "REP_HOWARD_FGASTHMA.rds"))$dat)
b_ <- as.data.table(readRDS(file.path(MRDIR, "REP_HOWARD_FGCOA.rds"))$dat)
m  <- merge(a[, .(SNP, bx = beta.exposure, b1 = beta.outcome, s1 = se.outcome)],
            b_[, .(SNP, b2 = beta.outcome, s2 = se.outcome)], by = "SNP")
rho <- 0.5; m[, `:=`(d = b1 - b2, sd = sqrt(s1^2 + s2^2 - 2*rho*s1*s2))]
fit <- lm(d ~ bx - 1, weights = 1/sd^2, data = m)
bb <- coef(summary(fit))[1, 1]; se <- coef(summary(fit))[1, 2] / min(1, summary(fit)$sigma)
s11 <- rbind(s11, data.table(Outcome = "Difference (all vs childhood)", `N IV` = nrow(m),
  `IVW OR (95% CI)` = ORc(bb, se), `IVW P` = fmtP(2*pnorm(-abs(bb/se)))), fill = TRUE)
fwrite(s11, file.path(TABO, "TableS11_FinnGen_replication.csv"))
print(s11[, .(Outcome, `N IV`, `IVW OR (95% CI)`, `IVW P`, `Outlier-corrected OR`, `Outlier-corrected P`, `Distortion P`, `Steiger P`)])

## =========================== 종합 ======================================= ##
message("\n### 종합")
message("  S8/S11/S12 저장: ", TABO)
message("  다음: source(\"", file.path(PROJ, "R script", "Step10_tables.R"), "\")  → 워크북 재생성")
writeLines(capture.output(sessionInfo()), file.path(MRDIR, paste0("_step13_sessionInfo_", Sys.Date(), ".txt")))
