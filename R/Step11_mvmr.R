###############################################################################
## Step 11 — Multivariable MR (MVMR): 우울(MDD) → 천식, BMI·흡연개시 교란 보정
## Project: 01_GWAS_Asthma_MDD
##
## 목적: Step4 역방향 MR에서 MDD→AOA IVW OR 1.19 (PRESSO 1.23). 우울의 효과가
##       BMI·흡연(우울과 천식 공통 위험요인)을 경유/교란한 것이 아닌지 검증.
##   주분석 : MDD + BMI + SmkInit → AOA
##   보조   : 동일 모형 → COA (음성 대조; univariable에서 null)
##   부분모형: MDD+BMI → AOA, MDD+SmkInit → AOA (어느 공변량이 효과를 바꾸는지)
## 판정: [MDD] 직접효과 OR>1 & P<0.05 유지 → BMI·흡연 교란 배제.
##
## 방법: IVW-MVMR(주, TwoSampleMR::mv_multiple) + MVMR-Egger + MVMR-median
##       (MendelianRandomization), 조건부 F(MVMR::strength_mvmr), Q(pleiotropy_mvmr).
## 표본중복: MDD(noUKBB)·BMI(GIANT2015, UKB없음)·SmkInit(WithoutUKB) vs 천식(Ferreira2019,
##       UKB포함) → 노출-결과 간 중복 없음.
##
## 데이터: raw/SNP_gwas_mc_merge_nogc.tbl.uniq.gz  (BMI, Locke 2015 GIANT; CHR/POS 없음→MDD/SMK로 보충)
##         raw/SmokingInitiation.WithoutUKB.txt.gz  (GSCAN Liu 2019; ALT=effect allele)
## 출력: harmonized/BMI.rds, SMK.rds | mr/mvmr_summary.csv, mr/MVMR_*.rds, mr/fig_mvmr_forest.png
###############################################################################

## --------------------------- 패키지 (없으면 자동 설치) --------------------- ##
.repo <- "https://cloud.r-project.org"
for (p in c("remotes","data.table","R.utils","ggplot2","dplyr","MendelianRandomization"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = .repo)
## MRCIEU 패키지: r-universe 바이너리(Rtools 불필요) → 실패 시 GitHub
.mrc <- c("TwoSampleMR", "ieugwasr", "genetics.binaRies")
.need <- .mrc[!sapply(.mrc, requireNamespace, quietly = TRUE)]
if (length(.need)) install.packages(.need, repos = c("https://mrcieu.r-universe.dev", .repo))
.gh <- c(TwoSampleMR = "MRCIEU/TwoSampleMR", ieugwasr = "MRCIEU/ieugwasr",
         genetics.binaRies = "MRCIEU/genetics.binaRies", MVMR = "WSpiller/MVMR")
for (p in names(.gh)) if (!requireNamespace(p, quietly = TRUE)) remotes::install_github(.gh[[p]], upgrade = "never")
suppressMessages({ library(data.table); library(TwoSampleMR); library(ggplot2)
                   library(MendelianRandomization); library(MVMR) })

## --------------------------- 경로 --------------------------------------- ##
PROJ  <- "PATH/TO/PROJECT"
RAW   <- file.path(PROJ, "raw")
HARM  <- file.path(PROJ, "harmonized")
MRDIR <- file.path(PROJ, "mr")
dir.create(MRDIR, showWarnings = FALSE, recursive = TRUE)
setwd(MRDIR)

## --------------------------- 설정 --------------------------------------- ##
IV_PVAL  <- 5e-8; CLUMP_R2 <- 0.001; CLUMP_KB <- 10000
MAF_MIN  <- 0.01; SEED <- 20260714; set.seed(SEED)
REHARMONIZE <- FALSE   # TRUE면 BMI/SMK.rds 있어도 다시 생성

## 로컬 clumping (Step4와 동일: ASCII 경로 EUR 패널; 없으면 ld_ref/에서 복사)
LOCAL_BFILE <- "PATH/TO/ldref_ascii/EUR"
{ ascii_dir <- "PATH/TO/ldref_ascii"; dir.create(ascii_dir, showWarnings = FALSE, recursive = TRUE)
  for (ext in c("bed","bim","fam")) {
    d <- file.path(ascii_dir, paste0("EUR.", ext)); s <- file.path(PROJ, "ld_ref", paste0("EUR.", ext))
    if (!file.exists(d) && file.exists(s)) { message("  LD패널 ASCII 복사: EUR.", ext); file.copy(s, d) }
  } }
LOCAL_PLINK <- if (requireNamespace("genetics.binaRies", quietly = TRUE)) genetics.binaRies::get_plink_binary() else NULL
USE_LOCAL   <- !is.null(LOCAL_PLINK) && file.exists(paste0(LOCAL_BFILE, ".bed"))
message("  로컬 clumping: ", if (USE_LOCAL) "준비됨" else "★참조패널 없음 → API 폴백")

## =========================================================================== ##
## PART 0 — BMI·흡연개시 harmonize (12컬럼 스키마: SNP CHR POS EA NEA EAF BETA SE P N TRAIT TYPE)
## =========================================================================== ##
finalize_lite <- function(dt, trait, type) {
  n0 <- nrow(dt)
  dt[, `:=`(EA = toupper(EA), NEA = toupper(NEA))]
  dt <- dt[EA %in% c("A","C","G","T") & NEA %in% c("A","C","G","T") & EA != NEA]
  dt[, CHR := suppressWarnings(as.integer(CHR))][, POS := suppressWarnings(as.integer(POS))]
  dt <- dt[CHR %in% 1:22 & !is.na(POS)]
  dt <- dt[is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) & P > 0 & P <= 1 & is.finite(N) & N > 0]
  dt <- dt[is.na(EAF) | (EAF >= MAF_MIN & EAF <= 1 - MAF_MIN)]
  dt <- dt[grepl("^rs", SNP)]
  dup <- dt$SNP[duplicated(dt$SNP)]; if (length(dup)) dt <- dt[!SNP %in% dup]
  dt[, `:=`(TRAIT = trait, TYPE = type)]
  dt <- dt[, .(SNP = as.character(SNP), CHR, POS, EA, NEA, EAF = as.numeric(EAF),
               BETA = as.numeric(BETA), SE = as.numeric(SE), P = as.numeric(P), N = as.numeric(N), TRAIT, TYPE)]
  message(sprintf("  [%s] %d -> %d SNP | EAF NA %d | top %s P=%.2e", trait, n0, nrow(dt),
                  sum(is.na(dt$EAF)), dt$SNP[which.min(dt$P)], min(dt$P)))
  dt
}
## EAF 결측 보충: 참조(MDD) 대립유전자 정렬 후 채움
fill_eaf <- function(dt, ref) {
  r <- ref[, .(SNP, rEA = EA, rNEA = NEA, rEAF = EAF)]
  dt <- merge(dt, r, by = "SNP", all.x = TRUE, sort = FALSE)
  dt[is.na(EAF) & EA == rEA & NEA == rNEA, EAF := rEAF]
  dt[is.na(EAF) & EA == rNEA & NEA == rEA, EAF := 1 - rEAF]
  dt[, c("rEA","rNEA","rEAF") := NULL]; dt
}
pick <- function(nm, cands) { h <- cands[cands %in% nm]; if (!length(h)) stop("컬럼 못 찾음: ", paste(cands, collapse="/")); h[1] }

MDD <- as.data.table(readRDS(file.path(HARM, "MDD.rds")))

## ---- 흡연개시 (GSCAN). 파일명 WithoutUKB 우선, 없으면 전체판 ----------------
harm_SMK <- function() {
  f <- c("SmokingInitiation.WithoutUKB.txt.gz", "SmokingInitiation.txt.gz")
  f <- f[file.exists(file.path(RAW, f))][1]; if (is.na(f)) stop("흡연개시 파일 없음(raw/)")
  message("  흡연개시 파일: ", f)
  d <- fread(file.path(RAW, f), showProgress = FALSE)
  nm <- names(d)
  c_chr <- pick(nm, c("CHROM","CHR","#CHROM")); c_pos <- pick(nm, c("POS","BP"))
  c_rs  <- pick(nm, c("RSID","SNP","rsID"));   c_ref <- pick(nm, c("REF")); c_alt <- pick(nm, c("ALT"))
  c_af  <- pick(nm, c("AF","EAF"));            c_p   <- pick(nm, c("PVALUE","P","PVAL"))
  c_b   <- pick(nm, c("BETA","beta"));         c_se  <- pick(nm, c("SE","se"))
  c_n   <- if ("EFFECTIVE_N" %in% nm) "EFFECTIVE_N" else pick(nm, c("N"))   # binary → effective N
  x <- d[, .(SNP = get(c_rs), CHR = get(c_chr), POS = get(c_pos),
             EA = get(c_alt), NEA = get(c_ref),                # GSCAN: ALT = effect allele
             EAF = as.numeric(get(c_af)), BETA = get(c_b), SE = get(c_se), P = get(c_p), N = get(c_n))]
  rm(d); gc()
  x <- finalize_lite(x, "SMK", "binary")
  x <- fill_eaf(x, MDD)
  message("  SMK EAF 보충 후 NA: ", sum(is.na(x$EAF))); x
}
## ---- BMI (GIANT Locke 2015; SNP A1 A2 Freq1.Hapmap b se p N; A1 = effect) ---
harm_BMI <- function(posmap) {
  d <- fread(file.path(RAW, "SNP_gwas_mc_merge_nogc.tbl.uniq.gz"), showProgress = FALSE)
  x <- d[, .(SNP, EA = A1, NEA = A2, EAF = as.numeric(Freq1.Hapmap), BETA = b, SE = se, P = p, N = N)]
  x <- merge(x, posmap, by = "SNP", all.x = FALSE, sort = FALSE)     # CHR/POS 보충 (없는 SNP 제외)
  message("  BMI 좌표 매칭: ", nrow(x), "/", nrow(d))
  x <- finalize_lite(x, "BMI", "continuous")
  fill_eaf(x, MDD)
}

f_smk <- file.path(HARM, "SMK.rds"); f_bmi <- file.path(HARM, "BMI.rds")
if (REHARMONIZE || !file.exists(f_smk)) { message("== SMK harmonize =="); SMK <- harm_SMK(); saveRDS(SMK, f_smk) } else SMK <- as.data.table(readRDS(f_smk))
if (REHARMONIZE || !file.exists(f_bmi)) {
  message("== BMI harmonize ==")
  posmap <- unique(rbind(MDD[, .(SNP, CHR, POS)], SMK[, .(SNP, CHR, POS)]), by = "SNP")
  BMI <- harm_BMI(posmap); saveRDS(BMI, f_bmi); rm(posmap)
} else BMI <- as.data.table(readRDS(f_bmi))

## =========================================================================== ##
## PART 1 — MVMR 함수
## =========================================================================== ##
EXPO <- list(MDD = MDD, BMI = BMI, SMK = SMK)
do_clump <- function(snp, pval, id = "x") {
  tb <- dplyr::tibble(rsid = snp, pval = pval, id = id)
  if (USE_LOCAL) {
    out <- try(ieugwasr::ld_clump(tb, clump_kb = CLUMP_KB, clump_r2 = CLUMP_R2, bfile = LOCAL_BFILE, plink_bin = LOCAL_PLINK), silent = TRUE)
    if (!inherits(out, "try-error")) return(out$rsid)
    message("    로컬 clumping 실패 → API")
  }
  out <- try(ieugwasr::ld_clump(tb, clump_kb = CLUMP_KB, clump_r2 = CLUMP_R2, pop = "EUR"), silent = TRUE)
  if (inherits(out, "try-error")) stop("clumping 실패"); out$rsid
}
fmt_exp <- function(dt, tr) {
  e <- format_data(as.data.frame(dt), type = "exposure", snp_col="SNP", beta_col="BETA", se_col="SE",
                   effect_allele_col="EA", other_allele_col="NEA", eaf_col="EAF", pval_col="P",
                   samplesize_col="N", chr_col="CHR", pos_col="POS")
  e$exposure <- tr; e$id.exposure <- tr; e
}
fmt_out <- function(dt, tr) {
  o <- format_data(as.data.frame(dt), type = "outcome", snp_col="SNP", beta_col="BETA", se_col="SE",
                   effect_allele_col="EA", other_allele_col="NEA", eaf_col="EAF", pval_col="P", samplesize_col="N")
  o$outcome <- tr; o$id.outcome <- tr; o
}

## 노출별 유의 SNP → clump → 합집합 → 합집합 재clump(min P 기준) → 전 노출·결과에 존재하는 SNP만
build_instruments <- function(exps) {
  common <- Reduce(intersect, lapply(EXPO[exps], function(d) d$SNP))
  ivs <- lapply(exps, function(tr) {
    d <- EXPO[[tr]][P < IV_PVAL & SNP %in% common]
    cl <- do_clump(d$SNP, d$P, tr); message(sprintf("    %s: p<5e-8 %d → clump %d", tr, nrow(d), length(cl))); cl
  })
  u <- unique(unlist(ivs))
  minp <- rbindlist(lapply(exps, function(tr) EXPO[[tr]][SNP %in% u, .(SNP, P)]))[, .(P = min(P)), by = SNP]
  keep <- do_clump(minp$SNP, minp$P, "union")
  message(sprintf("    합집합 %d → 재clump %d SNP", length(u), length(keep)))
  list(snps = keep, per_exposure = setNames(ivs, exps))
}

run_mvmr <- function(exps, outc, tag) {
  message(sprintf("\n  ── MVMR [%s]: %s → %s ──", tag, paste(exps, collapse="+"), outc))
  iv <- build_instruments(exps)
  exposure_dat <- do.call(rbind, lapply(exps, function(tr) fmt_exp(EXPO[[tr]][SNP %in% iv$snps], tr)))
  OUT <- as.data.table(readRDS(file.path(HARM, paste0(outc, ".rds"))))
  outcome_dat  <- fmt_out(OUT[SNP %in% iv$snps], outc)
  mvdat <- mv_harmonise_data(exposure_dat, outcome_dat, harmonise_strictness = 2)
  nsnp  <- nrow(mvdat$exposure_beta)
  message("    최종 도구변수: ", nsnp)

  ## (1) IVW-MVMR (TwoSampleMR)
  ivw <- mv_multiple(mvdat, intercept = FALSE, instrument_specific = FALSE, pval_threshold = IV_PVAL, plots = FALSE)$result
  ## (2) MVMR-Egger / median (MendelianRandomization)
  ord  <- colnames(mvdat$exposure_beta)
  mvi  <- mr_mvinput(bx = mvdat$exposure_beta, bxse = mvdat$exposure_se,
                     by = mvdat$outcome_beta, byse = mvdat$outcome_se, exposure = ord, outcome = outc)
  egg  <- try(mr_mvegger(mvi, orientate = 1), silent = TRUE)
  med  <- try(mr_mvmedian(mvi, iterations = 1000, seed = SEED), silent = TRUE)
  ## (3) 조건부 F, Q (MVMR 패키지; gencov=0 = 표본중복 없음 가정)
  rin  <- format_mvmr(BXGs = mvdat$exposure_beta, BY = mvdat$outcome_beta, seBXGs = mvdat$exposure_se,
                      seBY = mvdat$outcome_se, RSID = rownames(mvdat$exposure_beta))
  cF   <- try(strength_mvmr(rin, gencov = 0), silent = TRUE)
  Q    <- try(pleiotropy_mvmr(rin, gencov = 0), silent = TRUE)
  if (!inherits(cF, "try-error")) message("    조건부 F: ", paste(sprintf("%s=%.1f", ord, as.numeric(cF)), collapse=", "))
  if (!inherits(Q,  "try-error")) message(sprintf("    Q=%.1f (p=%.3g)", Q$Qstat, Q$Qpval))
  if (!inherits(egg,"try-error")) message(sprintf("    MVMR-Egger 절편 p=%.3g", egg@Pvalue.Int))

  row <- function(method, ex, b, se, p, extra = list()) data.table(
    TAG = tag, MODEL = paste(exps, collapse="+"), OUTCOME = outc, EXPOSURE = ex, METHOD = method, N_SNP = nsnp,
    B = round(b, 4), SE = round(se, 4), OR = round(exp(b), 3), LCI = round(exp(b - 1.96*se), 3),
    UCI = round(exp(b + 1.96*se), 3), P = signif(p, 3))
  res <- rbindlist(list(
    rbindlist(lapply(seq_len(nrow(ivw)), function(i) row("IVW-MVMR", ivw$exposure[i], ivw$b[i], ivw$se[i], ivw$pval[i]))),
    if (!inherits(egg,"try-error")) rbindlist(lapply(seq_along(ord), function(i) row("MVMR-Egger", ord[i], egg@Estimate[i], egg@StdError.Est[i], egg@Pvalue.Est[i]))),
    if (!inherits(med,"try-error")) rbindlist(lapply(seq_along(ord), function(i) row("MVMR-median", ord[i], med@Estimate[i], med@StdError[i], med@Pvalue[i])))
  ), fill = TRUE)
  if (!inherits(cF, "try-error")) res[, COND_F := round(as.numeric(cF)[match(EXPOSURE, ord)], 1)]
  if (!inherits(Q,  "try-error")) res[, `:=`(Q = round(Q$Qstat, 1), Q_P = signif(Q$Qpval, 3))]
  if (!inherits(egg,"try-error")) res[, EGGER_INT_P := signif(egg@Pvalue.Int, 3)]
  print(res[, .(EXPOSURE, METHOD, N_SNP, OR, LCI, UCI, P)])
  saveRDS(list(mvdat = mvdat, ivw = ivw, egger = egg, median = med, condF = cF, Q = Q, iv = iv),
          paste0("MVMR_", tag, "_", outc, ".rds"))
  res
}

## =========================================================================== ##
## PART 2 — 실행
## =========================================================================== ##
message("\n### 주분석: MDD + BMI + SMK → AOA")
r_main <- run_mvmr(c("MDD","BMI","SMK"), "AOA", "main")
message("\n### 보조: MDD + BMI + SMK → COA (음성 대조)")
r_coa  <- run_mvmr(c("MDD","BMI","SMK"), "COA", "main")
message("\n### 부분모형 → AOA")
r_p1   <- run_mvmr(c("MDD","BMI"), "AOA", "partial_BMI")
r_p2   <- run_mvmr(c("MDD","SMK"), "AOA", "partial_SMK")

summ <- rbindlist(list(r_main, r_coa, r_p1, r_p2), fill = TRUE)
fwrite(summ, "mvmr_summary.csv")

## --------------------------- 그림: forest (IVW-MVMR) ---------------------- ##
fp <- summ[METHOD == "IVW-MVMR"]
fp[, label := paste0(EXPOSURE, "  [", MODEL, " → ", OUTCOME, "]")][, label := factor(label, levels = rev(unique(label)))]
p <- ggplot(fp, aes(OR, label, colour = EXPOSURE == "MDD")) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.2) + geom_point(size = 2.5) +
  scale_x_log10() + scale_colour_manual(values = c("grey40", "#b2182b"), guide = "none") +
  labs(x = "OR per SD/log-odds (95% CI), IVW-MVMR", y = NULL, title = "Multivariable MR: direct effects on asthma") +
  theme_minimal(base_size = 11)
ggsave("fig_mvmr_forest.png", p, width = 8.5, height = 0.45 * nrow(fp) + 1.5, dpi = 150)

## --------------------------- 판정 ---------------------------------------- ##
m <- r_main[EXPOSURE == "MDD" & METHOD == "IVW-MVMR"]
message("\n### 판정")
message(sprintf("  MDD → AOA 직접효과 (BMI·흡연 보정): OR %.3f (%.3f–%.3f), P=%.3g  [Step4 univariable IVW OR 1.19]",
                m$OR, m$LCI, m$UCI, m$P))
message("  → ", if (m$OR > 1 && m$P < 0.05) "유지: BMI·흡연 교란 배제 (논문1 방어 완료)" else "★약화/소실: 교란 가능성 → Discussion에서 다룰 것")
message("\n저장: ", MRDIR, " (mvmr_summary.csv, MVMR_*.rds, fig_mvmr_forest.png)")
writeLines(capture.output(sessionInfo()), file.path(MRDIR, paste0("_step11_sessionInfo_", Sys.Date(), ".txt")))
