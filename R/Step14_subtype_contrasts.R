## Step 14 — Subtype contrasts in one table (Table E9, two panels)
##   Panel A: genetic-correlation contrasts  rg(AOA,X) − rg(COA,X)  and differences of these differences
##            (delta method with the full LDSC sampling covariance matrix from Step1's ldsc_result.rds)
##   Panel B: MR contrasts
##            (1) depression → AOA vs COA (UK Biobank): IVW regression of per-variant differences,
##                correlated sampling errors r = 0.10 (standardized LDSC cross-trait intercept), multiplicative random effects
##            (2) depression → FinnGen asthma (all) vs childhood asthma: same, assumed error correlation 0.5 (exploratory; as Step12/13)
##            (3) depression vs anxiety → AOA: z-test assuming independent estimates (exploratory)
## 실행: source("R/Step14_subtype_contrasts.R")
## 출력: Result/tables_csv/TableS9_subtype_contrasts.csv  (Step10 이 그대로 읽어 Table E9 로 넣음)
suppressPackageStartupMessages(library(data.table))
PROJ  <- "PATH/TO/PROJECT"
LDSC  <- file.path(PROJ, "ldsc"); MRDIR <- file.path(PROJ, "mr")
TABO  <- file.path(PROJ, "Result", "tables_csv"); dir.create(TABO, showWarnings = FALSE, recursive = TRUE)
fmtP  <- function(p) ifelse(is.na(p), "—", ifelse(p == 0, "<1e-300", ifelse(p < 0.001, formatC(p, format = "e", digits = 2), sprintf("%.3g", p))))
RHO_UKB <- 0.10   # standardized LDSC cross-trait intercept COA–AOA (Methods)
RHO_FG  <- 0.50   # assumed for the FinnGen outcomes (shared controls; Step12/13)

## ======================= Panel A: rg contrasts (delta method) ======================= ##
L <- readRDS(file.path(LDSC, "ldsc_result.rds"))
S <- L$S; V <- L$V; k <- nrow(S)
TRAITS <- c("COA", "AOA", "MDD", "BIP", "SCZ", "ANX", "FEV1", "FVC", "LUNG")   # Step1 order
stopifnot(k == length(TRAITS)); rownames(S) <- colnames(S) <- TRAITS
vech <- function(M) M[lower.tri(M, diag = TRUE)]                     # column-major lower triangle = GenomicSEM V order
unvech <- function(v) { M <- matrix(0, k, k); M[lower.tri(M, diag = TRUE)] <- v; M <- M + t(M) - diag(diag(M)); dimnames(M) <- list(TRAITS, TRAITS); M }
rg_of <- function(v, a, b) { M <- unvech(v); M[a, b] / sqrt(M[a, a] * M[b, b]) }
## numerical gradient of a scalar function of vech(S), then delta-method SE
delta_se <- function(f, v = vech(S), eps = 1e-7) {
  g <- vapply(seq_along(v), function(i) { u <- v; u[i] <- u[i] + eps; d <- v; d[i] <- d[i] - eps; (f(u) - f(d)) / (2 * eps) }, numeric(1))
  as.numeric(sqrt(t(g) %*% V %*% g))
}
d_rg <- function(v, X) rg_of(v, "AOA", X) - rg_of(v, "COA", X)
psy <- c("MDD", "ANX", "BIP", "SCZ"); NAME <- c(MDD = "depression", ANX = "anxiety", BIP = "bipolar disorder", SCZ = "schizophrenia")
rowA <- function(label, est, se) data.table(Panel = "A. Genetic-correlation contrasts", Contrast = label,
  Quantity = if (grepl("^\u0394", label)) "Difference of \u0394r_g" else "\u0394r_g = r_g(AOA) \u2212 r_g(COA)",
  `Estimate (95% CI)` = sprintf("%.3f (%.3f to %.3f)", est, est - 1.96 * se, est + 1.96 * se),
  SE = sprintf("%.3f", se), P = fmtP(2 * pnorm(-abs(est / se))),
  `Method / assumption` = "Delta method, full LDSC sampling covariance (accounts for shared asthma controls)")
A1 <- rbindlist(lapply(psy, function(X) { f <- function(v) d_rg(v, X); rowA(sprintf("r_g(AOA, %s) \u2212 r_g(COA, %s)", NAME[X], NAME[X]), f(vech(S)), delta_se(f)) }))
A2 <- rbindlist(lapply(c("MDD", "ANX"), function(X) rbindlist(lapply(c("BIP", "SCZ"), function(Y) {
  f <- function(v) d_rg(v, X) - d_rg(v, Y); rowA(sprintf("\u0394r_g %s \u2212 \u0394r_g %s", NAME[X], NAME[Y]), f(vech(S)), delta_se(f)) }))))
panelA <- rbind(A1, A2)

## ======================= Panel B: MR contrasts ======================= ##
ivw_diff <- function(rds1, rds2, rho) {   # per-variant outcome differences regressed on exposure betas (as Step13)
  a <- as.data.table(readRDS(rds1)$dat); b <- as.data.table(readRDS(rds2)$dat)
  m <- merge(a[, .(SNP, bx = beta.exposure, b1 = beta.outcome, s1 = se.outcome)],
             b[, .(SNP, b2 = beta.outcome, s2 = se.outcome)], by = "SNP")
  m[, `:=`(d = b1 - b2, sd = sqrt(s1^2 + s2^2 - 2 * rho * s1 * s2))]
  fit <- lm(d ~ bx - 1, weights = 1 / sd^2, data = m)
  bb <- coef(summary(fit))[1, 1]; se <- coef(summary(fit))[1, 2] / min(1, summary(fit)$sigma)   # multiplicative random effects
  list(n = nrow(m), b = bb, se = se, p = 2 * pnorm(-abs(bb / se)))
}
ivw_est <- function(rds) { o <- readRDS(rds); r <- o$res[o$res$method == "Inverse variance weighted", ]; list(b = r$b, se = r$se) }
rowB <- function(label, qty, b, se, p, n, method) data.table(Panel = "B. MR contrasts", Contrast = label, Quantity = qty,
  `Estimate (95% CI)` = if (is.na(b)) "—" else sprintf("%.2f (%.2f\u2013%.2f)", exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se)),
  SE = if (is.na(b)) "—" else sprintf("%.3f", se), P = fmtP(p),
  `Method / assumption` = if (is.na(n)) method else sprintf("%s; %d shared instruments", method, n))

u <- ivw_diff(file.path(MRDIR, "MR_reverse_MDD_AOA.rds"), file.path(MRDIR, "MR_reverse_MDD_COA.rds"), RHO_UKB)
B1 <- rowB("Depression \u2192 adult-onset vs childhood-onset asthma (UK Biobank)", "Ratio of ORs (AOA / COA)", u$b, u$se, u$p, u$n,
           sprintf("IVW regression of per-variant differences, error correlation %.2f (LDSC cross-trait intercept), multiplicative random effects", RHO_UKB))
fg <- ivw_diff(file.path(MRDIR, "REP_HOWARD_FGASTHMA.rds"), file.path(MRDIR, "REP_HOWARD_FGCOA.rds"), RHO_FG)
B2 <- rowB("Depression \u2192 asthma (all) vs childhood asthma (FinnGen)", "Ratio of ORs (all / childhood)", fg$b, fg$se, fg$p, fg$n,
           sprintf("IVW regression of per-variant differences, assumed error correlation %.1f (shared controls); exploratory", RHO_FG))
m1 <- ivw_est(file.path(MRDIR, "MR_reverse_MDD_AOA.rds")); a1 <- ivw_est(file.path(MRDIR, "MR_reverse_ANX_AOA.rds"))
zd <- (m1$b - a1$b) / sqrt(m1$se^2 + a1$se^2)
B3 <- rowB("Depression vs anxiety \u2192 adult-onset asthma", "Ratio of ORs (depression / anxiety)", m1$b - a1$b, sqrt(m1$se^2 + a1$se^2), 2 * pnorm(-abs(zd)), NA_integer_,
           "z-test assuming independent estimates (different instruments, shared outcome); exploratory, approximate")
panelB <- rbind(B1, B2, B3)

out <- rbind(panelA, panelB)
fwrite(out, file.path(TABO, "TableS9_subtype_contrasts.csv"))
cat("\n=== Table E9 (two panels) ===\n"); print(out[, .(Panel, Contrast, `Estimate (95% CI)`, P)])
cat("\n확인: Panel A 1행 = 0.221 (SE 0.028, P 6.3e-15); Panel B 1행 = 1.07 (0.94–1.22, P .30); 2행 = 1.03 (0.86–1.22, P .77); 3행 P ≈ .36\n")
writeLines(capture.output(sessionInfo()), file.path(TABO, paste0("_step14_sessionInfo_", Sys.Date(), ".txt")))
