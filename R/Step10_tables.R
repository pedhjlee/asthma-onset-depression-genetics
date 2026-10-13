## Manuscript tables — reframed story. Emit formatted CSVs + assembled xlsx (R openxlsx).
## 실행: source("Step10_tables.R")  (프로젝트 결과 CSV들을 읽어 Result/ 에 표 저장)
suppressPackageStartupMessages({
  for(p in c("data.table","openxlsx")) if(!requireNamespace(p,quietly=TRUE)) install.packages(p,repos="https://cloud.r-project.org")
  library(data.table)})
B <- "PATH/TO/PROJECT"        # 프로젝트 루트
O <- file.path(B,"Result","tables_csv"); dir.create(O, showWarnings=FALSE, recursive=TRUE)
## ★ 표 버전: 내용이 바뀔 때마다 올릴 것 (큰 변화 = 정수, 소수정·재계산 = 소수점). 파일명에 붙음.
TAB_VER <- "5.8"
TAB_NOTE <- "v5.8 (2026-10-10): all P values and scientific notation in E-tables converted to JACI style at workbook assembly (jaci_cells); E9 internal note removed; tissue labels spelled out. v5.7 (2026-10-08): P values in JACI style (two decimals, three when < .01, exponent notation below .001); Table E12 analysis labels spelled out; Table E4 Domain filled for eosinophil/CRP rows; version history removed from Legend (see HANDOFF)."
## v5.7: JACI P-value style — .xx (two decimals), .00x when < .01, a.b × 10⁻ⁿ when < .001, 0 → <1 × 10⁻³⁰⁰
SUP <- c("0"="⁰","1"="¹","2"="²","3"="³","4"="⁴","5"="⁵","6"="⁶","7"="⁷","8"="⁸","9"="⁹","-"="⁻")
sup <- function(x) vapply(strsplit(as.character(x),""), function(ch) paste(SUP[ch],collapse=""), "")
fmtP <- function(p){ p <- as.numeric(p); out <- rep("—", length(p)); ok <- !is.na(p)
  e <- floor(log10(abs(p[ok & p>0]))); m <- p[ok & p>0]/10^e
  expo <- paste0(sprintf("%.1f",m)," × 10",sup(e))
  out[ok & p==0] <- "<1 × 10⁻³⁰⁰"
  out[ok & p>0]  <- ifelse(p[ok & p>0] < 0.001, expo,
                    ifelse(p[ok & p>0] < 0.01, sub("^0","",sprintf("%.3f",p[ok & p>0])),
                    ifelse(p[ok & p>0] > 0.995, ">.99", sub("^0","",sprintf("%.2f",p[ok & p>0])))))
  out }
ORci <- function(or,l,u) sprintf("%.2f (%.2f–%.2f)", or,l,u)
NAME <- c(COA="Childhood-onset asthma",AOA="Adult-onset asthma",MDD="Depression (MDD)",
  ANX="Anxiety",BIP="Bipolar disorder",SCZ="Schizophrenia",FEV1="FEV1",FVC="FVC",
  LUNG="FEV1/FVC ratio",EOS="Eosinophil count",CRP="C-reactive protein")

## ---------- Table 1: traits & heritability ---------- ##
h2 <- fread(file.path(B,"ldsc/h2_results.csv")); setkey(h2,TRAIT)
# 확정된 소스(Step0 스크립트 + readme + 문헌). effective N = harmonized medN.
META <- data.table(
 T=c("COA","AOA","MDD","BIP","SCZ","ANX","FEV1","FVC","LUNG"),
 Cat=c("Asthma","Asthma","Psychiatric","Psychiatric","Psychiatric","Psychiatric",
       "Lung function","Lung function","Lung function"),
 Ref=c("Ferreira 2019, Am J Hum Genet (PMID 30929738)",
       "Ferreira 2019, Am J Hum Genet (PMID 30929738)",
       "PGC MDD, Adams 2025, Cell (PMID 39814019); EUR, no-UKB, no-23andMe",
       "PGC3 BIP, Mullins 2021, Nat Genet (PMID 34002096); no-UKB",
       "PGC3 SCZ, Trubetskoy 2022, Nature (PMID 35396580); EUR public release",
       "PGC-ANX, Strom 2026, Nat Genet (doi 10.1038/s41588-025-02485-8); woUTAH; includes UK Biobank",
       "Shrine 2019, Nat Genet (PMID 30804560); SpiroMeta",
       "Shrine 2019, Nat Genet (PMID 30804560); SpiroMeta",
       "Shrine 2019, Nat Genet (PMID 30804560); SpiroMeta"),
 CC=c("13,962 / 300,671","26,582 / 300,671","357,636 / 1,281,936","40,463 / 313,436",
      "52,017 / 75,889 (+1,369 trios)","122,083 / 729,602","75,676","75,423","75,639"),
 Neff=c("53,370","97,691","967,078","96,289","58,749","375,193","—","—","—"))
ORD1 <- c("COA","AOA","MDD","ANX","BIP","SCZ","FEV1","FVC","LUNG")
META <- META[match(ORD1, T)]
t1 <- data.table(
  Trait=NAME[META$T], Category=META$Cat, `Reference`=META$Ref,
  `Cases / controls (or N)`=META$CC, `Effective N`=META$Neff,
  Scale=h2[META$T]$SCALE,
  `SNP-h2 (SE)`=sprintf("%.3f (%.3f)",h2[META$T]$h2,h2[META$T]$h2_SE),
  `LDSC intercept`=sprintf("%.3f",h2[META$T]$intercept))
t1 <- rbind(t1, data.table(Trait=c("Body mass index","Smoking initiation"), Category="MVMR covariate",
  Reference=c("GIANT, Locke 2015, Nature (PMID 25673413); no-UKB","GSCAN, Liu 2019, Nat Genet (PMID 30643251); no-UKB, no-23andMe"),
  `Cases / controls (or N)`=c("up to 322,154","—"), `Effective N`=c("—","up to 249,176"),
  Scale="—", `SNP-h2 (SE)`="—", `LDSC intercept`="—"), fill=TRUE)
t1 <- rbind(t1, data.table(
  Trait=c("Asthma (FinnGen)","Childhood asthma (FinnGen)","Depression (replication exposure)"),
  Category="Replication",
  Reference=c("FinnGen R13, J10_ASTHMA_EXMORE (PMID 36653562)","FinnGen R13, ASTHMA_CHILD_EXMORE, first record age<16","Howard 2019, Nat Neurosci (PMID 30718901); PGC + UKB, no FinnGen"),
  `Cases / controls (or N)`=c("61,196 / 250,433","8,428 / 250,433","170,756 / 329,443"),
  `Effective N`=c("—","—","—"), Scale="—", `SNP-h2 (SE)`="—", `LDSC intercept`="—"), fill=TRUE)
fwrite(t1, file.path(O,"Table1_traits.csv"))

## ---------- Table 2: LDSC rg, asthma subtype x trait ---------- ##
rg <- fread(file.path(B,"ldsc/rg_results.csv"))
r2 <- rg[focus=="asthma_x_trait"]
r2[,Subtype:=NAME[T1]][,Trait:=NAME[T2]]
r2[,`rg (SE)`:=sprintf("%.3f (%.3f)",rg,SE)]
r2[,Z:=sprintf("%.2f",Z)][,P:=fmtP(P)][,`P (FDR)`:=fmtP(P_FDR)]
r2[,tord:=match(T2,c("MDD","ANX","BIP","SCZ","FEV1","FVC","LUNG"))]
r2[,sord:=match(T1,c("AOA","COA"))]
setorder(r2,sord,tord)
t2 <- r2[,.(Subtype,Trait,`rg (SE)`,Z,P,`P (FDR)`)]
fwrite(t2, file.path(O,"Table2_LDSC_rg.csv"))

## ---------- Table 3: MR bidirectional (main) ---------- ##
fw <- fread(file.path(B,"mr/mr_summary_forward.csv")); fw[,Direction:="Asthma → psychiatric"]
rv <- fread(file.path(B,"mr/mr_summary_reverse.csv")); rv[,Direction:="Psychiatric → asthma"]
mm <- rbindlist(list(fw,rv),fill=TRUE)
## MR-PRESSO P 와 Steiger 는 Step13 최종값(TableS8)에서 가져옴 (mr_summary_*.csv 의 값은 Step4 당시 1,000회·연속형 Steiger)
s8 <- fread(file.path(O,"TableS8_MR_sensitivity.csv"), colClasses="character")
## E8 의 t-검정 P 표시: 본문·각주 모두 .050 (원 계산값 .0503 은 Step13 출력에 보존)
s8[Exposure=="MDD" & Outcome=="COA" & `Outlier-corrected P`=="0.0503", `Outlier-corrected P`:="0.050"]
fwrite(s8, file.path(O,"TableS8_MR_sensitivity.csv"))
mm <- merge(mm, s8[,.(EXPOSURE=Exposure, OUTCOME=Outcome, ORCI8=`IVW OR (95% CI)`, PRESSO_P13=`MR-PRESSO global P`, STEIGER_P13=`Steiger P`)],
            by=c("EXPOSURE","OUTCOME"), all.x=TRUE, sort=FALSE)
mm[,dord:=match(Direction,c("Asthma → psychiatric","Psychiatric → asthma"))][,eord:=match(EXPOSURE,c("COA","AOA","MDD","ANX"))][,oord:=match(OUTCOME,c("MDD","ANX","COA","AOA"))]; setorder(mm,dord,eord,oord)
## JACI: OR·CI 는 소수 둘째 자리까지 — E8 과 동일 소스에서 가져옴 (v5.3)
t3 <- mm[,.(Direction, Exposure=NAME[EXPOSURE], Outcome=NAME[OUTCOME], `No. of instruments`=N_IV,
  `F`=sprintf("%.0f",F_MEAN), `OR (95% CI)`=ORCI8,   # E8(전정밀도) 문자열 그대로; mr_summary_*.csv 의 3자리 반올림값 사용 금지
  `IVW P`=fmtP(IVW_P), `P (FDR)`=fmtP(IVW_P_FDR),
  `Q P`=fmtP(Q_P), `Steiger P (liability)`=fmtP(ifelse(grepl("^<",STEIGER_P13),0,STEIGER_P13)))]   # v5.7: JACI P style (Step13 string → numeric)   # v5.4: Egger int. P, MR-PRESSO P → Table E8
fwrite(t3, file.path(O,"Table3_MR_main.csv"))

## ---------- Table 4: immune tissue cis-MR + colocalization ---------- ##
tm <- fread(file.path(B,"mr/mr_tissue_cis.csv"))[MATCHED==TRUE & TISSUE%in%c("Lung","Brain_Cortex")]
co <- fread(file.path(B,"coloc/coloc_summary.csv"))       # abf full
cn <- fread(file.path(B,"coloc/coloc_abf_narrow.csv"))    # abf narrow
su <- fread(file.path(B,"coloc/coloc_susie_summary.csv")) # susie
tm[,tis:=ifelse(TISSUE=="Lung","Lung","Brain cortex")]
key <- function(g,t,o) paste(g,t,o)
co[,k:=key(GENE,ifelse(TISSUE=="Lung","Lung","Brain_Cortex"),OUTCOME)]
cn[,k:=key(GENE,ifelse(TISSUE=="Lung","Lung","Brain_Cortex"),OUTCOME)]
su[,k:=key(GENE,ifelse(TISSUE=="Lung","Lung","Brain_Cortex"),OUTCOME)]
tm[,k:=key(GENE,TISSUE,OUTCOME)]
tm <- merge(tm, co[,.(k,ph4_full=PP_H4,ph3_full=PP_H3)], by="k", all.x=TRUE)
tm <- merge(tm, cn[,.(k,ph4_narrow=PP_H4)], by="k", all.x=TRUE)
tm <- merge(tm, su[,.(k,su_eqtl=nCS_eqtl,su_dis=nCS_dis,su_h4=maxPP_H4)], by="k", all.x=TRUE)
# 핵심 hit만: IL33[뇌]→MDD/ANX + 폐유전자→COA/AOA (coloc 있는 10행)
tm <- tm[(TISSUE=="Brain_Cortex" & GENE=="IL33" & OUTCOME%in%c("MDD","ANX")) |
         (TISSUE=="Lung" & OUTCOME%in%c("COA","AOA"))]
## v5.4: 미시도(IL1RL1: 참조 LD 수렴 실패, su 표에 없음) vs 시도했으나 평가 불가 구분
tm[,attempted:=k %in% su$k]
tm[,susie:=fifelse(!is.na(su_h4), sprintf("H4=%.2f",su_h4),
   fifelse(!attempted, "not attempted (LD-reference non-convergence)",
   fifelse(!is.na(su_eqtl) & su_eqtl==0 & !is.na(su_dis) & su_dis>0, "attempted; no eQTL credible set",
           "attempted; no evaluable credible-set pair")))]
tm[,tord:=match(OUTCOME,c("COA","AOA","MDD","ANX"))]
setorder(tm, TISSUE, GENE, tord)
t4 <- tm[,.(Gene=GENE, Tissue=tis, Outcome=NAME[OUTCOME],
  `cis-MR OR (95% CI)`=ORci(OR,OR_LCI,OR_UCI), `MR P`=fmtP(P), `MR P (FDR)`=fmtP(P_FDR),
  `coloc PP.H4 (full)`=sprintf("%.2f",ph4_full), `coloc PP.H4 (±100kb)`=sprintf("%.2f",ph4_narrow),
  `SuSiE-coloc`=susie)]
fwrite(t4, file.path(O,"Table5_immune_coloc.csv"))

## ---------- Supplementary ---------- ##
# S1 full rg
s1 <- rg[,.(Trait1=NAME[T1],Trait2=NAME[T2],rg=sprintf("%.3f",rg),SE=sprintf("%.3f",SE),
  Z=sprintf("%.2f",Z),P=fmtP(P),`P (FDR)`=fmtP(P_FDR),`gcov intercept`=sprintf("%.3f",gcov_intercept))]
fwrite(s1, file.path(O,"TableS1_LDSC_rg_all.csv"))
# S2 LAVA
lv <- fread(file.path(B,"lava/lava_sig_annotated.csv")); bv <- fread(file.path(B,"lava/lava_bivar.csv"))
lv <- merge(lv, unique(bv[,.(LOC,pair,rho.lower,rho.upper,p.Bonf)],by=c("LOC","pair")), by=c("LOC","pair"), all.x=TRUE, sort=FALSE)
setcolorder(lv, intersect(c("pair","LOC","CHR","START","STOP","rho","rho.lower","rho.upper","p","p.FDR","p.Bonf"), names(lv)))
cat("S2 rho CI 채워진 행:", sum(!is.na(lv$rho.lower)), "/", nrow(lv), "\n")
fwrite(lv, file.path(O,"TableS2_LAVA_loci.csv"))
# S3 PLACO summary
pl <- fread(file.path(B,"placo/placo_summary.csv"))
## v5.4: 독자용 열만, 설명적 제목 (MIN = 계산 시간(분), N_NA = 전부 0 → 제외)
pl[,eord:=match(EXPOSURE,c("COA","AOA"))][,oord:=match(OUTCOME,c("MDD","ANX","BIP","SCZ","LUNG"))]; setorder(pl,eord,oord)
s3 <- pl[,.(Exposure=NAME[EXPOSURE], Outcome=NAME[OUTCOME], Method=METHOD,
  `Z-score correlation (CorZ)`=sprintf("%.4f",CORZ), `Variants tested`=format(N_TESTED,big.mark=","),
  `Variants with P < 5 × 10⁻⁸`=format(N_SIG_GW,big.mark=","), `Regions (500-kb rule)`=N_LOCI, `Minimum P`=fmtP(MIN_P))]
fwrite(s3, file.path(O,"TableS3_PLACO_summary.csv"))
# S4 blood immune / eos / crp MR
im <- fread(file.path(B,"mr/mr_immune_ALL.csv"))
## v5.7: EOS/CRP rows had empty Domain → "biomarker"
s4 <- im[,.(Exposure=EXPOSURE,Outcome=OUTCOME,Domain=ifelse(is.na(DOMAIN)|DOMAIN=="","biomarker",DOMAIN),`N IV`=N_IV,F=sprintf("%.0f",F_MEAN),
  Method=METHOD,`OR (95% CI)`=ORci(OR,OR_LCI,OR_UCI),P=fmtP(P),`P (FDR)`=fmtP(P_FDR))]
fwrite(s4, file.path(O,"TableS4_blood_immune_MR.csv"))
# S5 tissue cis-MR full
s5 <- tm2 <- fread(file.path(B,"mr/mr_tissue_cis.csv"))
fwrite(s5, file.path(O,"TableS5_tissue_cisMR_all.csv"))
# S6 coloc all
s6 <- merge(co[,.(GENE,TISSUE,OUTCOME,nSNP,PP_H0,PP_H1,PP_H2,PP_H3,PP_H4)],
   cn[,.(GENE,TISSUE,OUTCOME,PP_H3_narrow=PP_H3,PP_H4_narrow=PP_H4)],
   by=c("GENE","TISSUE","OUTCOME"), all=TRUE)
fwrite(s6, file.path(O,"TableS6_coloc_all.csv"))
# S7 power
pw <- fread(file.path(B,"mr/mr_power_MDE.csv")); fwrite(pw, file.path(O,"TableS7_MR_power.csv"))

# Table IV (본문): joint model 만 / Table E14: MVMR 전체 (부분모형 포함)
t4 <- fread(file.path(O,"Table4_MVMR.csv"), colClasses="character")
fwrite(t4, file.path(O,"TableE14_MVMR_full.csv"))
fwrite(t4[Model=="Depression + BMI + smoking initiation"], file.path(O,"Table4_MVMR_main.csv"))
cat("TABLES CSV DONE\n"); print(list.files(O))

## ---------- assemble formatted workbooks (openxlsx) ---------- ##
## 두 파일: (1) 본문 표 Table I–IV  (2) Online Repository 표 Table E1–E13. 파일명에 TAB_VER 가 붙음.
library(openxlsx)
MAIN <- list(
 c("Table I","Table1_traits.csv","GWAS datasets, heritability and LDSC intercepts for the primary traits, and the exposure datasets used for multivariable MR and replication."),
 c("Table II","Table2_LDSC_rg.csv","Genome-wide genetic correlations of asthma subtypes with psychiatric and lung-function traits (LDSC; FDR across 36 pairs)."),
 c("Table III","Table3_MR_main.csv","Bidirectional two-sample MR (IVW) between asthma subtypes and depression/anxiety. MR-PRESSO global P (seed 20260911, 10,000 simulations) and liability-scale Steiger P are the final values (Table E8)."),
 c("Table IV","Table4_MVMR_main.csv","Multivariable MR of depression, BMI and smoking initiation on asthma subtypes, joint model (IVW-MVMR, MVMR-Egger, MVMR-median, Q-minimization; conditional F; Q). Conditional F statistics assume zero covariance between exposure-association estimates, except for Q-minimization rows, for which conditional F is computed under the middle exposure-correlation scenario. Q-minimization does not provide a P value. Partial models are in Table E14."))
OR <- list(
 c("Table E1","TableS1_LDSC_rg_all.csv","All pairwise LDSC genetic correlations."),
 c("Table E2","TableS2_LAVA_loci.csv","Significant local genetic correlations (LAVA, FDR<0.05) with mapped genes."),
 c("Table E3","TableS3_PLACO_summary.csv","Exploratory PLACO/PLACO+ screen restricted to variants with P<0.001 in at least one trait. PLACO+ was used when |CorZ| > 0.02. Regions = distance-defined regions (500-kb rule) containing variants with PLACO P<5e-8."),
 c("Table E4","TableS4_blood_immune_MR.csv","Blood cis-eQTL / eosinophil / CRP MR."),
 c("Table E5","TableS5_tissue_cisMR_all.csv","All tissue cis-MR Wald-ratio estimates."),
 c("Table E6","TableS6_coloc_all.csv","Colocalization posteriors: PP.H0-H4 for the full region (gene +/-500 kb); PP.H3 and PP.H4 for the +/-100-kb window around the lead eQTL."),
 c("Table E7","TableS7_MR_power.csv","Minimum OR detectable with 80% power (alpha 0.05) for bidirectional MR, per unit log-odds and per doubling."),
 c("Table E8","TableS8_MR_sensitivity.csv","All MR estimators, heterogeneity, leave-one-out, MR-PRESSO (seed 20260911, 10,000 simulations; outlier-corrected CI and P on the t distribution with df = retained instruments − 1) and liability-scale Steiger for the primary bidirectional MR. Footnote: for depression → childhood-onset asthma the outlier-corrected t-test P is .050 (lower 95% limit 0.99987, displayed as 1.00); a post hoc normal approximation gives P = .048 (95% CI 1.001–1.255). The t-based inference is primary."),
 c("Table E9","TableS9_subtype_contrasts.csv","Subtype contrasts. Panel A: adult- minus childhood-onset differences in genetic correlation and differences of these differences (delta method with the full LDSC sampling covariance matrix, which accounts for the shared asthma controls). Panel B: MR contrasts — depression → adult-onset vs childhood-onset asthma (UK Biobank; IVW regression of per-variant differences, error correlation 0.10, multiplicative random effects), depression → asthma overall vs childhood asthma (FinnGen; assumed error correlation 0.5, exploratory; also in Table E11), and depression vs anxiety → adult-onset asthma (z-test assuming independence, exploratory)."),
 c("Table E10","TableS10_lung_MR.csv","Pairwise MR between asthma subtypes or psychiatric traits and FEV1, FVC and FEV1/FVC (IVW; FDR within direction-specific family)."),
 c("Table E11","TableS11_FinnGen_replication.csv","Replication in FinnGen R13 (J10_ASTHMA_EXMORE, ASTHMA_CHILD_EXMORE) using the Howard 2019 depression GWAS; palindromic variants removed; MR-PRESSO with seed 20260911 and 10,000 simulations; Steiger on the liability scale. The last row is a ratio of odds ratios from an exploratory comparison assuming an error correlation of 0.5; the two outcomes share controls and are not mutually exclusive."),
 c("Table E12","TableS12_steiger_binary.csv","Steiger directionality recomputed on the liability scale using case and control counts scaled to the effective sample size and assumed population prevalences, across nine exposure-outcome prevalence combinations fixed before the liability-scale recomputation. Analysis labels: asthma → psychiatric = asthma liability as exposure; psychiatric → asthma = depression or anxiety liability as exposure (the primary direction in the manuscript)."),
 c("Table E13","Table5_immune_coloc.csv","Tissue cis-MR of candidate immune genes and colocalization (coloc, full window and +/-100 kb). SuSiE-based colocalization was attempted for eight combinations and yielded no evaluable shared-signal posterior; it was not attempted for IL1RL1 because the LD reference did not converge."),
 c("Table E14","TableE14_MVMR_full.csv","Multivariable MR, all models: joint model and partial models (depression + BMI; depression + smoking initiation) for adult-onset and childhood-onset asthma."))
hdr <- createStyle(fgFill="#1F3864", fontColour="white", textDecoration="bold", fontName="Arial", fontSize=10,
                   border="TopBottomLeftRight", borderColour="#D9D9D9", wrapText=TRUE, halign="left", valign="top")
cel <- createStyle(fontName="Arial", fontSize=10, border="TopBottomLeftRight", borderColour="#D9D9D9", wrapText=TRUE, valign="top")
ttl <- createStyle(fontName="Arial", fontSize=11, textDecoration="bold")
nte <- createStyle(fontName="Arial", fontSize=9, fontColour="#555555", textDecoration="italic")

## v5.8: 워크북 조립 시 모든 셀을 JACI 표기로 통일 — 과학적 표기(5.75e-06, <1e-04, 1.0 × 10^-15)를 위첨자 지수로,
##       P 열의 선행 0 제거(0.29 → .29). Step11/13 등 상위 단계 CSV 는 건드리지 않는다.
sci2sup <- function(x, keep_int=FALSE){   # keep_int: 정수 가수(5e-8)는 "5 × 10⁻⁸" 로 (각주·제목용)   # "5.75e-06" / "<1e-04" / "1.0 × 10^-15" 토큰 → "5.8 × 10⁻⁶" / "<1 × 10⁻⁴" / "1.0 × 10⁻¹⁵"
  f <- function(m, e, lt=""){ e <- as.integer(e); m <- as.numeric(m)
    ms <- if((lt=="<" && m==1) || (keep_int && !grepl("\\.", m))) as.character(as.integer(m)) else sprintf("%.1f", as.numeric(m))
    paste0(lt, ms, " × 10", sup(e)) }
  x <- gsub("(<?)(\\d+(?:\\.\\d+)?)[eE]([-+]?\\d+)", "\\1\\2e\\3", x, perl=TRUE)   # normalise
  x <- gsub("(\\d+(?:\\.\\d+)?) × 10\\^(-?\\d+)", "\\1e\\2", x, perl=TRUE)
  m <- gregexpr("(<?)(\\d+(?:\\.\\d+)?)e([-+]?\\d+)", x, perl=TRUE)
  regmatches(x, m) <- lapply(regmatches(x, m), function(tok) vapply(tok, function(t){
    lt <- ifelse(startsWith(t,"<"),"<",""); t2 <- sub("^<","",t); p <- strsplit(t2,"e")[[1]]; f(p[1],p[2],lt) }, ""))
  x }
jaci_cells <- function(d){
  isP <- grepl("(^|[^A-Za-z])[Pp]([^A-Za-z]|$)", names(d)) & !grepl("PP\\.H|Pos|Panel|PLACO", names(d))   # P, p, P_FDR, IVW P, Cochran Q (P) …
  for(j in seq_along(d)){ x <- d[[j]]
    if(!is.character(x)) next
    sci <- grepl("\\de[-+]?\\d|× 10\\^", x, perl=TRUE)
    if(isP[j]){
      x[sci] <- sci2sup(x[sci])
      num <- !sci & grepl("^0\\.\\d+$", x);  x[num] <- fmtP(x[num])          # 0.29 → .29, 0.009 → .009
      x <- gsub("\\(0\\.(\\d+)\\)", "(.\\1)", x, perl=TRUE)                # "344.4 (0.29)" 형태
    } else {
      pure <- sci & grepl("^-?\\d+(\\.\\d+)?e[-+]?\\d+$", x, perl=TRUE)       # 효과크기 등 순수 숫자 → 일반 소수
      x[pure] <- format(as.numeric(x[pure]), scientific=FALSE, digits=3, trim=TRUE)
      x[sci & !pure] <- sci2sup(x[sci & !pure])
    }
    x <- gsub("Whole_Blood","Whole blood",x,fixed=TRUE); x <- gsub("Brain_Cortex","Brain cortex",x,fixed=TRUE)
    d[[j]] <- x }
  d }
build_wb <- function(SHEETS, title, file) {
  wb <- createWorkbook()
  addWorksheet(wb,"Legend", gridLines=FALSE)
  leg <- c(title, paste0("Version ", TAB_VER, " — generated ", format(Sys.Date(), "%Y-%m-%d")), "",   # v5.7: version history kept in HANDOFF, not in the workbook
   "Abbreviations",
   "rg = genetic correlation; h2 = SNP heritability; IVW = inverse-variance weighted MR; FDR = Benjamini-Hochberg; PP.H4 = posterior prob. of a shared causal variant.",
   "COA = childhood-onset asthma; AOA = adult-onset asthma; MDD = depression; ANX = anxiety.",
   "", "Notes",
   "Liability-scale h2 uses assumed population prevalences (see Methods).",
   "ORs for binary exposures are per unit increase in log-odds of exposure liability.")
  for(i in seq_along(leg)){ writeData(wb,"Legend",leg[i],startRow=i,startCol=1)
    addStyle(wb,"Legend", if(i==1 || leg[i] %in% c("Abbreviations","Notes")) ttl else nte, rows=i, cols=1)}
  setColWidths(wb,"Legend",cols=1,widths=110)
  for(s in SHEETS){ nm<-substr(s[1],1,31); f<-file.path(O,s[2]); if(!file.exists(f)) { warning("missing: ", f); next }
    d<-fread(f,colClasses="character"); setDF(d); d[is.na(d)]<-""
    d <- jaci_cells(d)   # v5.8
    if(nm=="Table E12" && "Analysis" %in% names(d)){   # v5.7: spell out direction labels
      a <- d$Analysis
      a <- sub("^MR_main_([A-Z]+)_([A-Z]+)$",    "Asthma \u2192 psychiatric: \\1 \u2192 \\2", a)
      a <- sub("^MR_reverse_([A-Z]+)_([A-Z]+)$", "Psychiatric \u2192 asthma: \\1 \u2192 \\2", a)
      a <- sub("^REP_HOWARD_FGASTHMA$", "Replication: depression \u2192 FinnGen asthma (all)", a)
      a <- sub("^REP_HOWARD_FGCOA$",    "Replication: depression \u2192 FinnGen childhood asthma", a)
      d$Analysis <- a }
    if(nm=="Table E11" && "Outcome" %in% names(d)) d$Outcome <- sub("FinnGen childhood-onset asthma", "FinnGen childhood asthma", d$Outcome)   # v5.7: match manuscript term
    addWorksheet(wb,nm,gridLines=FALSE)
    writeData(wb,nm,s[1],startRow=1,startCol=1); addStyle(wb,nm,ttl,rows=1,cols=1)
    writeData(wb,nm,sci2sup(s[3],keep_int=TRUE),startRow=2,startCol=1); addStyle(wb,nm,nte,rows=2,cols=1)
    writeData(wb,nm,d,startRow=4,startCol=1,headerStyle=hdr,borders="all",borderColour="#D9D9D9")
    addStyle(wb,nm,cel,rows=5:(4+nrow(d)),cols=1:ncol(d),gridExpand=TRUE)
    w<-pmin(pmax(nchar(names(d))+2, apply(d,2,function(x) max(nchar(x))+2)),52)
    setColWidths(wb,nm,cols=1:ncol(d),widths=w); freezePane(wb,nm,firstActiveRow=5) }
  saveWorkbook(wb, file, overwrite=TRUE); cat("WORKBOOK:", file, "\n"); invisible(file)
}
build_wb(MAIN, "Main tables (Table I-IV) — asthma onset x depression/anxiety shared genetics",
         file.path(B,"Result",paste0("Manuscript_Tables_Main_v",TAB_VER,".xlsx")))
build_wb(OR,   "Online Repository tables (Table E1-E14) — asthma onset x depression/anxiety shared genetics",
         file.path(B,"Result",paste0("Manuscript_Tables_OnlineRepository_v",TAB_VER,".xlsx")))
