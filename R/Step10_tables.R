## Manuscript tables — reframed story. Emit formatted CSVs + assembled xlsx (R openxlsx).
## 실행: source("Step10_tables.R")  (프로젝트 결과 CSV들을 읽어 Result/ 에 표 저장)
suppressPackageStartupMessages({
  for(p in c("data.table","openxlsx")) if(!requireNamespace(p,quietly=TRUE)) install.packages(p,repos="https://cloud.r-project.org")
  library(data.table)})
B <- "PATH/TO/PROJECT"        # 프로젝트 루트
O <- file.path(B,"Result","tables_csv"); dir.create(O, showWarnings=FALSE, recursive=TRUE)
fmtP <- function(p) ifelse(is.na(p),"—", ifelse(p<0.001, formatC(p,format="e",digits=2), sprintf("%.3g",p)))
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
t3 <- mm[,.(Direction, Exposure=NAME[EXPOSURE], Outcome=NAME[OUTCOME], `N IV`=N_IV,
  `F`=sprintf("%.0f",F_MEAN), `OR (95% CI)`=sprintf("%.3f (%.3f–%.3f)",IVW_OR,IVW_LCI,IVW_UCI),
  `IVW P`=fmtP(IVW_P), `P (FDR)`=fmtP(IVW_P_FDR),
  `Egger int. P`=fmtP(EGGER_INT_P), `Q P`=fmtP(Q_P),
  `MR-PRESSO P`=PRESSO_P, `Steiger OK`=STEIGER_OK)]
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
tm[,susie:=fifelse(!is.na(su_h4), sprintf("H4=%.2f",su_h4), "not evaluable")]
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
fwrite(pl[, setdiff(names(pl), c("N_SIG_FDR","LAMBDA_scr")), with=FALSE], file.path(O,"TableS3_PLACO_summary.csv"))
# S4 blood immune / eos / crp MR
im <- fread(file.path(B,"mr/mr_immune_ALL.csv"))
s4 <- im[,.(Exposure=EXPOSURE,Outcome=OUTCOME,Domain=DOMAIN,`N IV`=N_IV,F=sprintf("%.0f",F_MEAN),
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

cat("TABLES CSV DONE\n"); print(list.files(O))

## ---------- assemble formatted workbook (openxlsx) ---------- ##
library(openxlsx)
SHEETS <- list(
 c("Table 1. Traits","Table1_traits.csv","GWAS datasets, heritability and LDSC intercept for the nine primary traits, and the exposure datasets used for multivariable MR and for replication."),
 c("Table 2. LDSC rg","Table2_LDSC_rg.csv","Genetic correlations of asthma subtypes with psychiatric and lung-function traits."),
 c("Table 3. MR main","Table3_MR_main.csv","Bidirectional two-sample MR (IVW) between asthma subtypes and depression/anxiety."),
 c("Table 4. MVMR","Table4_MVMR.csv","Multivariable MR of depression, BMI and smoking initiation on asthma subtypes (IVW-MVMR, MVMR-Egger, MVMR-median; conditional F; Q)."),
 c("Table 5. Immune coloc","Table5_immune_coloc.csv","Tissue cis-MR of candidate immune genes and colocalization (coloc, full window and +/-100 kb). SuSiE-based colocalization was not evaluable for any pair."),
 c("S1. LDSC rg all","TableS1_LDSC_rg_all.csv","All pairwise LDSC genetic correlations."),
 c("S2. LAVA loci","TableS2_LAVA_loci.csv","Significant local genetic correlations (LAVA, FDR<0.05) with mapped genes."),
 c("S3. PLACO","TableS3_PLACO_summary.csv","Exploratory PLACO/PLACO+ screen restricted to variants with P<0.001 in at least one trait; N_LOCI = distance-defined regions (500 kb) containing variants with P<5e-8."),
 c("S4. Blood immune MR","TableS4_blood_immune_MR.csv","Blood cis-eQTL / eosinophil / CRP MR."),
 c("S5. Tissue cisMR all","TableS5_tissue_cisMR_all.csv","All tissue cis-MR Wald-ratio estimates."),
 c("S6. Coloc all","TableS6_coloc_all.csv","Full colocalization posteriors (H0-H4), full region and +/-100kb."),
 c("S7. MR power","TableS7_MR_power.csv","Minimum OR detectable with 80% power (alpha 0.05) for bidirectional MR, per unit log-odds and per doubling."),
 c("S8. MR sensitivity","TableS8_MR_sensitivity.csv","All MR estimators, heterogeneity, leave-one-out, MR-PRESSO outliers and Steiger for the primary bidirectional MR. Steiger was computed on the liability scale (see Table S12)."),
 c("S9. Subtype contrasts","TableS9_subtype_contrasts_rg.csv","Adult- minus childhood-onset differences in genetic correlation and differences of these differences (delta method)."),
 c("S10. Lung MR","TableS10_lung_MR.csv","Pairwise MR between asthma subtypes or psychiatric traits and FEV1, FVC and FEV1/FVC (IVW; FDR within direction-specific family)."),
 c("S11. FinnGen replication","TableS11_FinnGen_replication.csv","Replication in FinnGen R13 (J10_ASTHMA_EXMORE, ASTHMA_CHILD_EXMORE) using the Howard 2019 depression GWAS; palindromic variants removed; MR-PRESSO with seed 20260911 and 10,000 simulations; Steiger recomputed on the liability scale. The last row is a ratio of odds ratios from an exploratory comparison assuming an error correlation of 0.5; the two outcomes share controls and are not mutually exclusive."),
 c("S12. Steiger sensitivity","TableS12_steiger_binary.csv","Steiger directionality recomputed on the liability scale using case and control counts and assumed population prevalences, across all exposure-outcome prevalence combinations (pre-specified)."))
wb <- createWorkbook()
hdr <- createStyle(fgFill="#1F3864", fontColour="white", textDecoration="bold",
                   fontName="Arial", fontSize=10, border="TopBottomLeftRight",
                   borderColour="#D9D9D9", wrapText=TRUE, halign="left", valign="top")
cel <- createStyle(fontName="Arial", fontSize=10, border="TopBottomLeftRight",
                   borderColour="#D9D9D9", wrapText=TRUE, valign="top")
ttl <- createStyle(fontName="Arial", fontSize=11, textDecoration="bold")
nte <- createStyle(fontName="Arial", fontSize=9, fontColour="#555555", textDecoration="italic")
addWorksheet(wb,"Legend", gridLines=FALSE)
leg <- c("Manuscript tables — asthma x depression/anxiety shared genetics",
 "Summary: genetic correlations of depression and anxiety were stronger with adult- than childhood-onset asthma; depression liability was associated with adult-onset asthma in MR, and in FinnGen with both asthma overall and childhood asthma, so an onset-specific causal effect was not established; shared regulatory signals for candidate immune genes were not established.",
 "","Abbreviations",
 "rg = genetic correlation; h2 = SNP heritability; IVW = inverse-variance weighted MR; FDR = Benjamini-Hochberg; PP.H4 = posterior prob. of a shared causal variant.",
 "COA = childhood-onset asthma; AOA = adult-onset asthma; MDD = depression; ANX = anxiety.",
 "","Notes",
 "Liability-scale h2 uses assumed population prevalences (see Methods).",
 "ORs for binary exposures are per unit increase in log-odds of exposure liability; full MR sensitivity results are in Table S8.")
for(i in seq_along(leg)){ writeData(wb,"Legend",leg[i],startRow=i,startCol=1)
  addStyle(wb,"Legend", if(i==1) ttl else if(leg[i] %in% c("Abbreviations","Notes")) ttl else nte, rows=i, cols=1)}
setColWidths(wb,"Legend",cols=1,widths=110)
for(s in SHEETS){ nm<-substr(s[1],1,31); f<-file.path(O,s[2]); if(!file.exists(f)) next
  d<-fread(f,colClasses="character"); setDF(d); d[is.na(d)]<-""
  addWorksheet(wb,nm,gridLines=FALSE)
  writeData(wb,nm,s[1],startRow=1,startCol=1); addStyle(wb,nm,ttl,rows=1,cols=1)
  writeData(wb,nm,s[3],startRow=2,startCol=1); addStyle(wb,nm,nte,rows=2,cols=1)
  writeData(wb,nm,d,startRow=4,startCol=1,headerStyle=hdr,borders="all",borderColour="#D9D9D9")
  addStyle(wb,nm,cel,rows=5:(4+nrow(d)),cols=1:ncol(d),gridExpand=TRUE)
  w<-pmin(pmax(nchar(names(d))+2, apply(d,2,function(x) max(nchar(x))+2)),52)
  setColWidths(wb,nm,cols=1:ncol(d),widths=w); freezePane(wb,nm,firstActiveRow=5) }
XLSX <- file.path(B,"Result","Manuscript_Tables_v4.xlsx")
saveWorkbook(wb, XLSX, overwrite=TRUE)
cat("WORKBOOK:", XLSX, "\n")
