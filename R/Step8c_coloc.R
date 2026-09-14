###############################################################################
## Step 8c — Colocalization (eQTL Catalogue full sumstats × 질환 GWAS)
## Project: 01_GWAS_Asthma_MDD  (면역축: 조직특이 공유변이 확정)
##
## 목적: Step8b 조직 cis-MR 에서 유의했던 신호가 "같은 인과변이"를 공유하는지 검정.
##   cis-MR 단일도구(Wald)는 LD로 인한 가짜연관 가능성 → coloc(PP.H4)로 확정.
##
## 대상 hit (조직·질환 매칭):
##   IL33   [Brain_Cortex, QTD000171] → MDD, ANX     (혈액서 안보이던 신호)
##   TSLP   [Lung,         QTD000271] → COA, AOA
##   IL4R   [Lung,         QTD000271] → COA, AOA
##   IL1RL1 [Lung,         QTD000271] → COA, AOA
##   TYK2   [Lung,         QTD000271] → COA, AOA
##
## 입력:
##   raw/QTD000171.all.tsv.gz(+.tbi)  = Brain cortex (GTEx) 전영역 eQTL
##   raw/QTD000271.all.tsv.gz(+.tbi)  = Lung        (GTEx) 전영역 eQTL
##   harmonized/<TR>.rds              = 질환 GWAS (SNP=rsid, BETA/SE/EAF/P/N/EA/NEA)
##
## 방법: 유전자 cis영역(±PAD) tabix → gene_id 필터 → 질환과 rsid 매칭 → coloc.abf.
##   D1 eQTL: type="quant", beta/varbeta, MAF, sdY=1 (정규화 발현).
##   D2 질환: type="cc",   beta/varbeta, MAF, N, s(=case비율).
##   (coloc.abf 의 ABF는 z²기반 → allele flip이 PP.H4에 영향 없음: rsid 매칭이면 충분.)
##
## 산출: coloc/coloc_summary.csv, coloc/fig_coloc_pph4.png,
##       coloc/locus_<GENE>_<TR>.png (eQTL vs 질환 -log10P 산점)
## 실행: source("Step8c_coloc.R")   (Step8b 이후, QTD*.all.tsv.gz 로컬 필요)
###############################################################################

suppressPackageStartupMessages({
  for (p in c("data.table","coloc","ggplot2","BiocManager"))
    if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")
  ## Rsamtools/GenomicRanges: seqminer 은 대용량 bgzip 뒤쪽 오프셋(문자열정렬 chr2~9)을
  ## 못 읽는 윈도우 버그가 있어 htslib 기반 Rsamtools 로 tabix read (64bit 정상).
  for (p in c("Rsamtools","GenomicRanges","IRanges"))
    if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, update=FALSE, ask=FALSE)
  library(data.table); library(coloc); library(ggplot2)
  library(Rsamtools); library(GenomicRanges); library(IRanges)
})

PROJ <- "PATH/TO/PROJECT"
HARM <- file.path(PROJ,"harmonized")
RAW  <- file.path(PROJ,"raw")
OUTD <- file.path(PROJ,"coloc"); dir.create(OUTD, showWarnings=FALSE, recursive=TRUE)
setwd(OUTD)

## ---- hit 정의 (유전자 GRCh38 좌표: cis창 ±PAD) --------------------------- ##
##  좌표는 gene body (Ensembl GRCh38). eQTL Catalogue chromosome 은 "chr" 없이 "9" 등.
PAD <- 500000L
HITS <- list(
  list(gene="IL33",   ensg="ENSG00000137033", tissue="Brain_Cortex", qtd="QTD000171",
       chr=9L,  s=6215786L,   e=6257983L,   n_eqtl=205, dis=c("MDD","ANX")),
  list(gene="TSLP",   ensg="ENSG00000145777", tissue="Lung",         qtd="QTD000271",
       chr=5L,  s=110405682L, e=110413722L, n_eqtl=515, dis=c("COA","AOA")),
  list(gene="IL4R",   ensg="ENSG00000077238", tissue="Lung",         qtd="QTD000271",
       chr=16L, s=27313974L,  e=27364403L,  n_eqtl=515, dis=c("COA","AOA")),
  list(gene="IL1RL1", ensg="ENSG00000115602", tissue="Lung",         qtd="QTD000271",
       chr=2L,  s=102311533L, e=102352175L, n_eqtl=515, dis=c("COA","AOA")),
  list(gene="TYK2",   ensg="ENSG00000105397", tissue="Lung",         qtd="QTD000271",
       chr=19L, s=10350529L,  e=10380676L,  n_eqtl=515, dis=c("COA","AOA"))
)

## case 비율 s (coloc type="cc"): COA/AOA 는 이 프로젝트 case/(case+ctrl)
S_CASE <- c(COA=0.044, AOA=0.081, MDD=0.33, ANX=0.12)

## ---- eQTL Catalogue 스키마 (헤더 자동 확인, 실패시 표준 순서) ------------ ##
STD_COLS <- c("molecular_trait_id","chromosome","position","ref","alt","variant",
              "ma_samples","maf","pvalue","beta","se","type","ac","an","r2",
              "molecular_trait_object_id","gene_id","median_tpm","rsid")

get_schema <- function(f){
  h <- try(names(fread(f, nrows=0, showProgress=FALSE)), silent=TRUE)
  if (inherits(h,"try-error") || !length(h) || h[1]=="V1") return(NULL)
  h
}

## ---- 지정 cis영역 tabix 읽기 (Rsamtools::scanTabix) → data.table --------- ##
read_region <- function(qtd, chr, start, end){
  f <- file.path(RAW, paste0(qtd, ".all.tsv.gz"))
  stopifnot(file.exists(f), file.exists(paste0(f,".tbi")))
  schema <- get_schema(f); if (is.null(schema)) schema <- STD_COLS
  ## 인덱스에 든 염색체 이름에 맞춰 "5" / "chr5" 자동 선택
  seqs <- try(Rsamtools::seqnamesTabix(f), silent=TRUE)
  cc   <- as.character(chr)
  if (!inherits(seqs,"try-error") && !(cc %in% seqs) && paste0("chr",cc) %in% seqs) cc <- paste0("chr",cc)
  gr <- GenomicRanges::GRanges(cc, IRanges::IRanges(start=start, end=end))
  ln <- try(Rsamtools::scanTabix(Rsamtools::TabixFile(f), param=gr)[[1]], silent=TRUE)
  if (inherits(ln,"try-error") || !length(ln)) return(NULL)
  dt <- fread(text=paste(ln, collapse="\n"), header=FALSE, sep="\t", showProgress=FALSE)
  setnames(dt, schema[seq_len(ncol(dt))])
  dt
}

## ---- 한 hit(유전자) × 한 질환 coloc -------------------------------------- ##
run_coloc <- function(eq, h, tr){
  fr <- file.path(HARM, paste0(tr,".rds"))
  if (!file.exists(fr)) { message("    [", tr, "] harmonized 없음 → skip"); return(NULL) }
  dis <- as.data.table(readRDS(fr))

  e <- eq[!is.na(rsid) & rsid!="" & is.finite(beta) & is.finite(se) & se>0 &
          is.finite(maf) & maf>0 & maf<0.5,
          .(rsid, b_e=beta, v_e=se^2, maf_e=maf, p_e=pvalue)]
  e <- e[!duplicated(rsid)]

  d <- dis[SNP %in% e$rsid,
           .(rsid=SNP, b_d=BETA, v_d=SE^2, eaf=EAF, N_d=N, p_d=P)]
  d <- d[is.finite(b_d) & is.finite(v_d) & v_d>0 & is.finite(eaf) & eaf>0 & eaf<1]
  d[, maf_d := pmin(eaf, 1-eaf)]
  d <- d[!duplicated(rsid)]

  m <- merge(e, d, by="rsid")
  m <- m[maf_d>0 & maf_d<0.5]
  if (nrow(m) < 20) { message(sprintf("    [%s→%s] 공통 SNP %d개(<20) → skip", h$gene, tr, nrow(m))); return(NULL) }

  D1 <- list(beta=m$b_e, varbeta=m$v_e, snp=m$rsid, type="quant",
             sdY=1, MAF=m$maf_e, N=h$n_eqtl)
  D2 <- list(beta=m$b_d, varbeta=m$v_d, snp=m$rsid, type="cc",
             MAF=m$maf_d, N=round(median(m$N_d)), s=unname(S_CASE[tr]))

  res <- try(coloc.abf(D1, D2), silent=TRUE)
  if (inherits(res,"try-error")) { message("    coloc 실패: ", h$gene,"→",tr); return(NULL) }
  s <- as.list(res$summary)

  ## locus 산점도 (eQTL vs 질환 -log10P)
  m[, `:=`(logP_e = -log10(pmax(p_e,1e-300)), logP_d = -log10(pmax(p_d,1e-300)))]
  gp <- ggplot(m, aes(logP_e, logP_d)) +
    geom_point(alpha=.5, size=1.4, color="#2166AC") +
    labs(title=sprintf("%s [%s] eQTL vs %s  (PP.H4=%.2f, nSNP=%d)",
                       h$gene, h$tissue, tr, s$PP.H4.abf, s$nsnps),
         x=paste0("-log10 P  eQTL(",h$gene,")"), y=paste0("-log10 P  ",tr)) +
    theme_minimal(base_size=10)
  ggsave(sprintf("locus_%s_%s.png", h$gene, tr), gp, width=5, height=4.2, dpi=200, bg="white")

  data.table(GENE=h$gene, ENSG=h$ensg, TISSUE=h$tissue, QTD=h$qtd, OUTCOME=tr,
             nSNP=s$nsnps,
             PP_H0=round(s$PP.H0.abf,3), PP_H1=round(s$PP.H1.abf,3),
             PP_H2=round(s$PP.H2.abf,3), PP_H3=round(s$PP.H3.abf,3),
             PP_H4=round(s$PP.H4.abf,3),
             H4_over_H3=round(s$PP.H4.abf/max(s$PP.H3.abf,1e-6),2))
}

## ---- 실행 ---------------------------------------------------------------- ##
message("=== Step 8c: colocalization (coloc.abf) ===")
ALL <- list()
for (h in HITS){
  message(sprintf("-- %s [%s] cis %d:%d-%d (±%dkb)", h$gene, h$tissue,
                  h$chr, h$s, h$e, PAD/1000))
  eq <- read_region(h$qtd, h$chr, h$s - PAD, h$e + PAD)
  if (is.null(eq) || !nrow(eq)) { message("   tabix 결과 없음 → skip"); next }
  ## 대상 유전자만 (버전접미 무시하고 startsWith)
  eq <- eq[startsWith(as.character(gene_id), h$ensg)]
  message(sprintf("   %s eQTL SNP: %d", h$ensg, nrow(eq)))
  if (!nrow(eq)) { message("   해당 gene_id 없음 → skip"); next }
  for (tr in h$dis){
    r <- run_coloc(eq, h, tr)
    if (!is.null(r)) { ALL[[length(ALL)+1]] <- r
      message(sprintf("   ★ %s→%s : PP.H4=%.3f (nSNP=%d)", h$gene, tr, r$PP_H4, r$nSNP)) }
  }
}

if (length(ALL)){
  RES <- rbindlist(ALL, fill=TRUE)
  setorder(RES, -PP_H4)
  fwrite(RES, "coloc_summary.csv")
  message("\n=== coloc 요약 (PP.H4 내림차순) ===")
  print(RES[, .(GENE,TISSUE,OUTCOME,nSNP,PP_H3,PP_H4,H4_over_H3)])

  ## PP.H4 막대그림
  RES[, lab := paste0(GENE," [",substr(TISSUE,1,4),"] → ",OUTCOME)]
  RES[, lab := factor(lab, levels=rev(lab))]
  RES[, strong := PP_H4>=0.8]
  p <- ggplot(RES, aes(PP_H4, lab)) +
    geom_col(aes(fill=strong), width=.65, color="grey30") +
    geom_vline(xintercept=0.8, linetype=2, color="#B2182B") +
    geom_text(aes(label=sprintf("%.2f", PP_H4)), hjust=-0.15, size=3) +
    scale_fill_manual(values=c(`TRUE`="#B2182B",`FALSE`="grey75"), guide="none") +
    scale_x_continuous(limits=c(0,1.08), breaks=seq(0,1,.2)) +
    labs(title="Colocalization: shared causal variant (PP.H4)",
         subtitle="dashed = 0.8 threshold; red = colocalized",
         x="PP.H4 (posterior prob. shared causal variant)", y=NULL) +
    theme_minimal(base_size=10)
  ggsave("fig_coloc_pph4.png", p, width=7.5, height=4.5, dpi=300, bg="white")
  message("  fig_coloc_pph4.png 저장")

  strong <- RES[PP_H4>=0.8]
  message("\n★ 공유변이(PP.H4≥0.8): ",
          if(nrow(strong)) paste(sprintf("%s[%s]→%s(%.2f)", strong$GENE,
              substr(strong$TISSUE,1,4), strong$OUTCOME, strong$PP_H4), collapse=" · ") else "없음")
  message("  (0.5≤PP.H4<0.8 = 시사적 / <0.5 = 근거약함 → cis-MR 신호는 LD 주의)")
} else {
  message("\n[결과 없음] tabix/매칭 실패 — QTD 파일 경로·.tbi·염색체표기 확인 필요.")
}

message("\n================= Step 8c 완료 =================")
message("산출: coloc/coloc_summary.csv, coloc/fig_coloc_pph4.png, coloc/locus_*.png")
message("해석: PP.H4 높음 = eQTL과 질환이 '같은 인과변이' 공유 → cis-MR 신호 = 진짜 공유 면역기전.")
message("특히 IL33[뇌]↔우울/불안 이 PP.H4≥0.8 이면: 조직특이 IL33-ST2 alarmin 축이")
message("       천식(폐)·우울(뇌)에 걸쳐 공유되는 핵심 발견을 coloc가 확정.")
