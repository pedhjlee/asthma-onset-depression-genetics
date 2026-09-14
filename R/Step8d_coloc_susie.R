###############################################################################
## Step 8d — SuSiE-coloc (다중 인과변이) + 좁은창 coloc.abf 민감도
## Project: 01_GWAS_Asthma_MDD  (면역축: 조직특이 공유변이 엄밀검정)
##
## 배경(Step8c): coloc.abf(단일 인과변이 가정)가 IL1RL1 등에서 H3=1.0 →
##   이 loci는 독립신호 여러 개(allelic heterogeneity)라 abf가 자동으로 H3로 몰림.
##   → 다중 인과변이를 허용하는 coloc.susie 로 재검정(신호쌍별 PP.H4).
##   병행: lead eQTL ±NARROW 로 좁힌 창에서 coloc.abf (단일신호 근사) 민감도.
##
## 방법:
##   1) 유전자 cis (Rsamtools tabix) eQTL + 질환 harmonized 를 rsid로 교차.
##   2) EUR 참조패널(PATH/TO/ldref_ascii/EUR)에서 LD(부호 있는 r) 행렬(ieugwasr).
##   3) eQTL·질환 beta 와 LD 를 '동일 effect allele(=LD A1)'로 정렬(부호 맞춤).
##   4) runsusie(각 형질) → coloc.susie → 신호쌍별 PP.H4 (max 보고).
##   5) (민감도) lead eQTL ±NARROW coloc.abf.
##
## 대상: 폐 IL1RL1·TSLP·IL4R·TYK2 → COA/AOA (검정력 충분),
##       뇌 IL33 → MDD/ANX (검정력 한계 예상 — 신호 없으면 그대로 보고).
## 산출: coloc/coloc_susie_summary.csv, coloc/coloc_abf_narrow.csv,
##       coloc/fig_susie_pph4.png
## 실행: source("Step8d_coloc_susie.R")   (Step8c 이후)
###############################################################################

suppressPackageStartupMessages({
  need_cran <- c("data.table","coloc","susieR","ggplot2","ieugwasr","dplyr","BiocManager","remotes")
  for (p in need_cran) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")
  for (p in c("Rsamtools","GenomicRanges","IRanges")) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, update=FALSE, ask=FALSE)
  library(data.table); library(coloc); library(susieR); library(ggplot2)
  library(ieugwasr); library(Rsamtools); library(GenomicRanges); library(IRanges)
})

PROJ <- "PATH/TO/PROJECT"
HARM <- file.path(PROJ,"harmonized"); RAW <- file.path(PROJ,"raw")
OUTD <- file.path(PROJ,"coloc"); dir.create(OUTD, showWarnings=FALSE, recursive=TRUE); setwd(OUTD)

SUSIE_WIN <- 200000L    # SuSiE/LD 창 (lead eQTL ±)  — LD 계산 tractable
NARROW    <- 100000L    # 좁은창 abf 민감도 (lead eQTL ±)
MAX_SNP   <- 2500L      # LD 행렬 상한(너무 크면 thinning)
set.seed(20260819)

## ---- LD 로컬 인프라 (Step8b와 동일) ------------------------------------- ##
LOCAL_BFILE <- "PATH/TO/ldref_ascii/EUR"
{ ad<-"PATH/TO/ldref_ascii"; dir.create(ad,showWarnings=FALSE,recursive=TRUE)
  for(e in c("bed","bim","fam")){d<-file.path(ad,paste0("EUR.",e));s<-file.path(PROJ,"ld_ref",paste0("EUR.",e))
    if(!file.exists(d)&&file.exists(s)) file.copy(s,d)} }
if (!requireNamespace("genetics.binaRies",quietly=TRUE)) remotes::install_github("MRCIEU/genetics.binaRies")
LOCAL_PLINK <- genetics.binaRies::get_plink_binary()
stopifnot(file.exists(paste0(LOCAL_BFILE,".bed")))

## ---- hit 정의 (Step8c와 동일 좌표/파일) --------------------------------- ##
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
S_CASE <- c(COA=0.044, AOA=0.081, MDD=0.33, ANX=0.12)

STD_COLS <- c("molecular_trait_id","chromosome","position","ref","alt","variant",
              "ma_samples","maf","pvalue","beta","se","type","ac","an","r2",
              "molecular_trait_object_id","gene_id","median_tpm","rsid")
get_schema <- function(f){ h <- try(names(fread(f,nrows=0,showProgress=FALSE)),silent=TRUE)
  if (inherits(h,"try-error")||!length(h)||h[1]=="V1") return(STD_COLS); h }

## ---- Rsamtools tabix (64bit) — 유전자 cis 창 → gene eQTL data.table ------ ##
read_gene_eqtl <- function(h, pad){
  f <- file.path(RAW, paste0(h$qtd,".all.tsv.gz"))
  stopifnot(file.exists(f), file.exists(paste0(f,".tbi")))
  schema <- get_schema(f)
  seqs <- try(Rsamtools::seqnamesTabix(f), silent=TRUE); cc <- as.character(h$chr)
  if (!inherits(seqs,"try-error") && !(cc %in% seqs) && paste0("chr",cc)%in%seqs) cc<-paste0("chr",cc)
  gr <- GenomicRanges::GRanges(cc, IRanges::IRanges(h$s-pad, h$e+pad))
  ln <- try(Rsamtools::scanTabix(Rsamtools::TabixFile(f), param=gr)[[1]], silent=TRUE)
  if (inherits(ln,"try-error")||!length(ln)) return(NULL)
  dt <- fread(text=paste(ln,collapse="\n"), header=FALSE, sep="\t", showProgress=FALSE)
  setnames(dt, schema[seq_len(ncol(dt))])
  dt <- dt[startsWith(as.character(gene_id), h$ensg)]                 # 대상 유전자만
  dt <- dt[!is.na(rsid) & rsid!="" & is.finite(beta) & is.finite(se) & se>0 &
           is.finite(maf) & maf>0 & maf<0.5 &
           ref %in% c("A","C","G","T") & alt %in% c("A","C","G","T")]
  dt <- dt[!duplicated(rsid)]
  if (!nrow(dt)) return(NULL)
  dt[, .(rsid, pos=position, EA=alt, NEA=ref, b=beta, se=se, maf=maf, p=pvalue)]  # eQTL EA=alt
}

## ---- 질환 harmonized 로드 (region 서브셋) ------------------------------- ##
load_disease <- function(tr, snps){
  d <- as.data.table(readRDS(file.path(HARM, paste0(tr,".rds"))))
  d <- d[SNP %in% snps & is.finite(BETA) & is.finite(SE) & SE>0 & is.finite(EAF) & EAF>0 & EAF<1]
  d <- d[!duplicated(SNP)]
  if (!nrow(d)) return(NULL)
  d[, p := 2*pnorm(-abs(BETA/SE))]
  d[, .(rsid=SNP, EA=EA, NEA=NEA, b=BETA, se=SE, eaf=EAF, N=N, p=p)]
}

SIG_P <- 1e-5   # SuSiE 게이트: 창 안에서 양쪽 형질 최소 p 가 이보다 작아야 진행

## ---- LD(부호 있는 r) 행렬 + 정렬 keyed by rsid -------------------------- ##
get_ld <- function(rsids){
  rsids <- unique(rsids)
  if (length(rsids) > MAX_SNP) rsids <- rsids[seq_len(MAX_SNP)]
  ld <- try(ieugwasr::ld_matrix_local(rsids, bfile=LOCAL_BFILE, plink_bin=LOCAL_PLINK,
             with_alleles=TRUE), silent=TRUE)
  if (inherits(ld,"try-error") || is.null(ld) || nrow(ld)<10) return(NULL)
  ld
}

## ---- 한 형질을 LD의 A1 기준으로 정렬(부호 맞춤) ------------------------- ##
##  ld colnames = "rsid_A1_A2".  형질 EA==A1 → 유지 / EA==A2 → beta 부호반전.
align_trait <- function(df, ldinfo){
  m <- merge(df, ldinfo, by="rsid")
  same <- m$EA==m$A1 & m$NEA==m$A2
  flip <- m$EA==m$A2 & m$NEA==m$A1
  m <- m[same | flip]
  fl <- m$EA==m$A2 & m$NEA==m$A1                     # 서브셋에서 재계산
  m[fl, b := -b]                                     # A2가 effect면 부호 반전
  m
}

## ---- SuSiE-coloc 한 쌍(gene eQTL vs 질환) ------------------------------- ##
run_susie_pair <- function(eq, dis, h, tr){
  common <- intersect(eq$rsid, dis$rsid)
  if (length(common) < 50) return(list(row=NULL, note=sprintf("공통SNP %d<50", length(common))))
  ## lead eQTL 중심 SUSIE_WIN 창
  lead_pos <- eq[which.min(p), pos]
  win <- eq[abs(pos-lead_pos)<=SUSIE_WIN, rsid]
  eq2 <- eq[rsid %in% win & rsid %in% common]; di2 <- dis[rsid %in% eq2$rsid]
  if (nrow(eq2) < 50) return(list(row=NULL, note="창내 SNP<50"))

  ## 게이트: 창 안에서 양쪽 형질에 실제 신호(min p<SIG_P)가 있어야 SuSiE 진행.
  ##  (신호 없는 데이터에 susie_rss 돌리면 수치오류 + colocalize 자체 불가)
  eqP <- min(eq2$p, na.rm=TRUE); diP <- min(di2$p, na.rm=TRUE)
  if (eqP >= SIG_P) return(list(row=data.table(GENE=h$gene,TISSUE=h$tissue,OUTCOME=tr,
        nSNP=nrow(eq2), nCS_eqtl=NA_integer_, nCS_dis=NA_integer_, maxPP_H4=NA_real_, best_H3=NA_real_),
        note=sprintf("eQTL 신호 미약(min p=%.1e)", eqP)))
  if (diP >= SIG_P) return(list(row=data.table(GENE=h$gene,TISSUE=h$tissue,OUTCOME=tr,
        nSNP=nrow(eq2), nCS_eqtl=NA_integer_, nCS_dis=NA_integer_, maxPP_H4=NA_real_, best_H3=NA_real_),
        note=sprintf("질환 신호 미약(min p=%.1e, 검정력)", diP)))

  ld <- get_ld(eq2$rsid)
  if (is.null(ld)) return(list(row=NULL, note="LD 실패"))
  ## 주의: data.table(key=..) 의 'key'는 예약 인자 → 컬럼명 'vid' 사용
  info <- data.table(vid=colnames(ld))
  info[, c("rsid","A1","A2") := tstrsplit(vid,"_",fixed=TRUE)]
  ldinfo <- info[, .(rsid,A1,A2)]

  ea <- align_trait(eq2, ldinfo); da <- align_trait(di2, ldinfo)
  keep <- Reduce(intersect, list(ea$rsid, da$rsid, ldinfo$rsid))
  if (length(keep) < 50) return(list(row=NULL, note=sprintf("정렬후 SNP %d<50", length(keep))))
  ## LD 서브셋 + 순서 통일
  idx <- match(keep, info$rsid)
  L <- ld[idx, idx, drop=FALSE]; rownames(L)<-colnames(L)<-keep
  L[is.na(L)] <- 0; diag(L) <- 1
  ea <- ea[match(keep, rsid)]; da <- da[match(keep, rsid)]

  D1 <- list(beta=ea$b, varbeta=ea$se^2, snp=keep, MAF=ea$maf, N=h$n_eqtl,
             type="quant", sdY=1, LD=L)
  D2 <- list(beta=da$b, varbeta=da$se^2, snp=keep, MAF=pmin(da$eaf,1-da$eaf),
             N=round(median(da$N)), s=unname(S_CASE[tr]), type="cc", LD=L)

  s1 <- try(suppressWarnings(runsusie(D1)), silent=TRUE)
  s2 <- try(suppressWarnings(runsusie(D2)), silent=TRUE)
  if (inherits(s1,"try-error")) return(list(row=NULL, note="eQTL SuSiE 실패"))
  if (inherits(s2,"try-error")) return(list(row=NULL, note="질환 SuSiE 실패"))
  ncs1 <- length(s1$sets$cs); ncs2 <- length(s2$sets$cs)
  if (is.null(ncs1)||ncs1==0) return(list(row=data.table(GENE=h$gene,TISSUE=h$tissue,OUTCOME=tr,
        nSNP=length(keep), nCS_eqtl=0L, nCS_dis=as.integer(ncs2), maxPP_H4=NA_real_),
        note="eQTL 신뢰집합 0"))
  if (is.null(ncs2)||ncs2==0) return(list(row=data.table(GENE=h$gene,TISSUE=h$tissue,OUTCOME=tr,
        nSNP=length(keep), nCS_eqtl=as.integer(ncs1), nCS_dis=0L, maxPP_H4=NA_real_),
        note="질환 신뢰집합 0(검정력)"))

  cs <- try(suppressWarnings(coloc.susie(s1, s2)), silent=TRUE)
  if (inherits(cs,"try-error") || is.null(cs$summary) || !nrow(cs$summary))
    return(list(row=data.table(GENE=h$gene,TISSUE=h$tissue,OUTCOME=tr,nSNP=length(keep),
        nCS_eqtl=as.integer(ncs1),nCS_dis=as.integer(ncs2),maxPP_H4=NA_real_), note="coloc.susie 무결과"))
  sm <- as.data.table(cs$summary)
  best <- sm[which.max(PP.H4.abf)]
  list(row=data.table(GENE=h$gene, TISSUE=h$tissue, OUTCOME=tr, nSNP=length(keep),
        nCS_eqtl=as.integer(ncs1), nCS_dis=as.integer(ncs2),
        maxPP_H4=round(best$PP.H4.abf,3),
        best_H3=round(best$PP.H3.abf,3)),
       note=sprintf("신호쌍 %d개", nrow(sm)))
}

## ---- 좁은창 coloc.abf 민감도 -------------------------------------------- ##
run_abf_narrow <- function(eq, dis, h, tr){
  lead_pos <- eq[which.min(p), pos]
  e <- eq[abs(pos-lead_pos)<=NARROW]; d <- dis[rsid %in% e$rsid]
  m <- merge(e[,.(rsid,b_e=b,v_e=se^2,maf_e=maf)],
             d[,.(rsid,b_d=b,v_d=se^2,eaf,N)], by="rsid")
  m[, maf_d:=pmin(eaf,1-eaf)]; m <- m[maf_d>0 & maf_d<0.5]
  if (nrow(m) < 20) return(NULL)
  r <- try(coloc.abf(
    list(beta=m$b_e,varbeta=m$v_e,snp=m$rsid,type="quant",sdY=1,MAF=m$maf_e,N=h$n_eqtl),
    list(beta=m$b_d,varbeta=m$v_d,snp=m$rsid,type="cc",MAF=m$maf_d,N=round(median(m$N)),s=unname(S_CASE[tr]))),
    silent=TRUE)
  if (inherits(r,"try-error")) return(NULL)
  s <- as.list(r$summary)
  data.table(GENE=h$gene,TISSUE=h$tissue,OUTCOME=tr,nSNP=s$nsnps,
             PP_H3=round(s$PP.H3.abf,3), PP_H4=round(s$PP.H4.abf,3))
}

## ---- 실행 --------------------------------------------------------------- ##
message("=== Step 8d: SuSiE-coloc + 좁은창 abf 민감도 ===")
SUS <- list(); ABF <- list()
for (h in HITS){
  message(sprintf("\n-- %s [%s]", h$gene, h$tissue))
  eq <- read_gene_eqtl(h, pad=SUSIE_WIN + 50000L)   # 넉넉히 읽고 창 잘라씀
  if (is.null(eq) || !nrow(eq)) { message("   eQTL 없음 → skip"); next }
  message(sprintf("   eQTL SNP %d (lead p=%.1e)", nrow(eq), min(eq$p)))
  for (tr in h$dis){
    dis <- load_disease(tr, eq$rsid)
    if (is.null(dis)) { message(sprintf("   [%s] 질환 매칭 0 → skip", tr)); next }
    ## 좁은창 abf
    a <- run_abf_narrow(eq, dis, h, tr); if (!is.null(a)) ABF[[length(ABF)+1]] <- a
    ## SuSiE-coloc (한 쌍 실패해도 전체는 계속 — 에러는 note로 기록)
    rs <- tryCatch(run_susie_pair(eq, dis, h, tr),
                   error=function(e) list(row=NULL, note=paste("ERROR:", conditionMessage(e))))
    if (!is.null(rs$row)) { SUS[[length(SUS)+1]] <- rs$row
      message(sprintf("   ★ SuSiE %s→%s : maxPP.H4=%s (eQTL CS=%s, 질환 CS=%s) [%s]",
        h$gene, tr, ifelse(is.na(rs$row$maxPP_H4),"NA",sprintf("%.3f",rs$row$maxPP_H4)),
        rs$row$nCS_eqtl, rs$row$nCS_dis, rs$note))
    } else message(sprintf("   · SuSiE %s→%s : %s", h$gene, tr, rs$note))
  }
}

## ---- 저장/요약 ---------------------------------------------------------- ##
if (length(ABF)){ ABFdt<-rbindlist(ABF,fill=TRUE); fwrite(ABFdt,"coloc_abf_narrow.csv")
  message("\n=== 좁은창 coloc.abf (lead ±100kb) ===\n"); print(ABFdt[order(-PP_H4)]) }

if (length(SUS)){
  SUSdt <- rbindlist(SUS, fill=TRUE); setorder(SUSdt, -maxPP_H4, na.last=TRUE)
  fwrite(SUSdt, "coloc_susie_summary.csv")
  message("\n=== SuSiE-coloc 요약 (maxPP.H4 내림차순) ===\n"); print(SUSdt)

  pd <- SUSdt[!is.na(maxPP_H4)]
  if (nrow(pd)){
    pd[, lab:=paste0(GENE," [",substr(TISSUE,1,4),"] → ",OUTCOME)]
    setorder(pd, maxPP_H4); pd[, lab:=factor(lab,levels=lab)]; pd[, strong:=maxPP_H4>=0.8]
    p <- ggplot(pd, aes(maxPP_H4,lab)) +
      geom_col(aes(fill=strong),width=.65,color="grey30") +
      geom_vline(xintercept=0.8,linetype=2,color="#B2182B") +
      geom_text(aes(label=sprintf("%.2f",maxPP_H4)),hjust=-0.15,size=3) +
      scale_fill_manual(values=c(`TRUE`="#B2182B",`FALSE`="grey75"),guide="none") +
      scale_x_continuous(limits=c(0,1.08),breaks=seq(0,1,.2)) +
      labs(title="SuSiE colocalization (multiple causal variants)",
           subtitle="max PP.H4 over signal pairs; dashed=0.8", x="max PP.H4", y=NULL) +
      theme_minimal(base_size=10)
    ggsave("fig_susie_pph4.png", p, width=7.5, height=4.5, dpi=300, bg="white")
    message("  fig_susie_pph4.png 저장")
  }
  strong <- SUSdt[!is.na(maxPP_H4) & maxPP_H4>=0.8]
  message("\n★ SuSiE 공유변이(PP.H4≥0.8): ",
    if(nrow(strong)) paste(sprintf("%s[%s]→%s(%.2f)", strong$GENE, substr(strong$TISSUE,1,4),
      strong$OUTCOME, strong$maxPP_H4), collapse=" · ") else "없음")
} else message("\n[SuSiE 결과 없음]")

message("\n================= Step 8d 완료 =================")
message("산출: coloc/coloc_susie_summary.csv, coloc/coloc_abf_narrow.csv, coloc/fig_susie_pph4.png")
message("해석: maxPP.H4≥0.8 = 다중신호 상황에서도 공유변이 존재(확증).")
message("      질환 CS=0 = 그 구간 질환신호가 약해 판정불가(특히 MDD/ANX polygenic 검정력 한계).")
