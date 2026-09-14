###############################################################################
## Step 6 : 논문 그림 생성 (Nature급 composite figures)
## Project: 01_GWAS_Asthma_MDD
##
## 결과 CSV/RDS를 읽어 본문 4개 복합그림 자동 생성. 색,폰트,라벨은 여기서 조정.
##   Fig 1  Overview & heritability  : (a)파이프라인 (b)SNP-h2 (c)표본중복 intercept
##   Fig 2  Genetic correlation      : (a)rg 히트맵 (b)핵심 rg forest (c)LAVA 국소방향
##   Fig 3  Pleiotropy & immune      : (a)PLACO Manhattan (b)면역유전자 x 페어 dot matrix
##   Fig 4  Causal (MR)              : (a)forest 3패널 (b)MDD→AOA scatter (c)funnel
##
## 입력: ldsc/{rg_results,h2_results,intercept_matrix}.csv, lava/lava_bivar.csv,
##       placo/placo_lead_annotated.csv, mr/mr_summary_*.csv, mr/MR_reverse_MDD_AOA.rds
## 출력: fig/ (Fig1-4 png+pdf) , tab/ (Table1-3 csv)
## 실행: source("step6_figures.R")
###############################################################################

## ------------------------- 0. 패키지 ------------------------------------- ##
.need <- c("data.table","ggplot2","scales","patchwork","ggrepel")
.miss <- .need[!sapply(.need, requireNamespace, quietly=TRUE)]
if (length(.miss)) install.packages(.miss, repos="https://cloud.r-project.org")
suppressPackageStartupMessages({ library(data.table); library(ggplot2)
  library(scales); library(patchwork); library(ggrepel) })

## ------------------------- 1. 경로 & 공통 -------------------------------- ##
PROJ <- "PATH/TO/PROJECT"
LDSC <- file.path(PROJ,"ldsc"); PLACO <- file.path(PROJ,"placo")
LAVA <- file.path(PROJ,"lava"); MRD <- file.path(PROJ,"mr")
FIG  <- file.path(PROJ,"fig"); TAB <- file.path(PROJ,"tab")
dir.create(FIG, showWarnings=FALSE); dir.create(TAB, showWarnings=FALSE)

TRAITS <- c("COA","AOA","MDD","ANX","BIP","SCZ","FEV1","FVC","LUNG")
LAB    <- c(COA="COA",AOA="AOA",MDD="MDD",ANX="ANX",BIP="BIP",SCZ="SCZ",
            FEV1="FEV1",FVC="FVC",LUNG="FEV1/FVC")
## 도메인 팔레트 (CVD-safe, 절제된 톤)
DOMdt <- data.table(TRAIT=TRAITS,
  DOMAIN=c("Asthma","Asthma","Psychiatric","Psychiatric","Psychiatric","Psychiatric",
           "Lung","Lung","Lung"))
DOMCOL <- c(Asthma="#3B7DA8", Psychiatric="#B24745", Lung="#6E9B7A")
INK <- "#222222"; MUTE <- "#6b6b6b"

stars <- function(p) ifelse(is.na(p),"", ifelse(p<1e-3,"***", ifelse(p<1e-2,"**", ifelse(p<0.05,"*",""))))
theme_pub <- theme_minimal(base_size=10) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_line(color="grey92", linewidth=0.3),
        plot.title=element_text(face="bold", size=10.5, color=INK),
        plot.subtitle=element_text(size=8.5, color=MUTE),
        plot.tag=element_text(face="bold", size=13),
        axis.title=element_text(size=9, color=INK),
        axis.text=element_text(color=INK),
        legend.text=element_text(size=8), legend.title=element_text(size=8.5))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name,".png")), p, width=w, height=h, dpi=340, bg="white")
  ggsave(file.path(FIG, paste0(name,".pdf")), p, width=w, height=h, device=cairo_pdf)
  message("  저장: ", name, " (", w, "x", h, ")")
}
rd <- function(f) if (file.exists(f)) fread(f) else NULL

###############################################################################
## Fig 1 : Overview & heritability
###############################################################################
message("\n=== Fig 1 ===")
## (a) 분석 파이프라인 스키마 (boxes + arrows)
steps <- data.table(
  x = 1:5, lab = c("Traits","LDSC","LAVA","PLACO","MR"),
  sub = c("Asthma (COA/AOA)\nMDD, ANX, BIP, SCZ\nFEV1, FVC, ratio",
          "Global r_g\n+ sample overlap","Local r_g\n(882 asthma loci)",
          "Pleiotropic SNPs\n(shared association)","Bidirectional\ncausality"))
p1a <- ggplot(steps, aes(x, 1)) +
  geom_segment(data=steps[x<5], aes(x=x+0.34, xend=x+0.66, y=1, yend=1),
               arrow=arrow(length=unit(0.14,"cm"), type="closed"), color=MUTE, linewidth=0.5) +
  geom_tile(width=0.66, height=0.9, fill="white", color="grey45", linewidth=0.5) +
  geom_text(aes(label=lab), y=1.28, fontface="bold", size=3.1, color=INK) +
  geom_text(aes(label=sub), y=0.86, size=2.15, color=MUTE, lineheight=0.92) +
  scale_x_continuous(limits=c(0.5,5.5)) + ylim(0.55,1.5) +
  labs(title="Study design and analytical pipeline") +
  theme_void() + theme(plot.title=element_text(face="bold", size=10.5, color=INK, hjust=0))

## (b) SNP-h2 막대 (CI, 도메인 색)
h2 <- merge(rd(file.path(LDSC,"h2_results.csv")), DOMdt, by="TRAIT", sort=FALSE)
h2[, TRAIT := factor(TRAIT, levels=rev(TRAITS))]
p1b <- ggplot(h2, aes(h2, TRAIT, fill=DOMAIN)) +
  geom_col(width=0.68) +
  geom_errorbarh(aes(xmin=pmax(0,h2-1.96*h2_SE), xmax=h2+1.96*h2_SE), height=0.28, color="grey30", linewidth=0.4) +
  geom_text(aes(label=sprintf("Z=%.0f", h2_Z)), hjust=-0.15, size=2.5, color=MUTE) +
  scale_fill_manual(values=DOMCOL, name=NULL) +
  scale_y_discrete(labels=LAB[rev(TRAITS)]) +
  scale_x_continuous(expand=expansion(mult=c(0,0.30))) +
  labs(title="SNP heritability", subtitle="liability scale (binary), observed (lung)",
       x=expression(italic(h)^2~"(95% CI)"), y=NULL) +
  theme_pub + theme(legend.position="bottom", legend.margin=margin(t=-6))

## (c) 표본중복 = 교차 intercept 히트맵 (비대각)
im <- rd(file.path(LDSC,"intercept_matrix.csv"))
imM <- as.matrix(im[, ..TRAITS]); rownames(imM) <- im$TRAIT
imdt <- as.data.table(as.table(imM)); setnames(imdt, c("T1","T2","val"))
imdt <- imdt[T1 %in% TRAITS & T2 %in% TRAITS]
imdt[T1==T2, val := NA]
imdt[, `:=`(T1=factor(T1,levels=TRAITS), T2=factor(T2,levels=rev(TRAITS)))]
p1c <- ggplot(imdt, aes(T1,T2,fill=val)) +
  geom_tile(color="white", linewidth=0.5) +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0,
                       limits=c(-0.25,0.25), oob=squish, na.value="grey93",
                       name="cross\nintercept") +
  scale_x_discrete(labels=LAB[TRAITS], position="top") + scale_y_discrete(labels=LAB[rev(TRAITS)]) +
  coord_equal() +
  labs(title="Sample overlap", subtitle="LDSC cross-trait intercept (0 = no overlap)") +
  theme_pub + theme(axis.title=element_blank(), axis.text.x.top=element_text(angle=45,hjust=0),
                    panel.grid=element_blank())

fig1 <- p1a / (p1b | p1c) + plot_layout(heights=c(0.5,1)) +
  plot_annotation(tag_levels="a")
save_fig(fig1, "Fig1_overview", 10, 8.2)

###############################################################################
## Fig 2 : Genetic correlation (global + local)
###############################################################################
message("\n=== Fig 2 ===")
rg <- rd(file.path(LDSC,"rg_results.csv"))
## (a) rg 히트맵
g <- CJ(T1=TRAITS, T2=TRAITS, sorted=FALSE)
kf <- function(a,b) paste(pmin(a,b), pmax(a,b))
rg[, k := kf(T1,T2)]; g[, k := kf(T1,T2)]
g <- merge(g, rg[,.(k,rg,P_FDR)], by="k", all.x=TRUE)
g[T1==T2, `:=`(rg=NA, P_FDR=NA)]
g[, `:=`(T1=factor(T1,levels=TRAITS), T2=factor(T2,levels=rev(TRAITS)), star=stars(P_FDR))]
g[, lab := ifelse(is.na(rg),"",sprintf("%.2f",rg))]
p2a <- ggplot(g, aes(T1,T2,fill=rg)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(aes(label=lab), size=2.4, na.rm=TRUE) +
  geom_text(aes(label=star), nudge_y=-0.27, size=2.3, na.rm=TRUE) +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0,
                       limits=c(-0.4,0.4), oob=squish, na.value="#ededed",
                       name=expression(r[g]), breaks=c(-0.3,0,0.3)) +
  scale_x_discrete(labels=LAB[TRAITS], position="top") + scale_y_discrete(labels=LAB[rev(TRAITS)]) +
  geom_vline(xintercept=c(2.5,6.5), color="grey45", linewidth=0.4) +
  geom_hline(yintercept=c(3.5,7.5), color="grey45", linewidth=0.4) +
  coord_equal() +
  labs(title="Global genetic correlation", subtitle="* FDR<0.05  ** <0.01  *** <0.001  (color capped at |rg|=0.4)") +
  theme_pub + theme(axis.title=element_blank(), axis.text.x.top=element_text(angle=45,hjust=0),
                    panel.grid=element_blank())

## (b) 핵심 rg 페어 forest
key_rg <- rg[(T1%in%c("COA","AOA") & T2%in%c("MDD","ANX","FEV1","FVC","LUNG")) |
             (T1%in%c("MDD","ANX") & T2%in%c("FEV1","FVC","LUNG"))]
key_rg[, `:=`(lo=rg-1.96*SE, hi=rg+1.96*SE,
              pairlab=paste0(LAB[T1]," - ",LAB[T2]),
              grp=fifelse(T2%in%c("FEV1","FVC","LUNG") | T1%in%c("FEV1","FVC","LUNG"),
                          "with lung function","asthma x psychiatric"))]
setorder(key_rg, grp, rg)
key_rg[, pairlab := factor(pairlab, levels=pairlab)]
key_rg[, sig := P_FDR<0.05]
p2b <- ggplot(key_rg, aes(rg, pairlab)) +
  geom_vline(xintercept=0, linetype=2, color="grey55", linewidth=0.4) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.25, color="grey35", linewidth=0.5) +
  geom_point(aes(fill=sig), shape=21, size=2.7, color="grey20") +
  scale_fill_manual(values=c(`TRUE`="#B2182B",`FALSE`="white"), guide="none") +
  facet_grid(grp~., scales="free_y", space="free_y", switch="y") +
  labs(title="Key genetic correlations", subtitle="filled = FDR<0.05",
       x=expression(r[g]~"(95% CI)"), y=NULL) +
  theme_pub + theme(strip.placement="outside", strip.text.y.left=element_text(angle=0, face="bold", size=8),
                    panel.spacing=unit(0.4,"lines"))

## (c) LAVA 국소 방향: 천식 x 정신 페어에서 유의(p<0.05) 국소상관의 +/- locus 수
lv <- rd(file.path(LAVA,"lava_bivar.csv"))
lv <- lv[!is.na(rho) & !is.na(p)]
lv[, pair := ifelse(paste(phen1,phen2,sep="_") %in% c("COA_MDD","COA_ANX","AOA_MDD","AOA_ANX"),
                    paste(phen1,phen2,sep="_"),
                    ifelse(paste(phen2,phen1,sep="_") %in% c("MDD_COA","ANX_COA","MDD_AOA","ANX_AOA"),
                           paste(phen2,phen1,sep="_"), pair))]
psychpairs <- c("COA_MDD","COA_ANX","AOA_MDD","AOA_ANX")
lvp <- lv[pair %in% psychpairs & p<0.05]
lvc <- lvp[, .(pos=sum(rho>0), neg=sum(rho<0)), by=pair]
lvm <- melt(lvc, id.vars="pair", variable.name="dir", value.name="n")
lvm[dir=="neg", n := -n]
lvm[, pairlab := factor(gsub("_"," x ",pair), levels=gsub("_"," x ",psychpairs))]
p2c <- ggplot(lvm, aes(n, pairlab, fill=dir)) +
  geom_col(width=0.6) +
  geom_vline(xintercept=0, color="grey40", linewidth=0.4) +
  scale_fill_manual(values=c(pos="#B2182B", neg="#2166AC"),
                    labels=c(pos="positive local rg", neg="negative local rg"), name=NULL) +
  scale_x_continuous(labels=abs) +
  labs(title="Local genetic correlation direction",
       subtitle="significant local rg loci (p<0.05); COA mixes directions",
       x="number of loci", y=NULL) +
  theme_pub + theme(legend.position="top")

fig2 <- (p2a | (p2b / p2c)) + plot_layout(widths=c(1.05,1)) + plot_annotation(tag_levels="a")
save_fig(fig2, "Fig2_genetic_correlation", 12, 7.2)

###############################################################################
## Fig 3 : Pleiotropy & immune convergence
###############################################################################
message("\n=== Fig 3 ===")
pl <- rd(file.path(PLACO,"placo_lead_annotated.csv"))
psychP <- c("COA_MDD","COA_ANX","AOA_MDD","AOA_ANX")
pl <- pl[pair %in% psychP]
pl[, `:=`(CHR=as.integer(CHR), POS=as.numeric(POS), P_PLACO=as.numeric(P_PLACO))]
pl <- pl[is.finite(CHR)&is.finite(POS)&is.finite(P_PLACO)]
CHRLEN <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,
            141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,
            81195210,78077248,59128983,63025520,48129895,51304566)
off <- c(0, cumsum(as.numeric(CHRLEN))[-22]); names(off) <- 1:22
axisdf <- data.table(CHR=1:22, center=off+CHRLEN/2)
pl[, xpos := off[as.character(CHR)]+POS]; pl[, mlogp := -log10(P_PLACO)]
pl[, pairF := factor(pair, levels=psychP)]
IMMUNE <- "^IL[0-9]|IL[0-9]+R|TLR[0-9]|TYK2|^HLA|SMAD3|GATA3|RUNX[0-9]|IKZF|STAT6|P2RX7|TSLP|ORMDL3|GSDMB|CLEC16|IL1RL1|IL18R|IL2RA|IL21R|IL4R|IL33|IL13|IL6R|RORA|CTLA4|FADS|NFKB"
pl[, glist := strsplit(gsub("\\.\\.$","",gene), ",")]
pick <- function(gs){gs<-trimws(gs); m<-gs[grepl(IMMUNE,gs)]; if(length(m)) m[1] else NA_character_}
pl[, immgene := sapply(glist, pick)]
PAIRLAB <- c(COA_MDD="COA x MDD", COA_ANX="COA x ANX", AOA_MDD="AOA x MDD", AOA_ANX="AOA x ANX")
pal4 <- c(COA_MDD="#2f5f9e", COA_ANX="#7DA9D6", AOA_MDD="#B2182B", AOA_ANX="#E58267")
shp4 <- c(COA_MDD=16, COA_ANX=17, AOA_MDD=15, AOA_ANX=18)
p3a <- ggplot(pl, aes(xpos, mlogp)) +
  geom_hline(yintercept=-log10(5e-8), linetype=2, color="grey55", linewidth=0.4) +
  geom_point(data=pl[LAVA==TRUE], shape=5, size=3.4, stroke=0.8, color="grey15") +
  geom_point(aes(color=pairF, shape=pairF), size=2.1) +
  geom_text_repel(data=pl[!is.na(immgene)], aes(label=immgene, color=pairF),
                  size=2.4, fontface="italic", min.segment.length=0, max.overlaps=18, seed=1, show.legend=FALSE) +
  scale_color_manual(values=pal4, labels=PAIRLAB, name=NULL) +
  scale_shape_manual(values=shp4, labels=PAIRLAB, name=NULL) +
  scale_x_continuous(breaks=axisdf$center, labels=axisdf$CHR, expand=c(0.01,0)) +
  labs(title="Shared pleiotropic loci (PLACO)",
       subtitle="Open diamond = also significant in LAVA; italic = immune/allergy genes",
       x="Chromosome", y=expression(-log[10]~italic(P)[PLACO])) +
  theme_pub + theme(legend.position="top", panel.grid.major.x=element_blank())

## (b) 면역유전자 x 페어 dot matrix (반복 출현 면역유전자)
imm <- pl[!is.na(immgene), .(pair, immgene, mlogp)]
topg <- imm[, .(n=uniqueN(pair), maxp=max(mlogp)), by=immgene][order(-n,-maxp)][n>=1][1:min(14,.N)]
immM <- imm[immgene %in% topg$immgene]
immM[, immgene := factor(immgene, levels=rev(topg$immgene))]
immM[, pairF := factor(pair, levels=psychP, labels=PAIRLAB)]
p3b <- ggplot(immM, aes(pairF, immgene)) +
  geom_point(aes(size=mlogp, color=pairF)) +
  scale_color_manual(values=setNames(pal4,PAIRLAB), guide="none") +
  scale_size_continuous(range=c(2,5.5), name=expression(-log[10]~italic(P))) +
  labs(title="Immune-gene convergence", subtitle="recurrent immune/allergy genes across pairs",
       x=NULL, y=NULL) +
  theme_pub + theme(axis.text.x=element_text(angle=30,hjust=1),
                    axis.text.y=element_text(face="italic"), panel.grid.major=element_line(color="grey93"))

fig3 <- (p3a | p3b) + plot_layout(widths=c(2,1)) + plot_annotation(tag_levels="a")
save_fig(fig3, "Fig3_pleiotropy", 13, 5.4)

###############################################################################
## Fig 4 : Causal (MR)
###############################################################################
message("\n=== Fig 4 ===")
## (a) forest 3패널 (스케일별)
forest <- function(dt, eff, lo, hi, fdr, null, xlab, title, logx=FALSE) {
  d <- copy(dt); setnames(d, c(eff,lo,hi,fdr), c("eff","lo","hi","fdr"))
  d[, label := paste0(LAB[EXPOSURE]," \u2192 ",LAB[OUTCOME])]
  d[, label := factor(label, levels=rev(label))]
  d[, sig := !is.na(fdr) & fdr<0.05]
  pp <- ggplot(d, aes(eff,label)) +
    geom_vline(xintercept=null, linetype=2, color="grey55", linewidth=0.4) +
    geom_errorbarh(aes(xmin=lo,xmax=hi), height=0.2, color="grey35", linewidth=0.45) +
    geom_point(aes(fill=sig), shape=21, size=2.6, color="grey20") +
    scale_fill_manual(values=c(`TRUE`="#B2182B",`FALSE`="white"), guide="none") +
    labs(title=title, x=xlab, y=NULL) + theme_pub
  if (logx) pp <- pp + scale_x_log10()
  pp
}
mrOR <- rbindlist(list(rd(file.path(MRD,"mr_summary_forward.csv")), rd(file.path(MRD,"mr_summary_reverse.csv"))), fill=TRUE)
pA <- forest(mrOR,"IVW_OR","IVW_LCI","IVW_UCI","IVW_P_FDR",1,"OR (95% CI)","Asthma \u2194 Depression / Anxiety",TRUE)
mrB <- rbindlist(list(rd(file.path(MRD,"mr_summary_lung_asthma2lung.csv")), rd(file.path(MRD,"mr_summary_lung_mental2lung.csv"))), fill=TRUE)
pB <- forest(mrB,"IVW_B","IVW_B_LCI","IVW_B_UCI","IVW_P_FDR",0,"Lung function change (SD)","\u2192 Lung function",FALSE)
mrC <- rbindlist(list(rd(file.path(MRD,"mr_summary_lung_lung2asthma.csv")), rd(file.path(MRD,"mr_summary_lung_lung2mental.csv"))), fill=TRUE)
pC <- forest(mrC,"IVW_OR","OR_LCI","OR_UCI","IVW_P_FDR",1,"OR (95% CI) per SD","Lung function \u2192",TRUE)
p4a <- pA / pB / pC + plot_layout(heights=c(0.8,1.1,1.1))

## (b) MDD→AOA scatter + 방법별 기울기
o <- readRDS(file.path(MRD,"MR_reverse_MDD_AOA.rds"))
dat <- as.data.table(o$dat); res <- as.data.table(o$res)
sl_ivw <- res[method=="Inverse variance weighted", b]
sl_wm  <- res[method=="Weighted median", b]
egg_b  <- res[method=="MR Egger", b]; egg_i <- if(!is.null(o$plei)) o$plei$egger_intercept else 0
dat[, `:=`(bx=beta.exposure, by=beta.outcome, sx=se.exposure, sy=se.outcome)]
## 노출 부호 양수로 정렬
dat[bx<0, `:=`(bx=-bx, by=-by)]
p4b <- ggplot(dat, aes(bx, by)) +
  geom_hline(yintercept=0, color="grey85", linewidth=0.3) + geom_vline(xintercept=0, color="grey85", linewidth=0.3) +
  geom_errorbar(aes(ymin=by-sy, ymax=by+sy), color="grey80", linewidth=0.25, width=0) +
  geom_errorbarh(aes(xmin=bx-sx, xmax=bx+sx), color="grey80", linewidth=0.25, height=0) +
  geom_point(size=1.3, color="#444444", alpha=0.8) +
  geom_abline(aes(slope=sl_ivw, intercept=0, color="IVW"), linewidth=0.8) +
  geom_abline(aes(slope=sl_wm, intercept=0, color="Weighted median"), linewidth=0.7, linetype=1) +
  geom_abline(aes(slope=egg_b, intercept=egg_i, color="MR-Egger"), linewidth=0.7, linetype=2) +
  scale_color_manual(values=c(IVW="#B2182B",`Weighted median`="#3B7DA8",`MR-Egger`="#6E9B7A"), name=NULL) +
  labs(title="MDD \u2192 AOA: instrument effects",
       subtitle=sprintf("IVW OR=1.19, P=2.5e-4; %d SNPs", nrow(dat)),
       x="SNP effect on MDD (log-OR)", y="SNP effect on AOA (log-OR)") +
  theme_pub + theme(legend.position=c(0.72,0.16), legend.background=element_rect(fill=alpha("white",0.7),color=NA))

## (c) funnel (다면발현 비대칭 점검)
dat[, wald := by/bx][, wse := sy/abs(bx)][, prec := 1/wse]
p4c <- ggplot(dat, aes(wald, prec)) +
  geom_vline(xintercept=sl_ivw, color="#B2182B", linewidth=0.6) +
  geom_point(size=1.3, color="#444444", alpha=0.75) +
  labs(title="MDD \u2192 AOA: funnel",
       subtitle="symmetric = no directional pleiotropy",
       x="Wald ratio (per-SNP causal estimate)", y=expression(1/SE)) +
  theme_pub
p4bc <- p4b / p4c + plot_layout(heights=c(1.25,1))

fig4 <- (p4a | p4bc) + plot_layout(widths=c(1,1)) + plot_annotation(tag_levels="a")
save_fig(fig4, "Fig4_MR", 13, 10.5)

###############################################################################
## Tables 1-3
###############################################################################
message("\n=== Tables ===")
meta <- data.table(TRAIT=TRAITS,
  Phenotype=c("Childhood-onset asthma","Adult-onset asthma","Major depression","Anxiety disorders",
              "Bipolar disorder","Schizophrenia","FEV1","FVC","FEV1/FVC ratio"),
  Source=c("Ferreira 2019","Ferreira 2019","PGC-MDD 2025 (noUKBB)","PGC-ANX 2026",
           "PGC3-BIP 2021 (noUKBB)","PGC3-SCZ w3","Shrine 2019 SpiroMeta","Shrine 2019 SpiroMeta","Shrine 2019 SpiroMeta"),
  Ancestry="European",
  Type=c("binary","binary","binary","binary","binary","binary","continuous","continuous","continuous"))
t1 <- merge(meta, rd(file.path(LDSC,"h2_results.csv"))[,.(TRAIT,h2,h2_SE,h2_Z,intercept)], by="TRAIT", sort=FALSE)
fwrite(t1, file.path(TAB,"Table1_traits.csv")); message("  Table1_traits.csv")
t2 <- rg[focus!="" | (T1%in%c("MDD","ANX")&T2%in%c("FEV1","FVC","LUNG")) |
         (T2%in%c("MDD","ANX")&T1%in%c("FEV1","FVC","LUNG"))][order(P)]
fwrite(t2[,.(T1,T2,rg,SE,P,P_FDR,gcov_intercept,focus)], file.path(TAB,"Table2_rg_focus.csv")); message("  Table2_rg_focus.csv")
mk <- function(dt,scale) if(is.null(dt)||!nrow(dt)) NULL else dt[, .(EXPOSURE,OUTCOME,DIRECTION=TAG,N_IV,F=F_MEAN,
  BETA=IVW_B,BETA_SE=IVW_SE, OR=if("IVW_OR"%in%names(dt)) IVW_OR else NA_real_,
  P=IVW_P,P_FDR=IVW_P_FDR,EGGER_INT_P,Q_P,STEIGER=STEIGER_OK,PRESSO_P,SCALE=scale)]
t3 <- rbindlist(list(mk(mrOR,"OR"), mk(mrB,"beta-SD"), mk(mrC,"OR")), fill=TRUE)
t3 <- unique(t3, by=c("EXPOSURE","OUTCOME","DIRECTION"))
fwrite(t3, file.path(TAB,"Table3_MR.csv")); message("  Table3_MR.csv")

message("\n================= Step 6 완료 =================")
message("그림: ", FIG, "  (Fig1-4 png+pdf)")
message("표:   ", TAB)
