## Manuscript figures — reframed story (robust results)
## 실행: source("Step9_figures.R")  (프로젝트 폴더 결과 CSV들을 읽어 Result/ 에 Fig1-4 저장)
suppressPackageStartupMessages({
  for(p in c("data.table","ggplot2","patchwork","ragg","svglite")) if(!requireNamespace(p,quietly=TRUE)) install.packages(p,repos="https://cloud.r-project.org")
  library(data.table); library(ggplot2); library(patchwork)})
if (packageVersion("patchwork") < "1.3.0") { install.packages("patchwork", repos="https://cloud.r-project.org"); library(patchwork) }  # free(type="label") 필요
B <- "PATH/TO/PROJECT"   # 프로젝트 루트
O <- file.path(B,"Result"); dir.create(O, showWarnings=FALSE, recursive=TRUE)
## 최종 제출용 세트(PNG 300 dpi + TIFF 600 dpi + PDF 벡터 + SVG): Submission/01 JACI/preview_png/Final_v{FIG_VER}
FIG_VER <- 14   # v14: Fig 1A 셀 라벨 반올림(0.055→0.06, 본문·Table II와 일치), Fig 3B 패널명 Childhood-onset. v13: 제출용 TIFF 600 dpi(JACI 선/하프톤 혼합 ≥500 dpi) + 벡터 PDF 추가; 그림 내용은 v12와 동일
## v9 부터 v0/Figures 아래에 저장
FINAL <- file.path(B,"Submission","01 JACI","v0","Figures",paste0("Final_v",FIG_VER)); dir.create(FINAL, showWarnings=FALSE, recursive=TRUE)
ARROW <- "\u2192"   # 논문용 화살표 (→); "->" 대신 사용
save_fig <- function(name, plot, w, h){
  ## Result/ (작업용) + Final_v{n}/ (제출용) 에 저장. ragg 로 유니코드(→) 렌더링 보장
  ggsave(file.path(O,paste0(name,".png")),  plot, width=w, height=h, dpi=300, bg="white", device=ragg::agg_png)
  ggsave(file.path(O,paste0(name,".tiff")), plot, width=w, height=h, dpi=300, bg="white", device=ragg::agg_tiff, compression="lzw")
  ggsave(file.path(FINAL,paste0(name,".png")),  plot, width=w, height=h, dpi=300, bg="white", device=ragg::agg_png)
  ggsave(file.path(FINAL,paste0(name,".tiff")), plot, width=w, height=h, dpi=600, bg="white", device=ragg::agg_tiff, compression="lzw")   # 제출용 600 dpi
  ggsave(file.path(FINAL,paste0(name,".pdf")),  plot, width=w, height=h, bg="white", device=grDevices::cairo_pdf)                       # 벡터(EPS/PDF 허용)
  ggsave(file.path(FINAL,paste0(name,".svg")),  plot, width=w, height=h, bg="white", device=svglite::svglite)
}
RED<-"#B2182B"; BLU<-"#2166AC"; GREY<-"grey70"
## v12: 패널 안 축/쌍 라벨은 약어(PGC 관례), 범례·제목은 단어 (LAB)
ABB <- c(COA="COA", AOA="AOA", MDD="MDD", ANX="ANX", BIP="BIP", SCZ="SCZ", FEV1="FEV1", FVC="FVC", LUNG="FEV1/FVC")
LAB <- c(COA="Childhood-onset asthma", AOA="Adult-onset asthma", MDD="Depression",
         ANX="Anxiety", BIP="Bipolar", SCZ="Schizophrenia",
         FEV1="FEV1", FVC="FVC", LUNG="FEV1/FVC ratio")
## JACI: 그림 내 글꼴 Times New Roman, 형식 TIFF/JPG ≥300 dpi
FF <- "serif"   # Windows 에서 serif = Times New Roman
## JACI 인쇄 폭(2단 = 7.2 in)에 맞춰 실제 크기로 저장 → 축소 시 글자가 작아지지 않음
FW <- 7.2; BS <- 8
th <- theme_minimal(base_size=BS, base_family=FF) +
  theme(panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold",size=BS+1,family=FF),
        plot.tag=element_text(face="bold",size=BS+3,family=FF),
        axis.text=element_text(color="grey20"))

## ================= FIG 1 — Shared genetic architecture (LDSC) ============= ##
rg  <- fread(file.path(B,"ldsc/rg_results.csv"))
ord <- c("COA","AOA","MDD","ANX","BIP","SCZ","FEV1","FVC","LUNG")
## 1a heatmap (symmetric build from pairs)
M <- matrix(NA_real_, 9, 9, dimnames=list(ord,ord)); diag(M)<-1
Pm<- matrix(NA_real_, 9, 9, dimnames=list(ord,ord))
for(i in 1:nrow(rg)){a<-rg$T1[i];b<-rg$T2[i]
  M[a,b]<-M[b,a]<-rg$rg[i]; Pm[a,b]<-Pm[b,a]<-rg$P_FDR[i]}
hm <- as.data.table(as.table(M)); setnames(hm,c("V1","V2","rg"))
hp <- as.data.table(as.table(Pm)); setnames(hp,c("V1","V2","pf"))
hm <- merge(hm,hp,by=c("V1","V2"),all.x=TRUE)
hm[,V1:=factor(V1,levels=ord)][,V2:=factor(V2,levels=rev(ord))]
hm[,star:=ifelse(!is.na(pf)&pf<0.05,"*","")]
r2h <- function(x) sign(x)*floor(abs(x)*100+0.5+1e-9)/100   # v14: 반올림(half away from zero) — sprintf 는 0.055→0.05
hm[,rglab:=ifelse(is.na(rg),"",sprintf("%.2f%s",r2h(rg),star))]
f1a <- ggplot(hm,aes(V1,V2,fill=rg))+
  geom_tile(color="white",linewidth=.6)+
  geom_text(aes(label=rglab),size=BS*0.28,family=FF)+
  scale_fill_gradient2(low=BLU,mid="white",high=RED,midpoint=0,limits=c(-1,1),
     na.value="grey92",name=expression(r[g]))+
  scale_x_discrete(labels=ABB)+scale_y_discrete(labels=ABB)+
  labs(title="Global genetic correlations (LDSC)",x=NULL,y=NULL,
       subtitle="* FDR<0.05")+
  guides(fill=guide_colorbar(title.position="left",title.vjust=0.9,barwidth=unit(90,"pt"),barheight=unit(7,"pt")))+
  th+theme(axis.text.x=element_text(angle=0,hjust=0.5),legend.position="bottom",   # 약어 라벨이라 회전 불필요; 색상 막대는 아래 가로형
           legend.box.spacing=unit(2,"pt"), legend.margin=margin(0,0,0,0),
           panel.grid.major=element_blank())

## 1b forest — asthma x trait (subtype specificity)
fo <- rg[focus=="asthma_x_trait"]
fo[,lci:=rg-1.96*SE][,uci:=rg+1.96*SE]
fo[,trait:=factor(T2,levels=c("MDD","ANX","BIP","SCZ","LUNG","FVC","FEV1"))]
fo[,sub:=factor(T1,levels=c("AOA","COA"),labels=c("Adult-onset","Childhood-onset"))]
fo[,sig:=P_FDR<0.05]
## 성인기 발병 = 세모(24), 소아기 발병 = 네모(22); FDR<0.05 = 채움, 그 외 = 흰색
fo[,fillc:=ifelse(sig, ifelse(sub=="Adult-onset",RED,BLU), "white")]
## 범례 키: 그래프와 같은 T자 막힌 가로 오차막대
draw_key_errorbarh <- function(data, params, size) {
  gp <- grid::gpar(col = data$colour %||% "black", lwd = (data$linewidth %||% 0.5) * .pt, lineend = "butt")
  grid::grobTree(grid::segmentsGrob(0.08, 0.5, 0.92, 0.5, gp = gp),
                 grid::segmentsGrob(0.08, 0.28, 0.08, 0.72, gp = gp),
                 grid::segmentsGrob(0.92, 0.28, 0.92, 0.72, gp = gp))
}
f1b <- ggplot(fo,aes(rg,trait,color=sub,shape=sub,group=sub))+   # group=sub: 오차막대와 점의 dodge 위치 일치
  geom_vline(xintercept=0,linetype=2,color="grey55")+
  geom_errorbarh(aes(xmin=lci,xmax=uci),height=.25,position=position_dodge(.55),key_glyph=draw_key_errorbarh)+
  geom_point(aes(fill=fillc,size=sub),stroke=.6,position=position_dodge(.55))+
  scale_size_manual(values=c(`Adult-onset`=1.6,`Childhood-onset`=1.9),guide="none")+   # 세모는 같은 size 값에서 더 커 보여 축소
  scale_color_manual(values=c(`Adult-onset`=RED,`Childhood-onset`=BLU),name="Asthma subtype")+
  scale_shape_manual(values=c(`Adult-onset`=24,`Childhood-onset`=22),name="Asthma subtype")+
  scale_fill_identity()+
  scale_y_discrete(labels=ABB)+
  labs(title="By asthma subtype",
       subtitle="filled = FDR<0.05",
       x=expression(genetic~correlation~r[g]),y=NULL)+
  th+theme(legend.position="bottom", legend.title=element_blank(),
           legend.box.spacing=unit(2,"pt"), legend.margin=margin(0,0,0,0), legend.key.height=unit(9,"pt"), legend.key.width=unit(20,"pt"), legend.spacing.x=unit(8,"pt"),   # 범례 키(선) 길이 확대
           axis.title.x=element_text(margin=margin(t=3)),
           plot.margin=margin(5,12,5,5))
## B 를 free() 로 정렬 해제 → A 의 기울어진 x 라벨 높이만큼 B 아래가 비는 문제 해소
f1 <- (f1a|free(f1b))+plot_layout(widths=c(0.95,1))+plot_annotation(tag_levels="A")   # B 패널 폭 추가 확보   # B 를 정렬에서 풀어 x축 제목·범례가 축 바로 아래에 오도록 (B 패널이 A 의 x 라벨 여백까지 사용)
save_fig("Fig1_shared_architecture", f1, FW, FW*0.5)

## ================= FIG 2 — Causal architecture (MR) ====================== ##
fw <- fread(file.path(B,"mr/mr_summary_forward.csv"))
rv <- fread(file.path(B,"mr/mr_summary_reverse.csv"))
fw[,dir:=paste0("Asthma ",ARROW," psychiatric")]; rv[,dir:=paste0("Psychiatric ",ARROW," asthma")]
mm <- rbindlist(list(fw,rv),fill=TRUE)
mm[,pair:=paste0(ABB[EXPOSURE]," ",ARROW," ",ABB[OUTCOME])]
mm[,sig:=IVW_P_FDR<0.05]
setorder(mm,dir,IVW_OR)
mm[,pair:=factor(pair,levels=unique(pair))]
f2a <- ggplot(mm,aes(IVW_OR,pair,color=sig))+
  geom_vline(xintercept=1,linetype=2,color="grey55")+
  geom_errorbarh(aes(xmin=IVW_LCI,xmax=IVW_UCI),height=.22)+
  geom_point(size=2.2,shape=15)+   # v9: 네모
  facet_wrap(~dir,ncol=1,scales="free_y")+
  scale_color_manual(values=c(`TRUE`=RED,`FALSE`=GREY),breaks=c(TRUE,FALSE),labels=c("FDR < 0.05","Not significant"),name=NULL)+
  scale_x_log10()+
  labs(title="Bidirectional MR estimates (IVW)",
       x="OR (95% CI) per unit increase in log-odds of exposure",y=NULL)+
  th+theme(strip.text=element_text(face="bold",hjust=0),strip.background=element_blank(),
           legend.position=c(0.985,0.985),legend.justification=c(1,1),legend.direction="vertical",   # 범례: A 패널 안 오른쪽 위, 한 번만
           legend.background=element_rect(fill="white",color=NA),legend.margin=margin(2,4,2,4),legend.key.size=unit(9,"pt"))

## 2b non-mediation: depression<->lung null vs asthma<->lung strong
lg <- fread(file.path(B,"mr/mr_summary_lung_ALL.csv"))
sel <- lg[(TAG=="lung_mental2lung"&OUTCOME=="LUNG") |
          (TAG=="lung_lung2mental"&EXPOSURE=="LUNG") |
          (TAG=="lung_asthma2lung"&OUTCOME=="LUNG") |
          (TAG=="lung_lung2asthma"&EXPOSURE=="LUNG")]
sel[,b:=IVW_B][,lci:=IVW_B_LCI][,uci:=IVW_B_UCI]
sel[,pair:=paste0(ABB[EXPOSURE]," ",ARROW," ",ABB[OUTCOME])]
sel[,grp:=ifelse(OUTCOME!="LUNG",
                 "Outcome: asthma or psychiatric trait (log OR)","Outcome: FEV1/FVC ratio (beta, SD)")]
sel[,sig:=IVW_P_FDR<0.05]
setorder(sel,grp,b); sel[,pair:=factor(pair,levels=unique(pair))]
f2b <- ggplot(sel,aes(b,pair,color=sig))+
  geom_vline(xintercept=0,linetype=2,color="grey55")+
  geom_errorbarh(aes(xmin=lci,xmax=uci),height=.22)+
  geom_point(size=2.1,shape=15)+   # v9: 네모
  facet_wrap(~grp,ncol=1,scales="free")+
  scale_color_manual(values=c(`TRUE`=RED,`FALSE`=GREY),breaks=c(TRUE,FALSE),labels=c("FDR < 0.05","Not significant"),name=NULL)+
  labs(title="Pairwise MR with FEV1/FVC ratio",
       x="MR estimate (95% CI)",y=NULL)+
  th+theme(strip.text=element_text(face="bold",hjust=0),strip.background=element_blank(),legend.position="none")   # 범례는 A 에만
f2 <- (f2a/f2b)+plot_layout(heights=c(1,1.05))+plot_annotation(tag_levels="A")
save_fig("Fig2_causal_MR", f2, FW*0.9, FW*0.9*8.6/9.2)

## ================= FIG 3 — Pleiotropy & shared loci ====================== ##
pl <- fread(file.path(B,"placo/placo_summary.csv"))
pl[,pair:=paste0(LAB[EXPOSURE],"\n x ",LAB[OUTCOME])]
pl[,dom:=fifelse(OUTCOME%in%c("MDD","ANX"),"Depression/anxiety",
          fifelse(OUTCOME%in%c("BIP","SCZ"),"Other psychiatric","Lung function"))]
pl[,sub:=factor(EXPOSURE,levels=c("COA","AOA"),labels=c("Childhood-onset","Adult-onset"))]   # v14
pl[,olab:=factor(OUTCOME,levels=c("MDD","ANX","SCZ","BIP","LUNG"))]   # x 순서: Dep > Anx > SCZ > BIP > FEV1/FVC
f3a <- ggplot(pl,aes(olab,N_LOCI,fill=dom))+
  geom_col(width=.7,color=NA,position=position_dodge())+
  geom_text(aes(label=N_LOCI),vjust=-0.3,size=BS*0.3,family=FF,position=position_dodge(.7))+
  facet_wrap(~sub,nrow=1)+
  scale_fill_manual(values=c(`Depression/anxiety`=RED,`Lung function`="grey65",`Other psychiatric`=BLU),
     breaks=c("Depression/anxiety","Other psychiatric","Lung function"),name=NULL)+   # 범례 순서
  scale_x_discrete(labels=ABB)+ylim(0,max(pl$N_LOCI)*1.15)+
  labs(title="Pleiotropic regions with asthma (exploratory PLACO screen)",
       x=NULL,y=expression(N~regions~"(500-kb rule, PLACO"~italic(P)<5%*%10^-8*")"))+
  th+theme(axis.text.x=element_text(angle=0,hjust=0.5),legend.position="bottom",   # 약어라 회전 불필요
           legend.box.spacing=unit(6,"pt"),legend.margin=margin(0,0,0,0),   # A(축 제목 아래 간격)와 시각적으로 같은 간격
           legend.key.size=unit(11,"pt"),legend.text=element_text(size=BS))   # A·B 범례 크기 통일

## 3b LAVA local rg — immune-locus highlight
lv <- fread(file.path(B,"lava/lava_sig_annotated.csv"))
imm <- "(^|,)(IL[0-9]|IL[0-9][0-9]|HLA|TYK2|TNF|IKBK|IFI|IRF|STAT|JAK|TSLP|CRP|CD[0-9]|NOD2|NFKB|ICAM|TLR)"
lv[,immune:=grepl(imm,genes,ignore.case=TRUE)]
lv[,pairf:=factor(pair)]
cnt <- lv[,.(n=.N, n_imm=sum(immune)),by=pair][order(-n)]
cnt <- melt(cnt,id.vars="pair",measure.vars=c("n_imm"),value.name="n_imm")
tot <- lv[,.(n=.N),by=pair]
cnt <- merge(tot,lv[,.(n_imm=sum(immune)),by=pair],by="pair")
PRL <- function(x){ p<-strsplit(x,"_")[[1]]
  s<-ABB[p[1]]; t<-ABB[p[2]]
  paste0(s," \u00d7 ",t) }
cnt[,plab:=sapply(as.character(pair),PRL)]
cnt[,pair:=factor(plab,levels=cnt[order(n)]$plab)]
cnt2 <- melt(cnt,id.vars="pair",measure.vars=c("n","n_imm"),
             variable.name="type",value.name="val")
cnt2[,type:=factor(type,levels=c("n","n_imm"),labels=c("All sig loci","Immune-gene loci"))]
f3b <- ggplot(cnt2,aes(val,pair,fill=type))+
  geom_col(data=cnt2[type=="All sig loci"],width=.62)+
  geom_col(data=cnt2[type=="Immune-gene loci"],width=.62)+
  scale_fill_manual(values=c(`All sig loci`="grey80",`Immune-gene loci`=RED),
     breaks=c("All sig loci","Immune-gene loci"),
     labels=c("All FDR-significant blocks","Blocks containing immune-related genes"),name=NULL)+
  geom_text(data=cnt[,.(pair,n)],aes(n,pair,label=n),inherit.aes=FALSE,hjust=-0.3,size=BS*0.3,family=FF)+
  scale_x_continuous(expand=expansion(mult=c(0,0.12)))+
  labs(title="Local genetic correlations (LAVA, FDR < 0.05)",
       x=expression(N~FDR*"-significant"~local~r[g]~blocks),y=NULL)+
  th+theme(legend.position="bottom",legend.box.spacing=unit(2,"pt"),legend.margin=margin(0,0,0,0),
           legend.key.size=unit(11,"pt"),legend.text=element_text(size=BS))   # A·B 범례 크기 통일
## v10.1 (2026-09-24): 패널 순서 교체 — A = LAVA(주 분석), B = PLACO(탐색적). 본문·legend 인용 순서(Fig 3, A → B)와 일치
## B(PLACO)의 y축 제목이 A의 긴 y 라벨에 맞춰 멀리 떨어지지 않도록 free() 로 정렬 해제 (patchwork ≥ 1.2)
f3 <- (f3b/free(f3a, type="label", side="l"))+plot_layout(heights=c(1,1))+plot_annotation(tag_levels="A")   # 패널 폭은 A 와 정렬, y축 제목만 축에 붙임
save_fig("Fig3_pleiotropy_loci", f3, FW*0.9, FW*0.9*8.2/9.5)

## ================= FIG 4 — cis-MR vs colocalization ====================== ##
tm <- fread(file.path(B,"mr/mr_tissue_cis.csv"))
tm <- tm[MATCHED==TRUE & TISSUE%in%c("Lung","Brain_Cortex")]
tm[,lab:=paste0(GENE," [",ifelse(TISSUE=="Lung","Lung","Brain"),"] ",ARROW," ",OUTCOME)]
tm[,sig:=!is.na(P_FDR)&P_FDR<0.05]
setorder(tm,OR); tm[,lab:=factor(lab,levels=lab)]
f4a <- ggplot(tm,aes(OR,lab,color=sig))+
  geom_vline(xintercept=1,linetype=2,color="grey55")+
  geom_errorbarh(aes(xmin=OR_LCI,xmax=OR_UCI),height=.22)+
  geom_point(size=2.1,shape=15)+   # v9: 네모
  scale_color_manual(values=c(`TRUE`=RED,`FALSE`=GREY),breaks=c(TRUE,FALSE),labels=c("FDR < 0.05","Not significant"),name=NULL)+
  scale_x_log10()+
  labs(title="Tissue cis-MR (single-instrument)",
       x="OR (95% CI) per unit normalized expression",y=NULL)+
  th+theme(legend.position="bottom",legend.box.spacing=unit(2,"pt"),legend.margin=margin(0,0,0,0),legend.key.size=unit(9,"pt"))

co <- fread(file.path(B,"coloc/coloc_summary.csv"))
co[,lab:=paste0(GENE," [",ifelse(TISSUE=="Lung","Lung","Brain"),"] ",ARROW," ",OUTCOME)]
co[,lab:=factor(lab,levels=co[order(PP_H4)]$lab)]
co[,strong:=PP_H4>=0.8]
f4b <- ggplot(co,aes(PP_H4,lab))+
  geom_col(width=.62,fill="grey70",color="grey30")+
  geom_vline(xintercept=0.8,linetype=2,color=RED)+
  geom_text(aes(label=sprintf("%.2f",PP_H4)),hjust=-0.2,size=BS*0.29,family=FF)+
  scale_x_continuous(limits=c(0,1.05),breaks=seq(0,1,.2),expand=expansion(mult=c(0,0.05)))+
  labs(title="Colocalization (coloc, PP.H4)",
       x="PP.H4 (shared causal variant)",y=NULL)+th
f4 <- (f4a|plot_spacer()|f4b)+plot_layout(widths=c(1,0.2,1))+plot_annotation(tag_levels="A")   # A–B 사이 간격
save_fig("Fig4_cisMR_vs_coloc", f4, FW, FW*4.8/11.5)

cat("DONE\n"); print(list.files(O)); cat("FINAL:", FINAL, "\n"); print(list.files(FINAL))
