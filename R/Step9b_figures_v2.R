## Step9b — 제출용 그림 v2 (Fig 1에 C 패널 추가, Fig 3A를 LAVA ideogram으로)
## v2.1 (2026-10-10): Fig 3A 위/아래 방향 수정(성인기 = 위), Fig 1C P 값을 JACI 형식(.xx, 반올림 half-up)으로; 그래픽 초록 코드는 제거(사용자 제작본 사용).
## 실행: 콘솔에서 source("R/Step9b_figures_v2.R")
## Step9_figures.R 를 먼저 source 하여 f1a, f1b, f3a(PLACO), lv(LAVA) 등을 그대로 재사용한다 (Step9 출력도 함께 재생성됨).
source("R/Step9_figures.R")
for(p in c("ggrepel")) if(!requireNamespace(p,quietly=TRUE)) install.packages(p,repos="https://cloud.r-project.org")
library(ggrepel)
OUT2 <- file.path(B,"Submission","01 JACI","v3.1","Figures_v2"); dir.create(OUT2, showWarnings=FALSE, recursive=TRUE)
save2 <- function(name, plot, w, h){
  ggsave(file.path(OUT2,paste0(name,".png")),  plot, width=w, height=h, dpi=300, bg="white", device=ragg::agg_png)
  ggsave(file.path(OUT2,paste0(name,".tiff")), plot, width=w, height=h, dpi=600, bg="white", device=ragg::agg_tiff, compression="lzw")
  ggsave(file.path(OUT2,paste0(name,".pdf")),  plot, width=w, height=h, bg="white", device=grDevices::cairo_pdf)
}
ORA <- "#E08214"; BLU2 <- "#5AAE61"   # 보조색 (불안 = 주황, 조현병 = 진파랑은 BLU 재사용)

## ================= FIG 1 v2: A heatmap | B forest | C 아형 대비 (Δrg) ================= ##
s9 <- fread(file.path(B,"Result","tables_csv","TableS9_subtype_contrasts.csv"))
ca <- s9[grepl("^A\\.",Panel) & grepl("^r_g\\(AOA",Contrast)]
ca[,trait:=c("MDD","ANX","BIP","SCZ")]
ca[,est:=as.numeric(sub(" .*","",`Estimate (95% CI)`))]
ca[,lci:=as.numeric(sub(".*\\((-?[0-9.]+) to.*","\\1",`Estimate (95% CI)`))]
ca[,uci:=as.numeric(sub(".*to (-?[0-9.]+)\\).*","\\1",`Estimate (95% CI)`))]
ca[,P:=as.numeric(P)][,sig:=P<0.05]
r2h <- function(x) floor(x*100+0.5+1e-9)/100   # half-up (0.295 -> .30; sprintf 은 .29 로 내림)
ca[,plab:=ifelse(P<0.001, sprintf("italic(P)==%s%%*%%10^%d", formatC(P/10^floor(log10(P)),digits=1,format="f"), floor(log10(P))),
                 sprintf("italic(P)=='%s'", sub("^0","",sprintf("%.2f",r2h(P)))))]   # 따옴표: '.30' 을 숫자로 파싱해 0.3 으로 찍히지 않게
ca[,trait:=factor(trait,levels=rev(c("MDD","ANX","BIP","SCZ")))]
f1c <- ggplot(ca,aes(est,trait))+
  geom_vline(xintercept=0,linetype=2,color="grey55")+
  geom_errorbarh(aes(xmin=lci,xmax=uci,color=sig),height=.22)+
  geom_point(aes(color=sig),size=2.2,shape=15)+
  geom_text(aes(x=0.33,label=plab),parse=TRUE,hjust=0,size=BS*0.28,family=FF)+
  scale_color_manual(values=c(`TRUE`=RED,`FALSE`=GREY),guide="none")+
  scale_x_continuous(limits=c(-0.12,0.62),breaks=seq(-0.1,0.3,0.1))+
  labs(title="Subtype contrast",
       x=expression(Delta*r[g]~"(adult-onset − childhood-onset)"),y=NULL)+
  th+theme(plot.margin=margin(5,5,5,5))
## Step9 에서 free() 를 거친 객체(f1b, f3a)에 free 설정 속성이 남아 중첩 배치 시 patchwork 오류(simplify_gt.free_table) → 속성 제거
unfree <- function(p){ class(p) <- setdiff(class(p), c("free_plot","free_table")); for(a in grep("free", names(attributes(p)), value=TRUE)) attr(p,a) <- NULL; p }
f1a <- unfree(f1a); f1b <- unfree(f1b); f3a <- unfree(f3a)
right <- (f1b / f1c) + plot_layout(heights=c(1.55,1))
f1v2 <- (f1a | right) + plot_layout(widths=c(1.05,1)) + plot_annotation(tag_levels="A")
save2("Fig1_shared_architecture_v2", f1v2, FW, FW*0.64)

## ================= FIG 3 v2: A LAVA ideogram | B PLACO ================= ##
## hg19 염색체 길이 (LAVA 블록 좌표 기준)
chrlen <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,
            135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566)
chr <- data.table(CHR=1:22, len=chrlen)
TRC <- c(MDD=RED, ANX=ORA, SCZ=BLU, BIP="#80B1D3", LUNG="grey55")
lv2 <- copy(lv)
lv2[,c("sub","trt"):=tstrsplit(pair,"_")]
lv2[,mid:=(START+STOP)/2]
lv2[,trait:=factor(trt,levels=c("MDD","ANX","SCZ","BIP","LUNG"))]
lv2[,sublab:=factor(sub,levels=c("AOA","COA"),labels=c("Adult-onset","Childhood-onset"))]
lv2[,yoff:=ifelse(sub=="AOA",-0.28,0.28)]            # scale_y_reverse 이므로 음수가 위: 성인기 = 염색체 위, 소아기 = 아래 (v2.1 수정)
lv2[,y:=CHR+yoff]
lv2[,sign:=ifelse(rho>0,"Positive","Negative")]
lv2[,fillc:=ifelse(rho>0, TRC[trt], "white")]
## 면역 유전자 라벨: 블록 내 면역 관련 유전자만 추출
immg <- function(g){ x<-strsplit(g,",")[[1]]; x<-x[grepl(imm,paste0(",",x),ignore.case=TRUE)]; if(!length(x)) NA_character_ else paste(head(x,2),collapse=", ") }
lv2[,glab:=sapply(genes,immg)]
f3i <- ggplot()+
  geom_segment(data=chr,aes(x=0,xend=len/1e6,y=CHR,yend=CHR),linewidth=2.6,color="grey88",lineend="round")+
  geom_point(data=lv2,aes(mid/1e6,y,color=trait,shape=sublab,fill=fillc),size=1.9,stroke=.5)+
  geom_text_repel(data=lv2[!is.na(glab)],aes(mid/1e6,y,label=glab),size=BS*0.26,family=FF,fontface="italic",
                  min.segment.length=0,segment.size=.25,segment.color="grey50",box.padding=.45,point.padding=.35,force=2.5,max.overlaps=Inf,seed=3)+
  scale_color_manual(values=TRC,labels=c(MDD="Depression",ANX="Anxiety",SCZ="Schizophrenia",BIP="Bipolar disorder",LUNG="FEV1/FVC ratio"),name=NULL)+
  scale_shape_manual(values=c(`Adult-onset`=24,`Childhood-onset`=22),name=NULL)+
  scale_fill_identity()+
  scale_y_reverse(breaks=1:22,labels=paste0("chr",1:22),expand=expansion(add=0.6))+
  scale_x_continuous(expand=expansion(mult=c(0.01,0.03)))+
  guides(color=guide_legend(order=1,override.aes=list(shape=15,size=2.6)),
         shape=guide_legend(order=2,override.aes=list(fill="grey40",size=2.2)))+
  labs(title="Local genetic correlations with asthma (LAVA, FDR < 0.05)",
       subtitle="Above chromosome, adult-onset; below, childhood-onset. Filled, positive; open, negative local correlation. Immune-related genes labeled.",
       x="Position (Mb)",y=NULL)+
  th+theme(legend.position="bottom",legend.box="horizontal",legend.box.spacing=unit(2,"pt"),legend.margin=margin(0,0,0,0),
           legend.key.size=unit(10,"pt"),legend.text=element_text(size=BS-1),plot.subtitle=element_text(size=BS-1,color="grey30"),
           panel.grid.major.y=element_blank(),panel.grid.major.x=element_line(color="grey94"))
f3v2 <- (f3i/f3a)+plot_layout(heights=c(1.45,1))+plot_annotation(tag_levels="A")
save2("Fig3_pleiotropy_loci_v2", f3v2, FW*0.9, FW*0.9*1.25)

cat("Figures_v2 saved to:", OUT2, "\n")
