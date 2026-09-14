## Manuscript figures — reframed story (robust results)
## 실행: source("Step9_figures.R")  (프로젝트 폴더 결과 CSV들을 읽어 Result/ 에 Fig1-4 저장)
suppressPackageStartupMessages({
  for(p in c("data.table","ggplot2","patchwork")) if(!requireNamespace(p,quietly=TRUE)) install.packages(p,repos="https://cloud.r-project.org")
  library(data.table); library(ggplot2); library(patchwork)})
B <- "PATH/TO/PROJECT"   # 프로젝트 루트
O <- file.path(B,"Result"); dir.create(O, showWarnings=FALSE, recursive=TRUE)
RED<-"#B2182B"; BLU<-"#2166AC"; GREY<-"grey70"
LAB <- c(COA="Childhood-onset asthma", AOA="Adult-onset asthma", MDD="Depression",
         ANX="Anxiety", BIP="Bipolar", SCZ="Schizophrenia",
         FEV1="FEV1", FVC="FVC", LUNG="FEV1/FVC ratio")
## JACI: 그림 내 글꼴 Times New Roman, 형식 TIFF/JPG ≥300 dpi
FF <- "serif"   # Windows 에서 serif = Times New Roman
th <- theme_minimal(base_size=11, base_family=FF) +
  theme(panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold",size=12,family=FF),
        plot.tag=element_text(face="bold",size=14,family=FF),
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
hm[,rglab:=ifelse(is.na(rg),"",sprintf("%.2f%s",rg,star))]
f1a <- ggplot(hm,aes(V1,V2,fill=rg))+
  geom_tile(color="white",linewidth=.6)+
  geom_text(aes(label=rglab),size=2.7,family=FF)+
  scale_fill_gradient2(low=BLU,mid="white",high=RED,midpoint=0,limits=c(-1,1),
     na.value="grey92",name=expression(r[g]))+
  scale_x_discrete(labels=LAB)+scale_y_discrete(labels=LAB)+
  labs(title="Global genetic correlations (LDSC)",x=NULL,y=NULL,
       subtitle="* FDR<0.05")+
  th+theme(axis.text.x=element_text(angle=40,hjust=1),legend.position="right",
           panel.grid.major=element_blank())

## 1b forest — asthma x trait (subtype specificity)
fo <- rg[focus=="asthma_x_trait"]
fo[,lci:=rg-1.96*SE][,uci:=rg+1.96*SE]
fo[,trait:=factor(T2,levels=c("MDD","ANX","BIP","SCZ","LUNG","FVC","FEV1"))]
fo[,sub:=factor(T1,levels=c("AOA","COA"),labels=c("Adult-onset","Childhood-onset"))]
fo[,sig:=P_FDR<0.05]
f1b <- ggplot(fo,aes(rg,trait,color=sub,shape=sig))+
  geom_vline(xintercept=0,linetype=2,color="grey55")+
  geom_errorbarh(aes(xmin=lci,xmax=uci),height=.25,position=position_dodge(.55))+
  geom_point(size=2.6,position=position_dodge(.55))+
  scale_color_manual(values=c(`Adult-onset`=RED,`Childhood-onset`=BLU),name="Asthma subtype")+
  scale_shape_manual(values=c(`TRUE`=19,`FALSE`=1),guide="none")+
  scale_y_discrete(labels=LAB)+
  labs(title="Genetic correlation by asthma subtype",
       subtitle="filled = FDR<0.05",
       x=expression(genetic~correlation~r[g]),y=NULL)+
  th+theme(legend.position="bottom")
f1 <- (f1a|f1b)+plot_layout(widths=c(1.25,1))+plot_annotation(tag_levels="a")
ggsave(file.path(O,"Fig1_shared_architecture.png"),f1,width=13,height=5.6,dpi=300,bg="white")
ggsave(file.path(O,"Fig1_shared_architecture.tiff"),f1,width=13,height=5.6,dpi=300,bg="white",compression="lzw")

## ================= FIG 2 — Causal architecture (MR) ====================== ##
fw <- fread(file.path(B,"mr/mr_summary_forward.csv"))
rv <- fread(file.path(B,"mr/mr_summary_reverse.csv"))
fw[,dir:="Asthma -> psychiatric"]; rv[,dir:="Psychiatric -> asthma"]
mm <- rbindlist(list(fw,rv),fill=TRUE)
mm[,pair:=paste0(LAB[EXPOSURE]," -> ",LAB[OUTCOME])]
mm[,sig:=IVW_P_FDR<0.05]
setorder(mm,dir,IVW_OR)
mm[,pair:=factor(pair,levels=unique(pair))]
f2a <- ggplot(mm,aes(IVW_OR,pair,color=sig))+
  geom_vline(xintercept=1,linetype=2,color="grey55")+
  geom_errorbarh(aes(xmin=IVW_LCI,xmax=IVW_UCI),height=.22)+
  geom_point(size=2.7)+
  facet_wrap(~dir,ncol=1,scales="free_y")+
  scale_color_manual(values=c(`TRUE`=RED,`FALSE`=GREY),guide="none")+
  scale_x_log10()+
  labs(title="Bidirectional MR estimates (IVW)",
       subtitle="red = FDR<0.05 within direction; OR per unit log-odds of exposure",
       x="OR (95% CI) per unit increase in log-odds of exposure",y=NULL)+
  th+theme(strip.text=element_text(face="bold",hjust=0),strip.background=element_blank())

## 2b non-mediation: depression<->lung null vs asthma<->lung strong
lg <- fread(file.path(B,"mr/mr_summary_lung_ALL.csv"))
sel <- lg[(TAG=="lung_mental2lung"&OUTCOME=="LUNG") |
          (TAG=="lung_lung2mental"&EXPOSURE=="LUNG") |
          (TAG=="lung_asthma2lung"&OUTCOME=="LUNG") |
          (TAG=="lung_lung2asthma"&EXPOSURE=="LUNG")]
sel[,b:=IVW_B][,lci:=IVW_B_LCI][,uci:=IVW_B_UCI]
sel[,pair:=paste0(LAB[EXPOSURE]," -> ",LAB[OUTCOME])]
sel[,grp:=ifelse(OUTCOME!="LUNG",
                 "Outcome: asthma or psychiatric trait (log OR)","Outcome: FEV1/FVC ratio (beta, SD)")]
sel[,sig:=IVW_P_FDR<0.05]
setorder(sel,grp,b); sel[,pair:=factor(pair,levels=unique(pair))]
f2b <- ggplot(sel,aes(b,pair,color=sig))+
  geom_vline(xintercept=0,linetype=2,color="grey55")+
  geom_errorbarh(aes(xmin=lci,xmax=uci),height=.22)+
  geom_point(size=2.6)+
  facet_wrap(~grp,ncol=1,scales="free")+
  scale_color_manual(values=c(`TRUE`=RED,`FALSE`=GREY),guide="none")+
  labs(title="Pairwise MR with FEV1/FVC ratio",
       subtitle="red = FDR<0.05; panels separate continuous and binary outcomes",
       x="MR estimate (95% CI)",y=NULL)+
  th+theme(strip.text=element_text(face="bold",hjust=0),strip.background=element_blank())
f2 <- (f2a/f2b)+plot_layout(heights=c(1,1.05))+plot_annotation(tag_levels="a")
ggsave(file.path(O,"Fig2_causal_MR.png"),f2,width=9.2,height=8.6,dpi=300,bg="white")
ggsave(file.path(O,"Fig2_causal_MR.tiff"),f2,width=9.2,height=8.6,dpi=300,bg="white",compression="lzw")

## ================= FIG 3 — Pleiotropy & shared loci ====================== ##
pl <- fread(file.path(B,"placo/placo_summary.csv"))
pl[,pair:=paste0(LAB[EXPOSURE],"\n x ",LAB[OUTCOME])]
pl[,dom:=fifelse(OUTCOME%in%c("MDD","ANX"),"Depression/anxiety",
          fifelse(OUTCOME%in%c("BIP","SCZ"),"Other psychiatric","Lung function"))]
pl[,sub:=factor(EXPOSURE,levels=c("COA","AOA"),labels=c("Childhood","Adult-onset"))]
pl[,olab:=factor(OUTCOME,levels=c("MDD","ANX","BIP","SCZ","LUNG"))]
f3a <- ggplot(pl,aes(olab,N_LOCI,fill=dom))+
  geom_col(width=.7,color="grey30",position=position_dodge())+
  geom_text(aes(label=N_LOCI),vjust=-0.3,size=3,family=FF,position=position_dodge(.7))+
  facet_wrap(~sub,nrow=1)+
  scale_fill_manual(values=c(`Depression/anxiety`=RED,`Other psychiatric`="#E39A9A",
     `Lung function`=BLU),name=NULL)+
  scale_x_discrete(labels=LAB)+ylim(0,max(pl$N_LOCI)*1.15)+
  labs(title="Pleiotropic regions with asthma (exploratory PLACO screen, 500-kb regions)",
       x=NULL,y="N regions")+
  th+theme(axis.text.x=element_text(angle=30,hjust=1),legend.position="bottom")

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
  s<-c(COA="Childhood",AOA="Adult-onset")[p[1]]
  t<-c(MDD="Depression",ANX="Anxiety",BIP="Bipolar",SCZ="Schizophrenia",LUNG="FEV1/FVC")[p[2]]
  paste0(s," x ",t) }
cnt[,plab:=sapply(as.character(pair),PRL)]
cnt[,pair:=factor(plab,levels=cnt[order(n)]$plab)]
cnt2 <- melt(cnt,id.vars="pair",measure.vars=c("n","n_imm"),
             variable.name="type",value.name="val")
cnt2[,type:=factor(type,levels=c("n","n_imm"),labels=c("All sig loci","Immune-gene loci"))]
f3b <- ggplot(cnt2,aes(val,pair,fill=type))+
  geom_col(data=cnt2[type=="All sig loci"],width=.62,fill="grey80")+
  geom_col(data=cnt2[type=="Immune-gene loci"],width=.62,fill=RED)+
  geom_text(data=cnt[,.(pair,n)],aes(n,pair,label=n),inherit.aes=FALSE,hjust=-0.3,size=3,family=FF)+
  scale_x_continuous(expand=expansion(mult=c(0,0.12)))+
  labs(title="Local genetic correlations (LAVA)",
       subtitle="grey = FDR<0.05 local rg; red = blocks with annotated immune-related genes (descriptive)",
       x="N FDR-significant local rg blocks",y=NULL)+th
f3 <- (f3a/f3b)+plot_layout(heights=c(1,1))+plot_annotation(tag_levels="a")
ggsave(file.path(O,"Fig3_pleiotropy_loci.png"),f3,width=9.5,height=8.2,dpi=300,bg="white")
ggsave(file.path(O,"Fig3_pleiotropy_loci.tiff"),f3,width=9.5,height=8.2,dpi=300,bg="white",compression="lzw")

## ================= FIG 4 — cis-MR vs colocalization ====================== ##
tm <- fread(file.path(B,"mr/mr_tissue_cis.csv"))
tm <- tm[MATCHED==TRUE & TISSUE%in%c("Lung","Brain_Cortex")]
tm[,lab:=paste0(GENE," [",ifelse(TISSUE=="Lung","Lung","Brain"),"] -> ",OUTCOME)]
tm[,sig:=!is.na(P_FDR)&P_FDR<0.05]
setorder(tm,OR); tm[,lab:=factor(lab,levels=lab)]
f4a <- ggplot(tm,aes(OR,lab,color=sig))+
  geom_vline(xintercept=1,linetype=2,color="grey55")+
  geom_errorbarh(aes(xmin=OR_LCI,xmax=OR_UCI),height=.22)+
  geom_point(size=2.6)+
  scale_color_manual(values=c(`TRUE`=RED,`FALSE`=GREY),guide="none")+
  scale_x_log10()+
  labs(title="Tissue cis-MR (single-instrument)",
       subtitle="red = FDR<0.05 within matched tissue-outcome comparisons",
       x="OR (95% CI) per unit normalized expression",y=NULL)+th

co <- fread(file.path(B,"coloc/coloc_summary.csv"))
co[,lab:=paste0(GENE," [",ifelse(TISSUE=="Lung","Lung","Brain"),"] -> ",OUTCOME)]
co[,lab:=factor(lab,levels=co[order(PP_H4)]$lab)]
co[,strong:=PP_H4>=0.8]
f4b <- ggplot(co,aes(PP_H4,lab))+
  geom_col(width=.62,fill="grey70",color="grey30")+
  geom_vline(xintercept=0.8,linetype=2,color=RED)+
  geom_text(aes(label=sprintf("%.2f",PP_H4)),hjust=-0.2,size=2.9,family=FF)+
  scale_x_continuous(limits=c(0,1.05),breaks=seq(0,1,.2),expand=expansion(mult=c(0,0.05)))+
  labs(title="Colocalization (coloc, PP.H4)",
       subtitle="dashed line = PP.H4 0.8; no tested pair met the threshold",
       x="PP.H4 (shared causal variant)",y=NULL)+th
f4 <- (f4a|f4b)+plot_annotation(tag_levels="a")
ggsave(file.path(O,"Fig4_cisMR_vs_coloc.png"),f4,width=11.5,height=4.8,dpi=300,bg="white")
ggsave(file.path(O,"Fig4_cisMR_vs_coloc.tiff"),f4,width=11.5,height=4.8,dpi=300,bg="white",compression="lzw")

cat("DONE\n"); print(list.files(O))
