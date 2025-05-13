rm(list=ls());gc()
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper/')
setwd('GitHub_repo/tracing/')
library(data.table);library(ggplot2);library(openxlsx);library(tidyr);library(dplyr)
dat<-read.xlsx('20240604 Updated Fig 1D Details for Ling.xlsx');setDT(dat)
dat[,source:=gsub(' (','(',fixed = T,source)]
dat%>%separate(source,into=c('origin','type'),remove = F,sep = ' ')%>%setDT->dat
dat[type=='PDX(P2)',type:='PDX(P2/3)']
dat[,source:=factor(source,levels=unique(source))]
dat

dat[,origin:=factor(origin,levels=unique(dat$origin))]
# compare P2(/3) to P6 within each line
stat.df<-dat[type!='Pt',list(pv=signif(t.test(ratio~type)$p.value,2)),by=origin]
ori.col<-readRDS('../batch2/ori_col.rds')

summary(with(dat,aov(ratio~origin*type)))
summary(with(dat,aov(ratio~origin+type)))
aov.pv<-signif(as.matrix(summary(with(dat,aov(ratio~origin+type)))[[1]])[1:2,5],2)
names(aov.pv)<-trimws(names(aov.pv))
pdf('Ext2A.pdf',w=5.1,h=2.9)
ggplot(dat,aes(x=source,y=ratio,color=origin))+facet_grid(~origin,space = 'free',scales = 'free_x')+geom_point()+  
  stat_summary(fun = "mean",geom = "crossbar", width = 0.5)+ylab('Citrate (M+2) /Pyruvate (M+3)')+
  scale_alpha_manual(values=c(Pt=1,`PDX(P2)`=.5,`PDX(P2/3)`=.5,`PDX(P6)`=.3))+
  scale_color_manual(values=ori.col[levels(dat$origin)])+theme_bw()+
  geom_segment(x=2,xend=3,y=1.6,color='black')+geom_text(stat.df,x=2.5,y=1.75,color='black',mapping=aes(label=paste0('pv = ',pv)))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5),strip.background = element_rect(fill = NA),axis.title.x = element_blank(),
        panel.spacing = unit(0,'lines'),legend.position = 'none')+ylim(.3,2)+
  ggtitle(paste0("Two-way ANOVA: labeling ratio ~ sample type + origin\nsample type pv = ",aov.pv['type'],', origin pv = ',aov.pv['origin']))
dev.off()
