rm(list=ls());gc()
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper/')
TCA.df<-readRDS('work/isotope/tca_df_with_rep.rds')
fwrite(TCA.df,'work/isotope/S1_input.csv')
TCA.df<-TCA.df[sample!='Patient plasma'][grep("P6|P3",source,invert = T)][grepl("MP10|MP4A|MP8|MP9",source)][source!='MP9_Tumor']
TCA.m<-c('GLUCOSE_M6','3PG_M3','PYRUVATE_M3','LACTATE_M3','ALANINE_M3','CITRATE_M2','GLUTAMATE_M2','SUCCINATE_M2','FUMARATE_M2','MALATE_M2','ASPARTATE_M2')
ori.col<-readRDS('work/batch2/ori_col.rds')
# TCA.df<-TCA.df[!grepl('3h',sample_name)]

TCA.df[,sample_name:=gsub("[1-3]$",'',as.character(sample_name))]
TCA.df[,.N,by=.(source,sample_name)]
TCA.df<-TCA.df[,median(value),by=.(source,feature,metabolite,labeling,sample_name,sample)]
TCA.df[,value:=V1][,V1:=NULL]
TCA.df[,origin:=sapply(source,function(x) unlist(strsplit(x,split = '_'))[1])]
TCA.df[,type:=sample][,sample:=NULL]
library(dplyr)
library(ggplot2)

TCA.df.norm.div<-TCA.df[feature!='GLUCOSE_M6']
ref<-TCA.df[feature=='GLUCOSE_M6']
for(f in unique(TCA.df.norm.div$feature)) TCA.df.norm.div[feature==f,value:=value/ref$value]

TCA.df.div.m<-TCA.df.norm.div[,median(value),by=.(source,feature,origin,type)]
TCA_summary <- TCA.df.norm.div %>%
  filter(type == "PDX") %>%
  group_by(origin, feature) %>%
  summarise(
    median_value = median(value, na.rm = TRUE),
    lower_quartile = quantile(value, 0.25, na.rm = TRUE),
    upper_quartile = quantile(value, 0.75, na.rm = TRUE)
  ) %>%
  ungroup()
setDT(TCA_summary);TCA_summary[,type:='PDX']
# TCA.m<-c(TCA.m,'CitM2/PyrM3')
TCA.df.div.m[,origin:=factor(origin,levels=intersect(names(ori.col),origin))]
TCA.df.div.m[,feature:=factor(feature,levels = TCA.m[-1])]
# t-test for paired difference
TCA.df.div.m.pt<-TCA.df.div.m[type=='Patient tumor']
TCA.df.div.m.pdx<-TCA.df.div.m[type=='PDX']
TCA.df.div.m.pt[,diff:=V1-TCA.df.div.m.pdx[match(with(TCA.df.div.m.pt,paste(origin,feature)),with(TCA.df.div.m.pdx,paste(origin,feature)))][['V1']]]
stat.df<-TCA.df.norm.div[,as.list(structure(signif(as.matrix(summary(aov(value~origin*type))[[1]])[1:3,5],2),names=c('origin','type','origin:type'))),by=feature]
stat.df<-TCA.df.norm.div[,as.list(structure(signif(as.matrix(summary(aov(value~origin+type))[[1]])[1:2,5],2),names=c('origin','type'))),by=feature]
stat.df
pdf('output/isotope/FigS1A',w=12,h=3.5)
ggplot()+geom_point(TCA.df.div.m,mapping=aes(x=type,y=V1,color=origin))+geom_line(TCA.df.div.m,mapping=aes(x=type,y=V1,group=origin,color=origin))+
  geom_errorbar(TCA_summary[feature!='GLUCOSE_M6'],mapping=aes(x=type,ymin = lower_quartile, ymax = upper_quartile,color=origin), width = 0.2)+
  geom_text(stat.df,mapping=aes(x=.5,y=1.3,label=paste('type pv:',type)),hjust=0,size=3)+
  geom_text(stat.df,mapping=aes(x=.5,y=1.5,label=paste('origin pv:',origin)),hjust=0,size=3)+
  scale_color_manual(values=ori.col[levels(TCA.df.div.m$origin)])+theme_bw()+guides(color = guide_legend(nrow = 1))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5),legend.position = 'top',legend.direction = 'horizontal',panel.spacing = unit(0,'lines'),axis.title.x = element_blank(),strip.background = element_rect(fill = NA))+
  facet_grid(~feature)+ylab("Normalized frational enrichment\n(divided by tumor GlcM6)")+ggtitle("Difference between patient tumor and PDX P2")
dev.off()

