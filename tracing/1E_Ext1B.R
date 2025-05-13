rm(list=ls());gc()
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper/')
setwd('GitHub_repo/tracing/')
TCA.df<-fread('S1_input.csv')
TCA.df<-TCA.df[sample!='Patient plasma'][grep("P6|P3",source,invert = T)][grepl("MP10|MP4A|MP8|MP9",source)][source!='MP9_Tumor'][!grepl("ALANINE|ASPARTATE",metabolite)]
TCA.m<-c('GLUCOSE_M6','3PG_M3','PYRUVATE_M3','LACTATE_M3','CITRATE_M2','GLUTAMATE_M2','SUCCINATE_M2','FUMARATE_M2','MALATE_M2')
ori.col<-readRDS('..//batch2/ori_col.rds')
TCA.df[,sample_name:=gsub("[1-3]$",'',as.character(sample_name))]
TCA.df[,.N,by=.(source,sample_name)]
TCA.df<-TCA.df[,median(value),by=.(source,feature,metabolite,labeling,sample_name,sample)]
TCA.df[,value:=V1][,V1:=NULL]
TCA.df[,origin:=sapply(source,function(x) unlist(strsplit(x,split = '_'))[1])]
TCA.df[,type:=sample][,sample:=NULL]
library(dplyr);library(ggplot2)

TCA.df.norm.div<-TCA.df[feature!='GLUCOSE_M6']
ref<-TCA.df[feature=='GLUCOSE_M6']
for(f in unique(TCA.df.norm.div$feature)) TCA.df.norm.div[feature==f,value:=value/ref$value]

TCA.df.norm.div[,feature:=factor(feature,levels = TCA.m[-1])]
tmp<-str_to_title(paste0(gsub('_M',' (M+',levels(TCA.df.norm.div[['feature']])),')'))
tmp<-gsub('3pg','3PG',tmp)
levels(TCA.df.norm.div$feature)<-tmp

TCA.df.div.m<-TCA.df.norm.div[,median(value),by=.(source,feature,origin,type)]
TCA.df.div.m[,origin:=factor(origin,levels=intersect(names(ori.col),origin))]
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



# t-test for paired difference
TCA.df.div.m.pt<-TCA.df.div.m[type=='Patient tumor']
TCA.df.div.m.pdx<-TCA.df.div.m[type=='PDX']
TCA.df.div.m.pt[,diff:=V1-TCA.df.div.m.pdx[match(with(TCA.df.div.m.pt,paste(origin,feature)),with(TCA.df.div.m.pdx,paste(origin,feature)))][['V1']]]
stat.df<-TCA.df.norm.div[,as.list(structure(signif(as.matrix(summary(aov(value~origin*type))[[1]])[1:3,5],2),names=c('origin','type','origin:type'))),by=feature]
stat.df<-TCA.df.norm.div[,as.list(structure(signif(as.matrix(summary(aov(value~origin+type))[[1]])[1:2,5],2),names=c('origin','type'))),by=feature]
stat.df[,sig:='ns']
stat.df[type<.05,sig:='*']
stat.df[type<.005,sig:='**']
stat.df[type<.001,sig:='***']

pdf('Ext1B.pdf',w=12,h=3.5)
ggplot()+geom_point(TCA.df.div.m,mapping=aes(x=type,y=V1,color=origin))+
  geom_line(TCA.df.div.m,mapping=aes(x=type,y=V1,group=origin,color=origin))+
  geom_errorbar(TCA_summary[feature!='GLUCOSE_M6'],mapping=aes(x=type,ymin = lower_quartile, ymax = upper_quartile,color=origin), width = 0.2)+
  geom_text(stat.df,mapping=aes(x=.5,y=1.3,label=paste('type pv:',type)),hjust=0,size=3)+
  geom_text(stat.df,mapping=aes(x=.5,y=1.5,label=paste('origin pv:',origin)),hjust=0,size=3)+
  scale_color_manual(values=ori.col[levels(TCA.df.div.m$origin)])+theme_bw()+guides(color = guide_legend(nrow = 1))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5),legend.position = 'top',legend.direction = 'horizontal',panel.spacing = unit(0,'lines'),axis.title.x = element_blank(),strip.background = element_rect(fill = NA))+
  facet_grid(~feature)+ylab("Normalized frational enrichment\n(divided by tumor GlcM6)")#+ggtitle("Difference between patient tumor and PDX P2")
dev.off()

# stats
TCA.df.norm.div.tumor<-TCA.df.norm.div[type=='Patient tumor']
TCA.df.norm.div.tumor[,lab.tumor:=paste(origin,feature)]
TCA.df.norm.div.pdx<-TCA.df.norm.div[type=='PDX']
TCA.df.norm.div.pdx[,lab.pdx:=paste(origin,feature)]
TCA.df.norm.div.pdx[,tumor.val:=TCA.df.norm.div.tumor[match(lab.pdx,lab.tumor)][['value']]]
TCA.df.norm.div.pdx[,diff:=value-tumor.val]
TCA.df.norm.div.pdx[,origin:=factor(origin,levels=intersect(names(ori.col),origin))]
pdf('1E.pdf',w=4,h=4.5)
ggplot(TCA.df.norm.div.pdx[!grepl("Alanine|Aspartate",as.character(feature))],aes(x=feature,y=diff))+
  geom_text(stat.df[!grepl("Alanine|Aspartate",as.character(feature))],mapping=aes(x=feature,y=2,label=sig))+
  geom_beeswarm(TCA.df.norm.div.pdx[!grepl("Alanine|Aspartate",as.character(feature))],mapping=aes(color=origin),alpha=.5)+
  stat_summary(fun = "mean",geom = "crossbar", width = 0.5,size=.2)+
  geom_hline(yintercept = 0,linetype='dashed')+ylab("PDX- Pt Fractional Enrichment\n(Normalised to Glucose M+6)")+
  theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size=11),axis.title.x = element_blank(),legend.position = 'bottom',axis.line = element_line(linewidth = .2))+
  guides(color = guide_legend(nrow = 2))
dev.off()

# n statement
TCA.df.norm.div[type=='Patient tumor',.N,by=sample_name] # 6 pt tumor
TCA.df.norm.div[type=='PDX',.N,by=sample_name] # 6 pt tumor
stat.df[,feature:=factor(as.character(feature),levels=levels(TCA.df.div.m$feature))]
levels(stat.df$feature)<-tmp
TCA.df.norm.div.pdx[type=='PDX',.N,by=sample_name]
