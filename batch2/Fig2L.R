rm(list=ls());gc()
library(data.table);library(patchwork);library(dendextend);library(ComplexHeatmap);library(RColorBrewer);library(ggplot2);library(ggpubr);library(nlme)
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper/')
setwd('GitHub_repo/batch2/')
library(limma);library(openxlsx)
dat.list<-readRDS('cleaned_data_refined.rds')
metab<-dat.list$data
# metab<-metab[,grep('hydroxygl',colnames(metab),ignore.case = T,invert = T)]

s.info<-dat.list$s.info
s.info[,in_mouse:=as.numeric(Passage!=0)]

pv_df.mixedLinear<-data.table(metabolite=colnames(metab),N_zeros=colSums(metab==0),
                    t(data.table(s.info,metab)[,lapply(.SD,function(x) tryCatch({
                      tmp<-data.table(x,in_mouse,Tumor.Name)[!is.na(x)]
                      fit<-lme(x ~ in_mouse, random=~1|Tumor.Name,data = tmp)
                      c(fit$coefficients$fixed[['in_mouse']],anova(fit)[c('in_mouse'),'p-value'])
                    },error=function(e) rep(as.numeric(NA),2))),.SDcols=-(1:ncol(s.info))]))
colnames(pv_df.mixedLinear)[3:4]<-c('beta','pv')
pv_df.mixedLinear[,padj:=p.adjust(pv,'BH')]
pv_df.mixedLinear[padj<.05]
stat.1<-pv_df.mixedLinear
##################### 
keep.ind<-s.info[,.I[!Tumor.Name%in%c('MP4A','MP5')]]
dat<-metab[keep.ind,];s.info<-s.info[keep.ind]
sel.metab<-stat.1[signif(pv,2)>.1][['metabolite']]
dat.sub<-dat[,sel.metab]
dat.sub[]<-apply(dat.sub,2,scale)
dist_df<-as.data.table(melt(as.matrix(dist(dat.sub))))[Var1!=Var2]

# retain shortest distance pairs
tmp1<-s.info[match(dist_df$Var1,Mouse.ID),c('Tumor.Name','Passage')];colnames(tmp1)<-paste0('Var1.',colnames(tmp1))
tmp2<-s.info[match(dist_df$Var2,Mouse.ID),c('Tumor.Name','Passage')];colnames(tmp2)<-paste0('Var2.',colnames(tmp2))
dist_df<-cbind(dist_df,tmp1,tmp2);dist_df<-dist_df[,sort(colnames(dist_df)),with=F]
dist_df<-dist_df[Var1.Passage!=0&Var2.Passage==0]
close_df<-dist_df[dist_df[,.I[value==min(value)],by=Var1][['V1']]]
close_df[,PDX.passage:=paste(Var1.Tumor.Name,Var1.Passage)]
tmp<-unique(close_df$PDX.passage)
pdx_levels<-tmp[order(as.numeric(sapply(gsub("[A-Z]",'',tmp),function(x) unlist(strsplit(x,split = ' ',fixed = T))[[1]])),
                      as.numeric(sapply(gsub("[A-Z]",'',tmp),function(x) unlist(strsplit(x,split = ' ',fixed = T))[[2]])))]
close_df[,PDX.passage:=factor(PDX.passage,levels=pdx_levels)]
close_df[,Var2.Tumor.Name:=factor(Var2.Tumor.Name,levels=unique(close_df$Var2.Tumor.Name)[order(as.numeric(gsub("[A-Z]",'',unique(close_df$Var2.Tumor.Name))))])]
close_df[,matched_to_origin:=Var1.Tumor.Name==Var2.Tumor.Name]
close_df[,patient1:=gsub("[A-Z]$",'',Var1.Tumor.Name)];close_df[,patient2:=gsub("[A-Z]$",'',Var2.Tumor.Name)]
close_df[,matched_to_patient:=patient1==patient2]
close_df[,patient2:=factor(patient2,levels=unique(close_df$patient2)[order(as.numeric(gsub("[A-Z]",'',unique(close_df$patient2))))])]
close_df
by<-'patient'
pdf('2L_PDX_patient_match.pdf',w=9.5,h=2.2)
ggplot(close_df,aes(y=patient2,x=PDX.passage,color=matched_to_patient))+geom_point(alpha=.3,size=3)+theme_bw()+ylab('Patient')+xlab('PDX passage')+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=.5))+labs(color=paste('matched to',by))+scale_color_manual(values=c('TRUE'='seagreen','FALSE'='orange'))+
  ggtitle(paste0('PDX - patient match by ',length(sel.metab),' metabolites that do not differ between patient and PDX, ',nrow(close_df[patient1==patient2]),'/',nrow(close_df),' matched'))
dev.off()
