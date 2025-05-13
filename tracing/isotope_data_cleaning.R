rm(list=ls());gc()
library(openxlsx);library(data.table);library(ggplot2)
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper/')
x.file<-'data/tracing/20240601 Stable Isotope Tracing Data for Paper.xlsx'
ss<-sheets(loadWorkbook(x.file))
dat.list<-structure(lapply(1:length(ss),function(i) read.xlsx(x.file,sheet=i)),names=ss)
lapply(dat.list, colnames)
for(i in which(unlist(lapply(dat.list, function(x) any(colnames(x)=='3h'))))){
  x<-dat.list[[i]]
  colnames(x)[-(1:2)]<-paste(as.character(x[1,-(1:2)]),substr(colnames(x)[-(1:2)],1,2),sep='_')
  x<-x[-1,]
  dat.list[[i]]<-x
}
f.df<-do.call(rbind,lapply(ss, function(s){
  x<-dat.list[[s]]
  x<-x[,1:2];colnames(x)[1:2]<-c('c1','c2')
  data.table(source=s,x)
}))
unique(f.df[grep(' ',c1,fixed = T)][[2]])
unique(f.df[grep(' ',c1,fixed = T,invert = T)][[2]])
f.df[,c1:=gsub('Isoleucine(','Isoleucine (',c1,fixed = T)]
f.df[,o:=1:nrow(f.df)]
f.df[,nf:=sapply(c1,function(x) toupper(unlist(strsplit(x,split = ' ',fixed = T))[1]))]
sort(unique(f.df[['nf']]))
f.df[grep('DHAP',nf)]
f.df[nf%in%c("3-PGLYCERATE","3-PHOSPHOGLYCERATE","3P-GLYCERATE"),nf:='3PG']
f.df[grep("IC$",nf),nf:=gsub('IC','ATE',nf)]
f.df[nf%in%c('DHAP#1'),nf:='DHAP']
f.df[nf%in%c('TYRSOINE'),nf:='TYROSINE']
f.df[,r:=rank(o)-1,by=.(source,nf)]
f.df[,feature:=paste0(nf,'_M',r)]
unique(f.df[['feature']])
for(s in ss){
  x<-dat.list[[s]]
  y<-setdiff(as.numeric(f.df[source==s][['c2']]),c(NA,NaN))
  if(sum(sum(y-floor(y)))>0) x<-x[,-1,drop=F] else x<-x[,-(1:2),drop=F]
  rownames(x)<-f.df[source==s][['feature']]
  dat.list[[s]]<-x
}
dat.list[]<-lapply(dat.list, function(x){x<-as.matrix(x);class(x)<-'numeric';as.data.frame(x)})
dat.df<-do.call(rbind,lapply(ss,function(s) data.table(s,melt(as.matrix(dat.list[[s]])))))
colnames(dat.df)<-c('source','feature','sample_name','value')
dat.df[,value:=as.numeric(value)]
library(dplyr);library(tidyr)
dat.df%>%separate(feature,into = c('metabolite','labeling'),sep = '_',remove = F)->dat.df;setDT(dat.df)
dat.df[grepl('Tumor',source),sample:='Patient tumor'][grepl('Plasma',source),sample:='Patient plasma'][grepl('PDX',source),sample:='PDX']
dat.df<-dat.df[!sample_name%in%dat.df[sample=='Patient plasma' & feature=='GLUCOSE_M6'][value<.1][['sample_name']]]
dat.df[sample_name=='9D.Tumour',source:='MP9D_Tumor'][sample_name=='9F.Tumour',source:='MP9F_Tumor']
dat.df[sample_name=='8A.Tumour',source:='MP8A_Tumor'][sample_name=='8B.Tumour',source:='MP8B_Tumor']
dat.df[sample_name=='MP4A.Tumour1_1004004',source:='MP4A_Tumor'][grepl('MP4a',source),source:=gsub('MP4a','MP4A',source)]
m.ref<-f.df[!duplicated(paste(source,nf)),.N,by=nf][order(-N)][N>10]
for(i in 1:nrow(m.ref)){
  pdf(paste0('work/isotope/',m.ref$nf[[i]],'.pdf'),w=6,h=m.ref$N[[i]]*.7)
  print(ggplot(dat.df[metabolite==m.ref$nf[[i]]],aes(y=source,x=value))+geom_point()+facet_grid(labeling+sample~.,space='free',scales = 'free_y')+theme_bw())
  dev.off()
}
ggplot(dat.df[metabolite=='PYRUVATE'],aes(y=source,x=value))+geom_point()+facet_grid(labeling+sample~.,space='free',scales = 'free_y')+theme_bw()
ggplot(dat.df[metabolite=='CITRATE'],aes(y=source,x=value))+geom_point()+facet_grid(labeling+sample~.,space='free',scales = 'free_y')+theme_bw()
ggplot(dat.df[metabolite=='MALATE'],aes(y=source,x=value))+geom_point()+facet_grid(labeling+sample~.,space='free',scales = 'free_y')+theme_bw()
ggplot(dat.df[metabolite=='GLUCOSE'],aes(y=source,x=value))+geom_point()+facet_grid(labeling+sample~.,space='free',scales = 'free_y')+theme_bw()
ggplot(dat.df[metabolite=='UREA'],aes(y=source,x=value))+geom_point()+facet_grid(labeling+sample~.,space='free',scales = 'free_y')+theme_bw()



TCA.m<-c('GLUCOSE_M6','3PG_M3','PYRUVATE_M3','LACTATE_M3','ALANINE_M3','CITRATE_M2','GLUTAMATE_M2','SUCCINATE_M2','FUMARATE_M2','MALATE_M2','ASPARTATE_M2')
# do not collapse to median
TCA.with_rep.df<-dat.df[feature%in%TCA.m&grepl("MP4A|MP10|MP8A|MP8B|MP9D|MP9F",source)]
saveRDS(TCA.with_rep.df,'work/isotope/tca_df_with_rep.rds')
# collapse to median
# calculate correlations
m.df<-dat.df[,median(value,na.rm = T),by=.(source,feature,metabolite,labeling,sample)]
library(reshape2)
m.mat<-acast(m.df[metabolite%in%m.ref[['nf']]],source~feature,value.var = 'V1',fill = NA)
library(Hmisc)
r.df<-data.table(melt(rcorr(m.mat)$r),melt(rcorr(m.mat)$P)[,3])
colnames(r.df)[3:4]<-c('r','pv')
r.df<-r.df[!is.na(pv)]
r.df[,m1:=sapply(as.character(Var1),function(x) unlist(strsplit(x,split = '_'))[1])][,m2:=sapply(as.character(Var2),function(x) unlist(strsplit(x,split = '_'))[1])]
r.df[m1!=m2][pv<.05][Var1=='LACTATE_M3']
n.mat<-m.mat;n.mat[]<-apply(m.mat,2,function(x) x/m.mat[,'GLUCOSE_M6'])
ggplot(dat.df[feature=='GLUCOSE_M6'][sample=='Patient tumor'],aes(x=value,y=source))+geom_point()
ggplot(dat.df[feature=='GLUCOSE_M6'][sample=='Patient plasma'],aes(x=value,y=source))+geom_point()
plot(m.df[feature=='GLUCOSE_M6'][sample=='Patient tumor'][['V1']],m.df[feature=='GLUCOSE_M6'][sample=='Patient plasma'][['V1']],xlim = c(0,.5),ylim = c(0,.5))
abline(a=0,b=1)

tmp.df<-dat.df[grepl("PDX_P|Tumor",source)&labeling=='M0']
tmp.df[,labeled.frac:=1-value]
pdf('labeled_fractions.pdf',w=12,h=9)
ggplot(tmp.df[metabolite%in%m.ref$nf],aes(y=source,x=labeled.frac))+geom_point()+facet_wrap(~metabolite,scales = 'free_x',nrow = 4)+theme_bw()+geom_vline(xintercept = 0,color='navy')
dev.off()

tca.df<-data.table(sample=rownames(m.mat),type=dat.df[match(rownames(m.mat),source)][['sample']],m.mat[,TCA.m])
library(ggrepel)
g.list<-lapply(1:(length(TCA.m)-1),function(i){
  ggplot(tca.df,aes_string(x=paste0('`',TCA.m[i],'`'),y=paste0('`',TCA.m[i+1],'`')))+geom_point(mapping=aes(color=type))+theme_bw()+theme(legend.position = 'none')+
    geom_point(tca.df[sample%in%c('MP9D_PDX_P6','MP10_PDX_P6')],mapping=aes_string(x=paste0('`',TCA.m[i],'`'),y=paste0('`',TCA.m[i+1],'`')),color='black')+
    geom_text_repel(tca.df[sample%in%c('MP9D_PDX_P6','MP10_PDX_P6')],mapping=aes(label=sample))+geom_abline(slope = 1,intercept = 0)+geom_smooth(method = 'lm',alpha = .2)
})
g.list[[4]]<-g.list[[4]]+theme(legend.position = 'right')
library(patchwork)
pdf('scatter_outliers_highlited.pdf',w=12,h=7);wrap_plots(g.list,nrow = 3);dev.off()

library(ggbeeswarm)
# reproduce 1D
for(m in TCA.m){
  print(m)
  print(identical(dat.df[feature=='PYRUVATE_M3' & sample!='Patient plasma'][['sample_name']],
                  dat.df[feature==m & sample!='Patient plasma'][['sample_name']]))
}
TCA.mat<-do.call(cbind,lapply(TCA.m,function(m){
  dat.df[feature==m & sample!='Patient plasma'][['value']]
}))
colnames(TCA.mat)<-TCA.m
tmp.df<-data.table(dat.df[feature=='PYRUVATE_M3' & sample!='Patient plasma',-match('value',colnames(dat.df)),with=F][,-2,with=F],TCA.mat)
tmp.df[,origin:=sapply(source,function(x) unlist(strsplit(x,split = '_'))[1])]
tmp.df<-tmp.df[grep("MP4A|MP10|MP8A|MP8B|MP9D|MP9F",origin)]#[GLUCOSE_M6>.2]
tmp.df[,source:=gsub('_3h','',source)];tmp.df[source%in%c('MP4A_PDX_P2','MP4A_PDX_P3'),source:='MP4A_PDX_P2/3']
tmp.df[,source:=factor(source,levels = unique(tmp.df$source))]
tmp.df[,sample_name2:=sapply(gsub("NB",'@NB',gsub("AB",'@AB',sample_name)),function(x) substr(unlist(strsplit(x,split = '@'))[2],1,6))]
tmp.df[is.na(sample_name2),sample_name2:=source]
tmp.df<-tmp.df[,lapply(.SD,median),by=.(source,sample_name2,origin,sample),.SDcols = match(TCA.m,colnames(tmp.df))]
# set.seed(7464);Heatmap(tmp.df[,TCA.m,with=F],right_annotation = HeatmapAnnotation(which='row',df=tmp.df[,c('source','origin','sample')]))

tmp.df.sub<-tmp.df[grep("MP4A|MP10|MP8A|MP9D",source)]
polish.g<-function(g) g+geom_beeswarm(alpha=.3)+geom_smooth(method='lm',aes(group=origin))+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=.5))
g1<-polish.g(ggplot(tmp.df.sub,aes(x=source,y=GLUCOSE_M6)))
g2<-polish.g(ggplot(tmp.df.sub,aes(x=source,y=PYRUVATE_M3)))
g3<-polish.g(ggplot(tmp.df.sub,aes(x=source,y=LACTATE_M3)))
g4<-polish.g(ggplot(tmp.df.sub,aes(x=source,y=PYRUVATE_M3/GLUCOSE_M6)))
g5<-polish.g(ggplot(tmp.df.sub,aes(x=source,y=LACTATE_M3/GLUCOSE_M6)))
g6<-polish.g(ggplot(tmp.df.sub,aes(x=source,y=LACTATE_M3/PYRUVATE_M3)))
g7<-polish.g(ggplot(tmp.df.sub,aes(x=source,y=CITRATE_M2/PYRUVATE_M3)))
(g1|g2|g3)/(g4|g5|g6)
g7
fwrite(tmp.df,'work/isotope/1D_test_input.csv')
# 'GLUCOSE_M6','3PG_M3','PYRUVATE_M3','LACTATE_M3','ALANINE_M3','CITRATE_M2','GLUTAMATE_M2','SUCCINATE_M2','FUMARATE_M2','MALATE_M2','ASPARTATE_M2'
polish.g(ggplot(tmp.df.sub,aes(x=source,y=GLUTAMATE_M2/CITRATE_M2)))
polish.g(ggplot(tmp.df.sub,aes(x=source,y=SUCCINATE_M2/GLUTAMATE_M2)))
polish.g(ggplot(tmp.df.sub,aes(x=source,y=MALATE_M2/CITRATE_M2)))
polish.g(ggplot(tmp.df.sub,aes(x=source,y=ASPARTATE_M2/GLUTAMATE_M2)))



tmp.df.sub<-tmp.df[grep("MP4A|MP10|MP8|MP9",source)][!grepl('P6',source)]
# tmp.df.sub<-data.table(tmp.df.sub[,1:4,with=F],tmp.df.sub[,lapply(.SD,function(x) x/GLUCOSE_M6),.SDcols=-(1:4)])
tmp.df.sub<-data.table(tmp.df.sub[,1:4,with=F],tmp.df.sub[,lapply(.SD,function(x) x/PYRUVATE_M3),.SDcols=-(1:4)])
ggplot(melt(tmp.df.sub),aes(x=origin,y=value,alpha=sample,color=origin))+geom_boxplot()+facet_grid(~variable)+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=.5))

fwrite(dat.df,'work/isotope/dat_df.csv')
saveRDS(m.mat,'work/isotope/m_mat.rds')
saveRDS(tca.df,'work/isotope/tca_df.rds')
saveRDS(TCA.m,'work/isotope/TCA_metabolites.rds')
saveRDS(tmp.df,'work/isotope/merged_df.rds')
