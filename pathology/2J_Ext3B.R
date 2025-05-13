setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper/')
setwd('GitHub_repo/pathology/')
library(openxlsx);library(data.table)
pt.dat<-read.xlsx('Pt PDX Histology Data 012525.xlsx',startRow = 2)
setDT(pt.dat)
pt.dat[,`%Viable.Tumor`:=as.numeric(`%Viable.Tumor`)]
pt.dat[,`%Ki-67`:=as.numeric(`%Ki-67`)]
pt.dat[,`%Immune.Cells.(and.Type)`:=as.numeric(`%Immune.Cells.(and.Type)`)]
pt.sub<-pt.dat[grep('Pt',X3)][,3:8,with=F]
pdx.sub<-pt.dat[grep('PDX',X3)][,3:8,with=F]
pt.sub[,X3:=gsub(' Pt','',X3)]
pdx.sub[,X3:=gsub(' PDX','',X3)]
cbind(with(melt(pt.sub),paste(X3,variable)),with(melt(pdx.sub),paste(X3,variable)))
m.df<-data.table(melt(pdx.sub),melt(pt.sub)[,3])
colnames(m.df)<-c('origin','feature','PDX','Pt')
library(ggrepel)
# manual correction:
m.df[origin=='MP4',origin:='MP4a'][origin=='MP2a',origin:='MP2b']
pdf('Ext3B_histology_compare.pdf',w=11,h=6)
ggplot(m.df[!is.na(PDX)][!is.na(Pt)],aes(x=Pt,y=PDX))+geom_point(alpha=.3,size=5,aes(color=origin))+
  geom_text_repel(m.df[paste(feature,PDX,Pt)%in%with(m.df[duplicated(paste(feature,PDX,Pt))],paste(feature,PDX,Pt))],mapping=aes(x=Pt,y=PDX,label=origin,color=origin),
                  force = 200,force_pull = 5,max.overlaps = Inf,show.legend = F)+geom_abline(slope = 1,intercept = 0)+
  facet_wrap(~feature,scales='free')+theme_bw()+stat_cor()#+xlim(0,100)+ylim(0,100)
dev.off()
tmp<-as.data.table(melt(m.df[!is.na(PDX)][!is.na(Pt)][feature=='%Ki-67'][,-2,with=F]))
tmp[,variable:=factor(variable,levels=c('Pt','PDX'))]
tmp[,origin:=factor(origin,levels=tmp[variable=='Pt'][order(value)][['origin']])][variable=='Pt',variable:='Human']
tmp[,variable:=factor(variable,levels=c('Human','PDX'))]
pdf('2J_histology_compare_Ki67_bar.pdf',w=2.8,h=3.1)
ggplot(tmp,aes(y=origin,x=value,fill=variable))+geom_bar(stat='identity',position = 'dodge',color='gray')+theme_bw()+xlab('%Ki-67')+scale_fill_manual(values=c(Human='#A50F15',PDX="#FEE5D9"))+
  ggtitle(paste0("Paired t-test pval = ",signif(with(m.df[!is.na(PDX)][!is.na(Pt)][feature=='%Ki-67'][,-2,with=F],t.test(x=PDX,y=Pt,paired=T))$p.value,2)))+theme(axis.title.y = element_blank())
dev.off()