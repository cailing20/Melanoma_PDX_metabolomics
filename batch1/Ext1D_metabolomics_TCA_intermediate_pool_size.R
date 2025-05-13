rm(list=ls());gc()
library(tidyr);library(dplyr);library(ggplot2);library(data.table)
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper')
setwd('GitHub_repo/batch1/')
TCA.df<-fread('../tracing/S1_input.csv')
TCA.df<-TCA.df[sample!='Patient plasma'][grep("P6|P3",source,invert = T)][grepl("MP10|MP4A|MP8|MP9",source)][source!='MP9_Tumor']

dat1<-readRDS('cleaned_data.rds')
sort(grep("Pyruvate|lact|Phosphoenolpyruvate|Glyceraldehyde|alanine|Citrate|Glutamate|glutar|conita|Succinate|Fumarate|Malate|aspart",ignore.case = T,rownames(dat1$data),value = T))
select.metab<-c("Glyceraldehyde 3-phosphate", "Phosphoenolpyruvate", "Pyruvate", "Lactate", "Alanine","Citrate/Isocitrate", "alpha-Ketoglutarate", "Glutamate", "Succinate", "Fumarate", "Malate", "Aspartate")
paste(unique(sapply(unique(TCA.df$source),function(x) unlist(strsplit(x,split = '_'))[1])),collapse = '|')
s.info.sub<-dat1$s.info[grep("MP4A|MP8A|MP8B|MP9D|MP9F|MP10",toupper(Sample))]
dat1.sub<-dat1$data[select.metab,s.info.sub$Sample]
colnames(dat1.sub)<-paste(toupper(sapply(s.info.sub$Sample,function(x) unlist(strsplit(x,split = ' ',fixed = T))[1])),s.info.sub$`Human/PDX`,sep = '_')
m.df1<-as.data.table(melt(dat1.sub))

m.df1%>%separate(Var2,into=c('Origin',"Sample"),sep = '_')%>%setDT()->m.df1
m.df1[Origin=='MP10A',Origin:='MP10']
m.df1.h<-m.df1[Sample=='Human',-3,with=F];m.df1.m<-m.df1[Sample=='PDX',-3,with=F]
m.df1.cbn<-data.table(m.df1.h,m.df1.m$value[match(with(m.df1.h,paste(Var1,Origin)),with(m.df1.m,paste(Var1,Origin)))])
colnames(m.df1.cbn)[3:4]<-c('human','pdx')
m.df1<-m.df1[!Var1%in%c('Glyceraldehyde 3-phosphate',"alpha-Ketoglutarate")]
m.df1.cbn<-m.df1.cbn[!Var1%in%c('Glyceraldehyde 3-phosphate',"alpha-Ketoglutarate")]

pdf('Ext1D.pdf',w=10.5,h=2.5)
ggplot()+geom_point(m.df1,mapping=aes(x=Sample,y=value,color=Origin,group=Origin))+
  geom_line(m.df1,mapping=aes(x=Sample,y=value,color=Origin,group=Origin))+
  geom_text(m.df1.cbn[,paste('pv =',signif(t.test(human,pdx,paired=T)$p.value,2)),by=Var1],mapping=aes(x=1,y=m.df1[, .(min_y = min(value)), by = Var1][[2]],label=V1),vjust=0)+
  facet_wrap(~Var1,nrow = 2,scales = 'free_y')+theme_bw()
dev.off()
