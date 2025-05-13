rm(list=ls());gc()
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper/')
setwd('GitHub_repo/batch2/')
dat.list<-readRDS('cleaned_data_refined.rds')
s.info<-dat.list$s.info
s.info[,Tumor.Name:=trimws(Tumor.Name)]
s.info[,pt:=gsub("[A-Z]$",'',Tumor.Name)]
s.info[,pt.id:=as.numeric(gsub("[A-Z]",'',pt))]
s.info[,Tumor.Name:=factor(Tumor.Name,levels=s.info[!duplicated(Tumor.Name)][order(pt.id,Tumor.Name)][['Tumor.Name']])]
s.info[,pt:=factor(pt,levels=s.info[!duplicated(pt)][order(pt.id)][['pt']])]
s.count<-s.info[,.N,by=.(pt,Tumor.Name,Passage)]
library(viridis)
s.count[,Tumor.Name:=factor(Tumor.Name,levels=rev(unique(sort(Tumor.Name))))]
pdf('Ext5_sample_summary.pdf',w=2.3,h=3.3)
ggplot(s.count,aes(x=Passage,y=Tumor.Name,fill=N))+geom_tile(color='white')+geom_text(aes(label=N),color='white')+
  facet_grid(pt~.,space = 'free',scales = 'free')+scale_fill_viridis_c(begin = 0.2,end = 0.8)+theme_bw()+
  theme(panel.spacing = unit(0,unit='lines'),strip.background = element_rect(fill = NA),axis.title.y = element_blank(),
        panel.grid = element_blank(),legend.position = 'none',strip.text.y = element_text(angle = 0))+
  scale_x_continuous(breaks = 0:6)
dev.off()
