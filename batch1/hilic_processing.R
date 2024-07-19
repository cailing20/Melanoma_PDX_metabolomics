rm(list=ls());gc()
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper/')
hilic.file<-'data/batch1/WORKING_Melanoma PDXHuman Global For Ling 20230202 v3.xlsx'
s.info<-read.xlsx(hilic.file,sheet = 3)
sheets(loadWorkbook(hilic.file))
# [1] "Updates"            "Introduction"       "SampleAnnotations"  "RawHILIC"           "TICAreas"           "Samples to include"
hilic<-read.xlsx(hilic.file,sheet = 4)
colnames(hilic)<-gsub('.',' ',colnames(hilic),fixed = T)
setDT(s.info)
# correction by Jennifer 5/22
s.info[Patient=='MP10',Pigment.Status:='None']
s.info[Patient%in%c('M481','M487','M405'),Sex:='Male']
s.info[Patient%in%c('M481','M487','M405')]
s.info[Sex==0]
clean.dat<-function(input.dat){
  input.dat<-input.dat[,grep("Blank|Pool",colnames(input.dat),invert = T,ignore.case = T)]
  updated.names<-gsub('M487 AB2641 (Non-Pigmented)','M487 AB2641 (NonPigmented)',gsub('TX-','TX',gsub(' -',' ',gsub(' - ',' ',colnames(input.dat)[-1],fixed = T),fixed = T),fixed = T),fixed = T)
  metabolites<-input.dat$X1
  tmp<-as.matrix(input.dat[,-1]);class(tmp)<-'numeric'
  rownames(tmp)<-metabolites
  zero.num<-rowSums(tmp==0)
  hist(zero.num) # decide to remove metabolites not detected in 33% of the samples
  tmp<-tmp[which(zero.num<ncol(tmp)*.9),]
  zero.num<-rowSums(tmp==0)
  tmp.sub<-tmp[which(zero.num==0 & log10(apply(tmp,1,median))>5 & log10(apply(tmp,1,median))<10),]
  norm.f<-apply(apply(tmp.sub,1,function(x) x/median(x[x!=0])),1,function(y) median(y[y!=0]))
  norm.df<-data.table(sample=names(norm.f),source=substr(names(norm.f),1,2),norm.f)
  plot(colSums(tmp),norm.f)
  plot(density(as.numeric(data.table(norm.f,t(tmp))[,lapply(.SD,function(x) cor.test(norm.f,x)$estimate),.SDcols=-1])))
  tmp2<-tmp;tmp2<-log10(tmp2);tmp2[is.infinite(tmp2)]<-0;colnames(tmp2)<-updated.names
  tmp2<-tmp;tmp2[]<-t(apply(tmp,1,function(x) x/norm.f));tmp2<-log10(tmp2);
  # impute
  tmp2[]<-t(apply(tmp2,1,function(x){x[is.infinite(x)]<-min(x[!is.infinite(x)]);x}))
  colnames(tmp2)<-updated.names
  tmp2
}

hilic<-clean.dat(hilic)
saveRDS(list(data=hilic,s.info=s.info),'work/batch1/cleaned_data.rds')