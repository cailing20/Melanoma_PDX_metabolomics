rm(list=ls());gc()
lapply(c('data.table','openxlsx','ggplot2','patchwork'),function(x) library(x,character.only = T))
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper/')
sheets(loadWorkbook('data/batch2/051424_Human PDX Metabolomics_for Ling_POSTGENEPRINT.xlsx'))
# [1] "Updates"                "Introduction"           "OLD_Sample Annotations" "NEW_Sample Annotations" "Sample Order Run"       "RawHILIC_Samples"       "TICAreas_Samples"      
# [8] "RawHILIC_QC POOL"       "TICAreas_QC POOL"  
s.info<-as.data.table(read.xlsx('data/batch2/051424_Human PDX Metabolomics_for Ling_POSTGENEPRINT.xlsx',3)[-1,])
# there is a sample misannotation problem, we need to take care of this after the signal correction
s.info.corrected<-as.data.table(read.xlsx('data/batch2/051424_Human PDX Metabolomics_for Ling_POSTGENEPRINT.xlsx',4)[-1,])
s.order<-as.data.table(read.xlsx('data/batch2/051424_Human PDX Metabolomics_for Ling_POSTGENEPRINT.xlsx',5))
s.dat<-read.xlsx('data/batch2/051424_Human PDX Metabolomics_for Ling_POSTGENEPRINT.xlsx',6)
qc.dat<-read.xlsx('data/batch2/051424_Human PDX Metabolomics_for Ling_POSTGENEPRINT.xlsx',8)

tmp<-rbind(data.table(source='original',s.info[,1:4,with=F]),data.table(source='corrected',s.info.corrected))
tmp<-tmp[paste(Tube,Tumor.Name,Passage,Mouse.ID)%in%tmp[,.N,by=paste(Tube,Tumor.Name,Passage,Mouse.ID)][N==1][['paste']]]
tmp
# what is unique in original and in corrected:
# source   Tube Tumor.Name Passage Mouse.ID
# <char> <char>     <char>  <char>   <char>
#  1:  original     26       TX36       6   AB4108
#  2:  original     34       TX16       4   AB3764
#  3:  original     36       TX16       3   AB3073
#  4:  original     59       TX16       0  TX16Hum
#  5:  original     64       TX36       6   AB4107
#  6:  original     88       TX16       6   AB6674
#  7:  original    120       TX36       6   AB4109
#  8:  original    123       TX16       2   AB2328
#  9:  original    124       TX16       5   AB4100
# 10:  original    125       TX16       1   NB4126
# 11:  original    140       TX16       5   AB4099
# 12:  original    141       TX16       6   AB6672
# 13:  original    145       TX16       2   AB2326
# 14:  original    165       TX16       6   AB6671
# 15:  original    171       TX16       2   AB2327
# 16:  original    176       TX16       5   AB4096
# 17:  original    180       TX16       3   AB3074
# 18:  original    191       TX16       3   AB3072
# 19: corrected    124       TX36       6        A
# 20: corrected    140       TX36       6        B
# 21: corrected    176       TX36       6        C

# basically all TX16 are removed, except 3 that will be renamed to TX36, passage 6
tmp[Tube%in%tmp[,.N,by=Tube][N==2][['Tube']]]
# source   Tube Tumor.Name Passage Mouse.ID
# <char> <char>     <char>  <char>   <char>
# 1:  original    124       TX16       5   AB4100
# 2:  original    140       TX16       5   AB4099
# 3:  original    176       TX16       5   AB4096
# 4: corrected    124       TX36       6        A
# 5: corrected    140       TX36       6        B
# 6: corrected    176       TX36       6        C

# Let's save this for use in later steps
fwrite(tmp,'work/PDX_sample_update_reference.csv')

qc.dat[1:3,1:3]
colnames(s.order)[1]<-'sample'
s.order[,order:=1:nrow(s.order)]
s.order[,in_run:=sample%in%colnames(cbind(s.dat,qc.dat))]

# Step 1. construct a data to be in the same order as the run order
setdiff(c(colnames(s.dat),colnames(qc.dat)),s.order[,1])
identical(s.dat[,1],qc.dat[,1])
dat<-t(cbind(s.dat,qc.dat)[,s.order[which(in_run)][['sample']]])
colnames(dat)<-s.dat[,1]
dat[1:4,1:4]
table(colSums(dat)==0) # 209 all zeros
hist(colSums(dat==0))
dat<-dat[,colSums(dat)!=0]

pool.ind<-s.order[which(in_run)][,.I[grep("Pool",sample)]];p.tmp<-data.table(ro=s.order[which(in_run)][['order']][pool.ind],dat[pool.ind,])
ro.stats1<-t(p.tmp[,lapply(.SD,function(x) unlist(cor.test(x,ro)[c('estimate','p.value')])),.SDcols=-1])

qc.ind<-s.order[which(in_run)][,.I[grep("QC",sample)]];q.tmp<-data.table(ro=s.order[which(in_run)][['order']][qc.ind],dat[qc.ind,])
ro.stats2<-t(q.tmp[,lapply(.SD,function(x) unlist(cor.test(x,ro)[c('estimate','p.value')])),.SDcols=-1])

ro.stats<-merge(ro.stats1,ro.stats2,by=0)
colnames(ro.stats)<-c('metabolite','pool.r','pool.pv','QC.r','QC.pv')
ggplot(ro.stats,aes(x=pool.r,QC.r))+geom_point()

s.ind<-s.order[which(in_run)][,.I[grep("QC|Pool",sample,invert = T)]];s.tmp<-data.table(ro=s.order[which(in_run)][['order']][s.ind],dat[s.ind,])
ro.stats3<-t(s.tmp[,lapply(.SD,function(x) unlist(cor.test(x,ro)[c('estimate','p.value')])),.SDcols=-1])
colnames(ro.stats3)<-c('sample.r','sample.pv')
ro.stats<-data.table(ro.stats,ro.stats3)
with(ro.stats,cor.test(pool.r,sample.r))
ggplot(ro.stats,aes(x=pool.r,y=sample.r))+geom_point()
ggplot(ro.stats,aes(x=QC.r,y=sample.r))+geom_point()
g1<-ggplot(p.tmp,aes(x=ro,y=Indole))+geom_point()
g2<-ggplot(s.tmp,aes(x=ro,y=Indole))+geom_point()
g3<-ggplot(q.tmp,aes(x=ro,y=Indole))+geom_point()
source('/project/CRI/DeBerardinis_lab/lcai/CommonFunctions/beautify_df.R')
ro.stats<-beautify.dt(ro.stats)
# fwrite(ro.stats,'../work/June2023/run_order_correlation.csv')

g1+g2+g3
library(ComplexHeatmap)
dat.zero<-dat;dat.zero[dat!=0]<-1
Heatmap(dat.zero,row_names_gp = gpar(fontsize=5),name='non-zero value')
dat.zero.sub<-dat.zero[grep("Pool|QC",rownames(dat.zero),invert = T),]
dat.zero.sub<-data.table(ro=s.order[match(rownames(dat.zero.sub),sample)][['order']],po=as.numeric(s.info[match(rownames(dat.zero.sub),Tube)][['Passage']]),dat.zero.sub)
ggplot(dat.zero.sub,aes(x=po,y=ro))+geom_point()
po_pv<-t(dat.zero.sub[,lapply(.SD,function(x) tryCatch({t.test(po~x)$p.value},error=function(e) as.numeric(NA))),.SDcols=-(1:2)])
po_pv<-data.table(metabolite=rownames(po_pv),pv=po_pv[,1])[!is.na(pv)][order(pv)];po_pv[,padj:=p.adjust(pv,'BH')]
po_pv2<-t(dat.zero.sub[,lapply(.SD,function(x) tryCatch({fisher.test(table((po==0),x))$p.value},error=function(e) as.numeric(NA))),.SDcols=-(1:2)])
po_pv2<-data.table(metabolite=rownames(po_pv2),pv=po_pv2[,1])[!is.na(pv)][order(pv)];po_pv2[,padj:=p.adjust(pv,'BH')]
ro_pv<-t(dat.zero.sub[,lapply(.SD,function(x) tryCatch({wilcox.test(ro~x)$p.value},error=function(e) as.numeric(NA))),.SDcols=-(1:2)])
ro_pv<-data.table(metabolite=rownames(ro_pv),pv=ro_pv[,1])[!is.na(pv)][order(pv)];ro_pv[,padj:=p.adjust(pv,'BH')]
# missing values by passage - how many are missing?
nrow(dat.zero.sub) # 214
colSums(dat.zero.sub==0)[-(1:2)][po_pv[padj<.05][['metabolite']]]
dat.sub<-data.table(dat.zero.sub[,1:2,with=F],dat[grep("Pool|QC",rownames(dat),invert = T),])
wrap_plots(lapply(po_pv[padj<.05][['metabolite']],function(m) ggplot(dat.sub,aes_string(x='po',y=paste0('`',m,'`')))+geom_jitter(width = .2,alpha=.3)))
wrap_plots(lapply(po_pv2[padj<.05][['metabolite']],function(m) ggplot(dat.sub,aes_string(x='po',y=paste0('`',m,'`')))+geom_jitter(width = .2,alpha=.3)))

wrap_plots(lapply(intersect(c(po_pv[pv<.05][['metabolite']],po_pv2[pv<.05][['metabolite']]),ro_pv[padj<.05][['metabolite']]),
                  function(m) ggplot(dat.sub,aes_string(x='po',y=paste0('`',m,'`')))+geom_jitter(width = .2,alpha=.3)))
saveRDS(ro_pv[padj<.05][['metabolite']],'work/ro_metab.rds')
Heatmap(dat.zero[grep("Pool|QC",rownames(dat.zero),invert = T),ro_pv[padj<.05][['metabolite']]],
        row_names_gp = gpar(fontsize=5),column_names_gp = gpar(fontsize=8),name='non-zero value',cluster_rows = F)

structure(lapply(grep('pv',colnames(ro.stats),value = T),function(x) length(which(ro.stats[[x]]<.05))),names=grep('pv',colnames(ro.stats),value = T))
