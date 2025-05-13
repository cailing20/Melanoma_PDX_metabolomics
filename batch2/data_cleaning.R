rm(list=ls());gc()
lapply(c('data.table','openxlsx','ggplot2','patchwork'),function(x) library(x,character.only = T))
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper/')
input.data.file<-'data/batch2/051424_Human PDX Metabolomics_for Ling_POSTGENEPRINT.xlsx'
sheets(loadWorkbook(input.data.file))
# [1] "Updates"                "Introduction"           "OLD_Sample Annotations" "NEW_Sample Annotations" "Sample Order Run"       "RawHILIC_Samples"       "TICAreas_Samples"      
# [8] "RawHILIC_QC POOL"       "TICAreas_QC POOL"  
s.info<-as.data.table(read.xlsx(input.data.file,3)[-1,])
s.order<-as.data.table(read.xlsx(input.data.file,5))
s.dat<-read.xlsx(input.data.file,6)
qc.dat<-read.xlsx(input.data.file,8)
qc.dat[1:3,1:3]
colnames(s.order)[1]<-'sample'
s.order[,order:=1:nrow(s.order)]
s.order[,in_run:=sample%in%colnames(cbind(s.dat,qc.dat))]

# Step 1. construct data to be in the same order as the run order
setdiff(c(colnames(s.dat),colnames(qc.dat)),s.order[,1])
identical(s.dat[,1],qc.dat[,1]) # same metabolite order
dat<-t(cbind(s.dat,qc.dat)[,s.order[which(in_run)][['sample']]])
colnames(dat)<-s.dat[,1]
dat[1:4,1:4]
table(colSums(dat)==0) # 209 all zeros
hist(colSums(dat==0))
# remove metabolites with all missing value
dat<-dat[,colSums(dat)!=0]

# Step 2. take out QC data since it is not as informative as the Pooled sample
dat<-dat[grep("QC",rownames(dat),invert = T),]
s.order[grepl('QC',sample),in_run:=F]
identical(s.order[which(in_run)][[1]],rownames(dat))
# plot to view
pv_list<-structure(lapply(colnames(dat),function(x){
  tryCatch({  tmp<-data.table(ro=s.order[which(in_run)][['order']],metabolite=dat[,x])[metabolite!=0]
  lq <- quantile(tmp$metabolite, 0.25);uq <- quantile(tmp$metabolite, 0.75)
  iqr <- uq - lq
  tmp<-tmp[metabolite >= lq - 1.5 * iqr & metabolite <= uq + 1.5 * iqr]
  cor.test(tmp$ro,tmp$metabolite)$p.value},
  error=function(e){NA})
}),names=colnames(dat))
g_list<-lapply(colnames(dat),function(x){
  tryCatch({  tmp<-data.table(ro=s.order[which(in_run)][['order']],metabolite=dat[,x])[metabolite!=0]
  lq <- quantile(tmp$metabolite, 0.25);uq <- quantile(tmp$metabolite, 0.75)
  iqr <- uq - lq
  tmp<-tmp[metabolite >= lq - 1.5 * iqr & metabolite <= uq + 1.5 * iqr]
  ggplot(tmp,aes(x=ro,y=metabolite))+geom_point(alpha=.3,color='navy')+
    annotate('text',Inf,Inf,label=with(with(tmp,cor.test(ro,metabolite)),paste0('r=',round(estimate,2),',pv=',signif(p.value,2))),hjust=1.1,vjust=1.2)+
    theme_bw()+ggtitle(x)},
  error=function(e){NULL})
})
pv_list<-unlist(pv_list)[!unlist(lapply(g_list,is.null))]
g_list<-g_list[!unlist(lapply(g_list,is.null))]
# show all plots with significant signal correlation to run order
wrap_plots(g_list[pv_list<.05],ncol = 11)
# show a few heavily affected examples
wrap_plots(g_list[match(c('Carnitine (18:1)','N-Acetylserine','Sulfolactate'),names(pv_list))])

# exclude metabolites that have discontinous trend in missing values or low values
ro.metab<-c(readRDS('../work/June2023/ro_metab.rds'),c('Carnitine (18:1)','N-Acetylserine','Sulfolactate'))
# perform correction in metabolites affected by run order
attention.metab<-setdiff(names(pv_list)[pv_list<.05],ro.metab)

dim(dat)
# a function for extrapolating missing values

# Define a function to impute missing values
impute_missing <- function(df) {
  # Find the indices of the missing values
  missing_indices <- which(is.na(df$y))
  
  # Loop through the missing indices
  for (i in missing_indices) {
    # Find the neighboring indices
    neighbors <- (i - 5):(i + 5)
    # Remove out of range indices
    neighbors <- neighbors[neighbors > 0 & neighbors <= length(df$x)]
    # Remove the current index
    neighbors <- setdiff(neighbors, i)
    # Fit a linear model using the neighboring values
    fit <- lm(y ~ x, data = data.frame(x=df$x[neighbors],y=df$y[neighbors]))
    
    # Predict the missing value using the fitted model
    df$y[i] <- predict(fit, newdata = data.frame(x = df$x[i]))
  }
  
  return(df$y)
}


# use the IQR samples to run LOESS correction
dat.corrected<-do.call(cbind,lapply(colnames(dat),function(m){
  x<-dat[,m]
  # only process if it has significant correlation with run order
  if(m%in%attention.metab){
    tmp<-data.table(ro=s.order[which(in_run)][['order']],sample=s.order[which(in_run)][['sample']],metabolite=x)
    x1<-tmp[metabolite!=0][['metabolite']]
    lq <- quantile(x1, 0.25);uq <- quantile(x1, 0.75)
    iqr <- uq - lq
    sel_ind<-tmp[,.I[metabolite >= lq - 1.5 * iqr & metabolite <= uq + 1.5 * iqr & metabolite!=0]]
    if(length(sel_ind)>=50){
      ll<-loess(data = tmp[sel_ind],metabolite~ro)
      # first interpolate
      aa <- approx(x = tmp[sel_ind]$ro, y = ll$fitted, xout = tmp$ro)
      # If there are NAs, then extrapolate NAs with neighboring values
      if(any(is.na(aa$y))) aa$y<-impute_missing(aa)
      # center the trend
      norm.f<-exp(mean(log(aa$y[which(grepl('Pool',tmp$sample))]),na.rm = T))
      fluctuation<-aa$y/norm.f
      tmp[,corrected:=metabolite/fluctuation]
      na.ind<-which(is.na(tmp$corrected))
      if(length(na.ind)==0)
        x<-tmp$corrected
    }
  }
  x
}))
# COMPARE CV before and after correction
pool.ind<-grep('Pool',rownames(dat.corrected))
colnames(dat.corrected)<-colnames(dat)
rsd_df<-data.table(metabolite=colnames(dat.corrected),
                   original=apply(dat[pool.ind,],2,function(x){sd(log10(x+1),na.rm = T)/mean(log10(x+1))}) ,
                   corrected=apply(dat.corrected[pool.ind,],2,function(x){sd(log10(x+1),na.rm = T)/mean(log10(x+1))}) )
ggplot(rsd_df,aes(x=original,y=corrected))+geom_point()+theme_bw()
ggplot(rsd_df,aes(x=original,y=corrected))+geom_point()+theme_bw()+xlim(0,.1)+ylim(0,.1)

# correct by stable IDs and log transformation
plot(density(rsd_df$corrected,na.rm = T),xlim=c(0,.1))
mm<-apply(dat.corrected,2,function(x) median(x,na.rm = T))
dat.sub<-dat.corrected[,which(mm>1e5 & mm<1e10 & colSums(dat.corrected==0)==0 & rsd_df$corrected<.05)]
norm.f<-apply(apply(dat.sub,2,function(x) x/median(x)),1,median)
dat.corrected[]<-apply(dat.corrected,2,function(x) log10(x/norm.f+1))
colnames(dat.corrected)<-colnames(dat)
# 
dat.corrected[1:4,1:4]
rownames(dat.corrected)
dat.corrected<-dat.corrected[grep('Pool',rownames(dat.corrected),invert = T),]
rownames(dat.corrected)<-s.info[match(rownames(dat.corrected),Tube)][['Mouse.ID']]

saveRDS(dat.corrected,'work/cleaned_data.rds')
saveRDS(s.info[match(rownames(dat.corrected),Mouse.ID)],'work/sample_info.rds')
