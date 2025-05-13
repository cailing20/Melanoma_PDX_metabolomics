rm(list=ls());gc()
library(data.table);library(patchwork);library(dendextend);library(ComplexHeatmap);library(RColorBrewer);library(ggplot2);library(ggpubr);library(nlme)
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper/')
setwd('GitHub_repo/batch2/')
source('../common/beautify_df.R')
dat.list<-readRDS('cleaned_data_refined.rds')
library(limma);library(openxlsx)

############################
# Data loading and cleaning
############################
metab<-dat.list$data
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
pv_df.mixedLinear<-beautify.dt(pv_df.mixedLinear)
pv_df.mixedLinear[padj<.05]
stat.1<-pv_df.mixedLinear
##################### 
keep.ind<-s.info[,.I[!Tumor.Name%in%c('MP4A','MP5')]]
dat<-metab[keep.ind,];s.info<-s.info[keep.ind]
sel.metab<-stat.1[pv>.1][['metabolite']]
length(sel.metab)
n.all<-integer()
n.matched<-integer()
for(i in 1:length(sel.metab)){
  dat.sub<-dat[,sel.metab[-i]]
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
  n.matched[i]<-nrow(close_df[patient1==patient2]);n.all[i]<-nrow(close_df)
  
}
table(n.matched)
names(n.matched)<-sel.metab
hist(n.matched,breaks = 7,main = "number of PDX (out of 171) matched to\npatient origin after leaving one metabolite out")

# Number of samples and metabolites
dat.sub<-dat[,sel.metab]
dat.sub[]<-apply(dat.sub,2,scale)
n.samples <- nrow(dat.sub)
n.metabolites <- ncol(dat.sub)
identical(s.info$Mouse.ID,rownames(dat.sub))

# Function to calculate pairwise distances and contributions
calculate_distance_contributions <- function(data) {
  # Initialize the 3D array: samples x samples x metabolites
  contribution_array <- array(0, dim = c(n.samples, n.samples, n.metabolites))
  
  # Compute pairwise distances
  for (i in 1:(n.samples - 1)) {
    for (j in (i + 1):n.samples) {
      # Calculate squared differences for each metabolite
      squared_diffs <- (data[i, ] - data[j, ])^2
      
      # Total Euclidean distance squared
      total_squared_distance <- sum(squared_diffs)
      
      # Avoid division by zero
      if (total_squared_distance > 0) {
        # Contribution fraction for each metabolite
        contributions <- squared_diffs #/ total_squared_distance
      } else {
        # If total distance is zero, contributions are zero
        contributions <- rep(0, n.metabolites)
      }
      
      # Store contributions in the array
      contribution_array[i, j, ] <- contributions
      contribution_array[j, i, ] <- contributions # Symmetric
    }
  }
  
  return(contribution_array)
}

# Run the function
contribution_array <- calculate_distance_contributions(dat.sub)

s.info2<-s.info
s.info2[,Tumor.Name.original:=Tumor.Name]
s.info2[,Tumor.Name:=gsub("[A-Z]$","",Tumor.Name)]

uniq.pt<-unique(s.info2$Tumor.Name)
check.match <- vector()
match.list<-list()
for (pt.name in uniq.pt) {
  pt.i<-s.info2[,.I[Tumor.Name==pt.name & Passage == 0]]
  other.pt <- s.info2[,.I[Tumor.Name%in%setdiff(uniq.pt, pt.name) & Passage == 0]]
  for (pdx.i in s.info2[,.I[Tumor.Name == pt.name & Passage != 0]]) {
    # Distance for the matched pair
    m <- matrix(contribution_array[pt.i, pdx.i, ], length(pt.i), n.metabolites)#data.table(matching = 'yes', )
    
    # Distances for the non-matched pairs
    nm <- contribution_array[other.pt, pdx.i, ]#data.table(matching = 'no', )
    
    # Ensure proper summation for comparison
    m_distance <- rowSums(m)  # Total distance for matched
    nm_distances <- rowSums(nm)  # Distances for non-matched
    
    # Update check.match with comparison result
    # check.match <- c(check.match, min(m_distance) <= min(nm_distances))
    if(min(m_distance) <= min(nm_distances)){
      match.list[[length(match.list)+1]]<-list(pt.name,pt.i,pdx.i,m,nm)
    }
  }
}
# Inspect results
table(check.match)
length(match.list)
hm.list<-list()
for(i in 1:30){
  m<-match.list[[i]][[4]]
  min.ind<-which.min(rowSums(m))
  pt<-with(s.info2[match.list[[i]][[2]]][min.ind],paste0(Tumor.Name.original,' P',Passage))
  pdx<-with(s.info2[match.list[[i]][[3]]],paste0(Tumor.Name.original,' P',Passage))
  m<-m[min.ind,,drop=F]
  dist.diff<-m-apply(match.list[[i]][[5]],2,mean)
  mm<-rbind(m,match.list[[i]][[5]])[,order(dist.diff)]
  set.seed(123)
  hm.list[[i]]<-Heatmap(mm,cluster_rows = F,cluster_columns = F,row_title = paste(pt,'-',pdx,'match'),col=circlize::colorRamp2(c(0, 15), c("white", "#8F0AC0")),
                        name="metabolite-specific\ndistance contribution",
                        top_annotation = HeatmapAnnotation(which='column',difference=sort(dist.diff),col=list(difference=circlize::colorRamp2(c(-10, 0,10), c("blue",'white',"red")))),
                        right_annotation = HeatmapAnnotation(which='row',match=rep(c(paste0('yes (',pdx,')'),'no'),c(1,nrow(match.list[[i]][[5]]))),
                                                             col=list(match=structure(c('seagreen','gold'),names=c(paste0('yes (',pdx,')'),'no'))),show_annotation_name = F))
}
hm.list[[1]]
hm.list[[3]]
hm.list[[30]]
dist.diff.list<-list()
match.names<-vector()
for(i in 1:length(match.list)){
  m<-match.list[[i]][[4]]
  pt<-with(s.info2[match.list[[i]][[2]]][min.ind],paste0(Tumor.Name.original,' P',Passage))
  pdx<-with(s.info2[match.list[[i]][[3]]],paste0(Tumor.Name.original,' P',Passage))
  min.ind<-which.min(rowSums(m))
  m<-m[min.ind,,drop=F]
  match.names[i]<-paste(pt,'-',pdx)
  dist.diff.list[[i]]<-m-apply(match.list[[i]][[5]],2,mean)
}
dist.mat<-do.call(rbind,dist.diff.list)
colnames(dist.mat)<-sel.metab;rownames(dist.mat)<-match.names
# Drop Itaconate since it failed Tom Matthew's manual check

pdf('Ext4A_dist_mat.pdf',w=12,h=9)
Heatmap(dist.mat[,which(colSums(dist.mat)!=0 & !grepl('Itaconate',colnames(dist.mat)))],row_names_gp = gpar(fontsize=6),column_names_gp = gpar(fontsize=6),clustering_method_columns = 'ward.D2',
        clustering_method_rows = 'ward.D2',circlize::colorRamp2(c(-10, 0,10), c("blue",'white',"red")),name='difference')
dev.off()

pdf('Ext4B_contrib_metab_hist.pdf',w=4.5,h=4)
hist(rowSums(dist.mat<0),
     xlab="number of metabolites with difference < 0\n(contributing to greater similarity)",main=NULL,col = scales::alpha('navy',.3))
dev.off()
