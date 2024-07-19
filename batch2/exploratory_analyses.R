rm(list=ls());gc()
library(data.table);library(patchwork);library(dendextend);library(ComplexHeatmap);library(RColorBrewer);library(ggplot2);library(ggpubr);library(nlme)
source('/project/CRI/DeBerardinis_lab/lcai/CommonFunctions/beautify_df.R')
library(limma);library(openxlsx)

############################
# Data loading and cleaning
############################
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper/work/batch2/')
# These are the data that have previously been cleaned
dat<-readRDS('cleaned_data.rds')
s.info<-readRDS('sample_info.rds')
# misspelled metabolite
colnames(dat)[colnames(dat)=='Ascrobate']<-'Ascorbate'
# Jennifer notes May 2024: remove TX16
correct.ref<-fread('PDX_sample_update_reference.csv')
# do a little cleaning of the sample data
s.info[,.N,by=Tumor.Name]
s.info[,Tumor.Name:=trimws(Tumor.Name)]
s.info[,Passage:=as.integer(Passage)]
# remove and update sample
identical(rownames(dat),s.info$Mouse.ID)
rm.ind<-s.info[,.I[Tube%in%correct.ref[,.N,by=Tube][N==1][['Tube']]]]
s.info<-s.info[-rm.ind];dat<-dat[-rm.ind,]
correct.ref<-correct.ref[Tube%in%correct.ref[,.N,by=Tube][N==2][['Tube']]]
correct.ref[source=='corrected',Mouse.ID:=c('AC001','AC002', 'AC003')]
correct.ref<-correct.ref[source=='corrected']
for(t in correct.ref[['Tube']]){
  s.info[match(t,Tube)]<-correct.ref[Tube==t,-1,with=F]
}
s.info[match(c(124,140,176),Tube)]
rownames(dat)<-s.info$Mouse.ID
# This set of metabolites were checked by Joao:
# From the list you sent, "3-(4-Hydroxyphenyl)pyruvate", "Trehalose", "Acetylacrylate", "Sphingosine", "Phenylacetate", "Melatonin", "2-Hydroxyphenylacetate", "4-Hydroxybenzoate” 
# have no signal in both QC (and pools) and samples.
# “AICAR” and “N-Acetylserine” have no signal in the QC (and pools).
# Therefore, we may consider including these as problematic metabolites
issue.metab.Joao<-c("3-(4-Hydroxyphenyl)pyruvate", "Trehalose", "Acetylacrylate", "Sphingosine", "Phenylacetate", "Melatonin","2-Hydroxyphenylacetate", "4-Hydroxybenzoate","AICAR", "N-Acetylserine") 
# These are the metabolites that had issues with signal abrupt changes by run order
issue.metab<-c('Carnitine (18:1)','Carnitine (18:0)','Choline',"N-Acetylserine")
# Before we remove them, we check to see if they display interesting trend by sample type
issue.metab<-unique(c(issue.metab.Joao,issue.metab))
# find metabolites with non-missing values in at least 15 samples (take a look before throwing them away)
hist(colSums(dat==0));s.info[,.N,by=Passage]
# to keep statistical analysis robust, we want to throw out metabolites with too many zeros, but we don't want to discard the interesting ones
too.many.zeros.metab<-colnames(dat)[which(colSums(dat!=0)<(nrow(dat)/2))]
tmp<-data.table(s.info,dat[,too.many.zeros.metab])[order(Passage)]
Heatmap(tmp[,-(1:4)],cluster_rows = F,name='signal',right_annotation = HeatmapAnnotation(which = 'row',passage=tmp$Passage))
tmp<-data.table(s.info,dat[,issue.metab])[order(Passage)]
Heatmap(tmp[,-(1:4)],cluster_rows = F,name='signal',right_annotation = HeatmapAnnotation(which = 'row',passage=tmp$Passage))
too.many.zeros.metab<-setdiff(too.many.zeros.metab,issue.metab)
# They do not seem to be correlated with passage, we can remove them.
dat<-dat[,-match(issue.metab,colnames(dat))]
hist(colSums(dat==0),main='Number of zeros in metabolite')
dat[1:4,1:4]
# The current data has lots of zero values
# The minimum of non-zero values is:
min(dat[dat!=0],na.rm = T)
# 3.91
hist(dat[dat!=0],na.rm = T)
# most of the log10 values are between 6 and 9
# We make an imputed dataset to assign the non-zero mimimum to each metabolite
# remove all zeros
dat<-dat[,which(colSums(dat!=0)>4)]
dat.imputed<-dat
dat.imputed[]<-apply(dat,2,function(x){
  x[x==0]<-min(x[x!=0]);x
})
# We calculate a residual data with PDX and patient tumor difference removed
too.many.zeros.metab<-colnames(dat)[which(colSums(dat!=0)<(nrow(dat)/2))]
rm.ind<-match(too.many.zeros.metab,colnames(dat))
fit <- lmFit(t(dat.imputed[,-rm.ind]), model.matrix(~(Passage==0),s.info) )
dat.passage.residuals <- t(residuals(fit, t(dat.imputed[,-rm.ind])))
fit <- lmFit(t(dat.imputed[,-rm.ind]), model.matrix(~Tumor.Name,s.info) )
dat.origin.residuals <- t(residuals(fit, t(dat.imputed[,-rm.ind])))

# Also, a dataset with missing values
dat.missing<-dat;dat.missing[dat==0]<-NA
##############
# Saving data
##############
wb<-createWorkbook()
addWorksheet(wb,'samples');writeDataTable(wb,'samples',s.info)
addWorksheet(wb,'with_zeros');writeDataTable(wb,'with_zeros',as.data.frame(dat),rowNames = T)
addWorksheet(wb,'with_imputation');writeDataTable(wb,'with_imputation',as.data.frame(dat.imputed),rowNames = T)
addWorksheet(wb,'with_NAs');writeDataTable(wb,'with_NAs',as.data.frame(dat.missing),rowNames = T)
addWorksheet(wb,'sample_type_residuals');writeDataTable(wb,'sample_type_residuals',as.data.frame(dat.passage.residuals),rowNames = T)
addWorksheet(wb,'origin_residuals');writeDataTable(wb,'origin_residuals',as.data.frame(dat.origin.residuals),rowNames = T)
saveWorkbook(wb,'../../output/batch2/processed_data.xlsx',overwrite = T)
saveRDS(list(data=dat.imputed,s.info=s.info),'cleaned_data_refined.rds')
#######################################################################
# General dot plot and heatmap functions that we will use along the way
#######################################################################
p.col<-structure(rev(brewer.pal(7,'Reds')),names=as.character(0:6))
plot.m<-function(dat_select,m,g.title='',scatter=T,includeP0=T,add.loess=F,compare=F){
  m_df.sub<-data.table(s.info[,c('Tumor.Name','Passage','in_mouse'),with=F],dat_select)
  tmp<-m_df.sub[,c('Passage','Tumor.Name','in_mouse',m),with=F];colnames(tmp)[4]<-'m'
  if(!includeP0) tmp<-tmp[Passage!=0]
  if(scatter){
    if(compare){
      ggplot(tmp,aes(x=Passage,y=m))+geom_point(color = 'navy',alpha=.3)+
        geom_smooth(tmp[Passage!=0],mapping=aes(x=Passage,y=m),method = 'lm',color='pink',se = F)+facet_wrap(~Tumor.Name,nrow = 2)+theme_bw()+
        ggtitle(g.title)+ylab(m)
    }else{
      facet(ggscatter(tmp,x='Passage',y='m',ylab = m,add = 'reg.line',color = 'navy',alpha=.3,title = g.title)+stat_cor(),facet.by = 'Tumor.Name',scales = 'free_y',nrow=2) 
    }
  }else{
    tmp[,in_mouse:=c('yes','no')[in_mouse+1]]
    ggplot(tmp,aes(y=in_mouse,x=m,fill=as.character(Passage)))+geom_point(pch=21,color='black',alpha=.5)+facet_grid(Tumor.Name~.)+theme_bw()+
      scale_fill_manual(values=p.col)+xlab(m)+ylab('PDX')+theme(strip.text.y = element_text(angle = 0))+labs(fill='Passage')+ggtitle(g.title)
  }
}


make_hm<-function(dat.sub,hm.title,setNA=F,include_p0=T,supplied.order=NULL,fs=7,keep.ind=NULL,col.clust='complete'){
  if(setNA) dat.sub[dat.sub==0]<-NA
  order.by.origin<-F
  # If we don't want to include the human tumors (Passage 0), we remove them from the keep.ind
  if(is.null(keep.ind)){
    keep.ind<-(if(!include_p0) which(s.info$Passage!=0 & rowSums(!is.na(dat.sub))!=0) else 1:nrow(dat.sub))
  }else{
    order.by.origin<-T
  }
  dat.sub<-dat.sub[keep.ind,]
  dat.sub[]<-apply(dat.sub,2,scale)
  if(is.null(supplied.order)){
    # We want the dendrogram to be sorted by passage
    hc <- hclust(dist(dat.sub),method = 'ward.D2')
    dd<-as.dendrogram(hc)
    if(order.by.origin){
      dd2 <- dendextend::rotate(dd, s.info[keep.ind][order(Tumor.Name,Passage,Mouse.ID)][['Mouse.ID']])
    }else{
      dd2 <- dendextend::rotate(dd, s.info[keep.ind][order(Passage,Tumor.Name,Mouse.ID)][['Mouse.ID']])
    }
    
    pt<-s.info[match(labels(dd2),Mouse.ID)]$Tumor.Name
  }else{
    dat.sub<-dat.sub[supplied.order,]
    pt<-s.info$Tumor.Name[supplied.order]
    dd2<-F
  }
  
  pt<-pt[!duplicated(pt)]
  patient<-do.call(cbind,lapply(pt,function(x) ifelse(s.info[match(rownames(dat.sub),Mouse.ID)][['Tumor.Name']]==x,
                                                      ifelse(s.info[match(rownames(dat.sub),Mouse.ID)][['Passage']]==0,'human','mouse'),F)))
  colnames(patient)<-pt
  passage.col<-structure(rev(brewer.pal(length(unique(s.info[keep.ind]$Passage)),'Reds')),names=unique(sort(s.info[keep.ind]$Passage)))
  origin.col<-c('human'=as.character(passage.col['0']),'mouse'='pink','FALSE'="azure")
  passage<-if(is.null(supplied.order)) s.info$Passage[keep.ind] else s.info$Passage[keep.ind][supplied.order]
  ha<-HeatmapAnnotation(which = 'row',Passage = passage,origin=patient,col=list(Passage=passage.col,origin=origin.col))
  if(ncol(dat.sub)>80){
    Heatmap(dat.sub,name='z score',show_row_names = F,show_column_names = F,column_title = hm.title,cluster_rows = dd2,
            clustering_method_columns = col.clust,right_annotation = ha)
  }else{
    Heatmap(dat.sub,name='z score',show_row_names = F,column_names_gp = gpar(fontsize=fs),column_title = hm.title,cluster_rows = dd2,
            clustering_method_columns = col.clust,right_annotation = ha)
  }
}

################################################
# Clustering by all metabolites for all samples
################################################
dev.off()
hist(colSums(dat==0),main='Number of zeros in metabolite')
# we do many many metabolites with zeros, therefore we set different cut-offs for the number of allowed zeros in our output, and we check how zeros (NAs) affect clustering
n<-100
# pdf('../output/heatmap_overall_cluster.pdf',w=10,h=9)
make_hm(dat.imputed[,which(colSums(dat==0)<=n)],setNA = F,
        hm.title=paste0("Cluster by ",length(which(colSums(dat==0)<=n)),' metabolites with ',ifelse(n==0,'no zeros',paste0(n,' or fewer zeros')),', imputed'))
# dev.off()

######################################################
# Clustering by  metabolite residuals for all samples
######################################################
# pdf('../heatmaps_residual_cluster.pdf',w=8.2,h=5)

make_hm(dat.origin.residuals,setNA=F,include_p0 = T,hm.title = 'residuals from regressing on tumor origin')
make_hm(dat.passage.residuals,setNA=F,include_p0 = T,hm.title = 'residuals from regressing on sample type')
make_hm(dat.passage.residuals,setNA=F,include_p0 = F,hm.title = 'residuals from regressing on sample type')
# dev.off()
#####################################################################################################
# Do PDX samples passaged over time have durable metabolic phenotypes (ie do they cluster together)?
#####################################################################################################
n<-50
# pdf('../output/heatmap_PDX_cluster.pdf',w=8.2,h=5)
make_hm(dat.imputed[,which(colSums(dat==0)<=n)],setNA = F,include_p0 = F,
        hm.title=paste0("Cluster by ",length(which(colSums(dat==0)<=n)),' metabolites with ',ifelse(n==0,'no zeros',paste0(n,' or fewer zeros')),', imputed'))
# dev.off()
#################################################################################################################################################################
# Do PDX samples from different sites but from the same patient (ie TX23A and B, MP8A and B, MP9D and F) converge/diverge in the PDX samples and over passaging?
# Focus on the pairs
#################################################################################################################################################################
s.info.tmp<-s.info;s.info.tmp[,pt:=gsub("[A-Z]$","",Tumor.Name)]
# pdf('../output/heatmaps_From_the_same_patients.pdf',w=5.4,h=5)
pair.ind<-s.info.tmp[,.I[pt%in%s.info.tmp[!duplicated(Tumor.Name)][,.N,by=pt][N==2][['pt']]]]
make_hm(dat.imputed[,-rm.ind],hm.title = 'Samples from the same patients',keep.ind = pair.ind)
pair.ind<-s.info.tmp[,.I[Passage!=0&pt%in%s.info.tmp[!duplicated(Tumor.Name)][,.N,by=pt][N==2][['pt']]]]
make_hm(dat.imputed[,-rm.ind],hm.title = 'Samples from the same patients',keep.ind = pair.ind)
# dev.off()
########################################################################
# Metabolites that differ between patient tumor and PDX
########################################################################
s.info[,in_mouse:=ifelse(Passage==0,0,1)]
m_df<-data.table(s.info[,c('Tumor.Name','Passage','in_mouse'),with=F],dat)
m_df.i<-data.table(s.info[,c('Tumor.Name','Passage','in_mouse'),with=F],dat.imputed)
m_df.m<-data.table(s.info[,c('Tumor.Name','Passage','in_mouse'),with=F],dat.missing)

# linear mixed-effect model
get_pv_mixed_linear<-function(m_df.sel){
  stat_df<-data.table(metabolite=colnames(dat),N_zeros=colSums(dat==0),
             t(m_df.sel[,lapply(.SD,function(x) tryCatch({
               tmp<-data.table(x,in_mouse,Tumor.Name)[!is.na(x)]
               fit<-lme(x ~ in_mouse, random=~1|Tumor.Name,data = tmp)
               c(fit$coefficients$fixed[['in_mouse']],anova(fit)[c('in_mouse'),'p-value'])
               },error=function(e) rep(as.numeric(NA),2))),.SDcols=-(1:3)]))
  colnames(stat_df)[3:4]<-c('beta','pv')
  stat_df
}
pv_df.mixedLinear<-rbind(data.table(input='with_zeros',get_pv_mixed_linear(m_df)),
                          data.table(input='with_imputation',get_pv_mixed_linear(m_df.i)),
                          data.table(input='with_NAs',get_pv_mixed_linear(m_df.m)))
pv_df.mixedLinear[,padj:=p.adjust(pv,'BH'),by=input]
pv_df.mixedLinear<-beautify.dt(pv_df.mixedLinear)
rbind(data.table(cutoff='pv<.05',pv_df.mixedLinear[pv<.05][,.N,by=input]),data.table(cutoff='padj<.05',pv_df.mixedLinear[padj<.05][,.N,by=input]))
pv_df1<-data.table(pv_with_zeros=pv_df.mixedLinear[input=='with_zeros'][['pv']],
                   pv_with_imputation=pv_df.mixedLinear[input=='with_imputation'][['pv']],
                   pv_with_NA=pv_df.mixedLinear[input=='with_NAs'][['pv']])
(ggplot(pv_df1,aes(x=-log10(pv_with_zeros),y=-log10(pv_with_imputation)))+geom_point(alpha=.3)+theme_bw())+
  (ggplot(pv_df1,aes(x=-log10(pv_with_zeros),y=-log10(pv_with_NA)))+geom_point(alpha=.3)+theme_bw())+
  (ggplot(pv_df1,aes(x=-log10(pv_with_imputation),y=-log10(pv_with_NA)))+geom_point(alpha=.3)+theme_bw())

# keep using results from imputed data
species.stats<-pv_df.mixedLinear[input=='with_imputation']
diff.df<-m_df.i[,lapply(.SD,mean),by=in_mouse,.SDcols=-(1:3)]
species.stats[,logFC:=(t(diff.df[in_mouse==1,-1])-t(diff.df[in_mouse==0,-1]))[species.stats$metabolite,]]
ggplot(species.stats,aes(x=beta,y=logFC))+geom_point()+geom_abline(intercept = 0,slope = 1) # beta is already log10FC
ggplot(species.stats,aes(x=beta*log(10)/log(2),y=-log10(padj)))+geom_point()
species.stats[,log2FC:=beta*log(10)/log(2)][,beta:=NULL][,logFC:=NULL]
library(ggrepel)
species.stats.lab<-species.stats;species.stats.lab[padj==0,padj:=2.2e-16]
pdf('../output/volcano_PDX_vs_human.pdf',w=6.5,h=3.6)
ggplot(species.stats.lab,aes(x=log2FC,y=-log10(padj)))+geom_point(color='lightgray')+
  geom_text_repel(species.stats.lab[padj<.05&log2FC<(-1)],mapping=aes(x=log2FC,y=-log10(padj),label=metabolite),nudge_y = -.5,nudge_x = -.2,color=p.col[1],size=3)+
  geom_text_repel(species.stats.lab[padj<.05&log2FC>1],mapping=aes(x=log2FC,y=-log10(padj),label=metabolite),nudge_y = -.5,nudge_x=.2,color=p.col[5],size=3)+
  geom_point(species.stats.lab[padj<.05&log2FC<(-1)],mapping=aes(x=log2FC,y=-log10(padj)),nudge_y = -1,color=p.col[1])+
  geom_point(species.stats.lab[padj<.05&log2FC>1],mapping=aes(x=log2FC,y=-log10(padj)),nudge_y = -1,color=p.col[5])+
  theme_bw()+
  geom_vline(xintercept = c(-1,1),linetype='dashed')+geom_hline(yintercept = -log10(.05),linetype='dashed')
dev.off()
m<-'7-Methylxanthine'
m<-'Quinolinate'
plot.m(dat_select = dat,m = m,g.title = 'with zeros',scatter = F)+theme(legend.position = 'none')|
  plot.m(dat_select = dat.imputed,m = m,g.title = 'with imputed values',scatter = F)+theme(legend.position = 'none')|
  plot.m(dat_select = dat.missing,m = m,g.title = 'with NAs',scatter = F)
species.stats<-beautify.dt(species.stats[,-1,with=F])
fwrite(species.stats,'../output/stats_compare_patient_tumor_vs_PDX_mixed_effect_linear_with_imputed_data.csv')

species_metab<-species.stats[padj<.05][['metabolite']]
# pdf('../output/heatmap_metabolites_differing_between_patient_tumors_and_PDX.pdf',w=7,h=4.5)
make_hm(dat.imputed[,species_metab],setNA = F,hm.title=paste0("Cluster by ",length(species_metab)," metabolites\ndiffering between patient tumors and PDX"))
# dev.off()
grep("phthoa|Methylglut|Guanidino|eoph|ethylx|reatinine|Urate|Allantoin",colnames(m_df.i),value = T)
selected.metab.list<-list(c("1-Hydroxy-2-naphthoate","N-Methylglutamate","4-Guanidinobutanoate"),
                          "7-Methylxanthine",
                          "Salicylate",
                          'Creatinine',
                          c("Urate", "Allantoin"))
plot.select<-function(i){
  plot.tmp<-melt(data.table(source=c('Human','PDX')[m_df.i$in_mouse+1],m_df.i[,match(selected.metab.list[[i]],colnames(m_df.i)),with=F]))
  ggplot(plot.tmp,aes(x=source,y=value,color=source))+
    geom_boxplot(outlier.shape = NA)+geom_jitter(width = .2,alpha=.3,size=1)+facet_wrap(~variable,scales = 'free_y')+theme_classic()+
    scale_color_manual(values=structure(p.col[c(1,5)],names=c('Human','PDX')))+stat_compare_means()+ylab('log10 abundance')+scale_y_continuous(expand = expansion(mult = .1)) +
    theme(plot.background = element_rect(fill = NA),strip.background = element_blank(),strip.text = element_text(size=11),legend.position = 'none',axis.title.x = element_blank())  
}
pdf('../output/boxplot_PDX_vs_Human_1.pdf',w=6.4,h=2.5);plot.select(1);dev.off()
pdf('../output/boxplot_PDX_vs_Human_2.pdf',w=2.4,h=2.5);plot.select(2);dev.off()
pdf('../output/boxplot_PDX_vs_Human_3.pdf',w=2.4,h=2.5);plot.select(3);dev.off()
pdf('../output/boxplot_PDX_vs_Human_4.pdf',w=2.4,h=2.5);plot.select(4);dev.off()
pdf('../output/boxplot_PDX_vs_Human_5.pdf',w=4.5,h=2.5);plot.select(5);dev.off()

ggplot(data.table(source=c('Human','PDX')[s.info$in_mouse+1],prcomp(dat.imputed[,species_metab],center = T,scale. = T)$x[,1:2]),
       aes(x=PC1,y=PC2,color=source))+geom_point()+theme_bw()+scale_color_manual(values=structure(p.col[c(1,5)],names=c('Human','PDX')))
#######################################
# NOT different between human and mouse
#######################################
co<-.25
non_species_metab<-setdiff(unique(pv_df.mixedLinear[!is.na(pv)][pv>co][['metabolite']]),too.many.zeros.metab)
# pdf('../../output/heatmap_metabolites_not_differing_between_patient_tumors_and_PDX.pdf',w=7.2,h=4.5)
make_hm(dat.imputed[,non_species_metab],setNA = F,hm.title=paste0("Cluster by ",length(non_species_metab)," metabolites\nnot differing between patient tumors and PDX"))
# dev.off()

########################################################################################################################
# What metabolic features or pathways change over the course of passaging? Is this the same for all tumors? Some tumors?
########################################################################################################################
m_df<-m_df[Passage!=0]
# i is for imputed
m_df.i<-m_df.i[Passage!=0]
# m is for missing
m_df.m<-m_df.m[Passage!=0]

# we use a random-effects linear model with or without an interaction term
get_nlm_pv<-function(dat_select,i.term=T){
  pv_tmp<-data.table(metabolite=colnames(dat),pv=t(dat_select[,lapply(.SD,function(x){
    tmp<-data.table(x,Passage,Tumor.Name)[!is.na(x)]
    tryCatch({
      if(i.term){
        fit<-lme(x ~ Passage*Tumor.Name, random=~1|Tumor.Name,data = tmp)
        c(fit$coefficients$fixed[['Passage']],anova(fit)[c('Passage','Passage:Tumor.Name'),'p-value'])
      }else{
        fit<-lme(x ~ Passage, random=~1|Tumor.Name,data = tmp)
        c(fit$coefficients$fixed[['Passage']],anova(fit)[c('Passage'),'p-value'])
      }
    },error=function(e) rep(as.numeric(NA),ifelse(i.term,2,1)))
  }
),.SDcols=-(1:3)]))
  colnames(pv_tmp)[2:3]<-c('beta','passage_pv')
  if(ncol(pv_tmp)==4) colnames(pv_tmp)[4]<-'interaction_pv'
  pv_tmp
}

pv_df.mxLinear_passage<-get_nlm_pv(m_df.i,i.term = F)
pv_df.mxLinear_passage[,passage_padj:=p.adjust(passage_pv,'BH')]
plot.m(dat_select = dat.imputed,m = 'Stachydrine')
plot.m(dat_select = dat.imputed,m = 'Oleamide')

non_species_metab2<-unique(pv_df.mxLinear_passage[!is.na(passage_pv)][passage_pv>co][['metabolite']])
# pdf('../../output/heatmap_metabolites_not_differing_between_patient_tumors_and_PDX_or_passage.pdf',w=7.2,h=4.5)
make_hm(dat.imputed[,intersect(non_species_metab,non_species_metab2)],setNA = F,
        hm.title=paste0("Cluster by ",length(non_species_metab)," metabolites\nnot differing between patient tumors and PDX\nor by passage"))
# dev.off()
################################################################################################################################################################
# PAUSED ## PAUSED ## PAUSED ## PAUSED ## PAUSED ## PAUSED #
save.image('exploratory_analyses_paused.RData')
# compare p0-p1 and p1-p6
sel_input<-'with_imputation'#'with NAs'
m.stat<-merge(pv_df.mixedLinear[input==sel_input][,c('metabolite','beta','padj')],
              pv_df.mxLinear_passage[input==sel_input&interaction=='no'][,c('metabolite','beta','passage_padj')],by='metabolite')
result.types<-c('significant in both','only significant in P0 vs. PDX','only significant in P1-P6','not significant in either','contradictory')
m.stat[padj<=.05&passage_padj<=.05&sign(beta.x)==sign(beta.y),result:=result.types[1]]
m.stat[padj<=.05&(passage_padj>.05),result:=result.types[2]]
m.stat[passage_padj<=.05&(padj>.05),result:=result.types[3]]
m.stat[passage_padj>.05&padj>.05,result:=result.types[4]]
m.stat[padj<=.05&passage_padj<=.05&sign(beta.x)!=sign(beta.y),result:=result.types[5]]
m.stat[,result:=factor(result,levels=result.types)]
m.stat[,.N,by=result][order(result)]

ggplot(m.stat,aes(x=beta.x,y=beta.y))+geom_point(alpha=.3,aes(color=result),size=1)+
  scale_color_manual(values=structure(c('red','brown','pink','gray','blue'),names=result.types))+
  annotate(geom = 'text',-Inf,Inf,label=with(with(m.stat,cor.test(beta.x,beta.y)),paste0('r = ',round(estimate,2),', pv = ',signif(p.value,2))),vjust=1.1,hjust=-.02)+
  xlab('non-P0 vs. P0 beta')+ylab('P1-P6 beta')+theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 3,alpha=1)))

ggplot(m.stat,aes(x=beta.x,y=beta.y))+geom_point(alpha=.7,aes(color=result),size=.6)+
  # geom_text_repel(m.stat[result==result.types[1]],mapping=aes(label=metabolite),max.overlaps = 100)+
  # geom_text_repel(m.stat[result==result.types[2]],mapping=aes(label=metabolite),max.overlaps = 100)+
  # geom_text_repel(m.stat[result==result.types[3]],mapping=aes(label=metabolite),max.overlaps = 100)+
  geom_text_repel(m.stat[result==result.types[5]],mapping=aes(label=metabolite),max.overlaps = 100)+
  scale_color_manual(values=structure(c('red','brown','pink','gray','blue'),names=result.types))+
  xlab('non-P0 vs. P0 beta')+ylab('P1-P6 beta')+theme_bw()+
  guides(color = guide_legend(override.aes = list(size = 3)))
m<-'Adenylsuccinate';plot.m(dat_select = dat.imputed,m=m,scatter = T,g.title = m,compare = T)
m<-'Guanine';plot.m(dat_select = dat.imputed,m=m,scatter = T,g.title = m,compare = T)
m<-'Deoxycholate';plot.m(dat_select = dat.imputed,m=m,scatter = T,g.title = m,compare = T)
m<-'Urocanate';plot.m(dat_select = dat.imputed,m=m,scatter = T,g.title = m,compare = T)


# Run PCA
source('../../script/June_2023/PCA.R')
pca_output(0)
sel_metab<-colnames(dat.imputed)[colSums(dat==0)==0]
pdf('../../output/June2023/PCA_plots/PCA_by_all_metabolites_wo_zeros_scree.pdf',w=2.3,h=2.8)
pca_output(sel_metab = sel_metab,1);dev.off()
pdf('../../output/June2023/PCA_plots/PCA_by_all_metabolites_wo_zeros_combined.pdf',w=11.75,h=4.5)
pca_output(sel_metab = sel_metab,6);dev.off()
pdf('../../output/June2023/PCA_plots/PCA_by_all_metabolites_wo_zeros_facet_PC12.pdf',w=10,h=3.4)
pca_output(sel_metab = sel_metab,7);dev.off()
pdf('../../output/June2023/PCA_plots/PCA_by_all_metabolites_wo_zeros_facet_PC34.pdf',w=10,h=3.4)
pca_output(sel_metab = sel_metab,8);dev.off()

# we use a random-effects linear model with interaction terms
sel_metab<-pv_df.mxLinear_passage[input=='with_imputation'&interaction=='no'][passage_padj<.5e-2][['metabolite']]
pdf('../../output/June2023/PCA_plots/PCA_by_passage_metabolites_scree.pdf',w=2.3,h=2.8)
pca_output(sel_metab = sel_metab,1);dev.off()
pdf('../../output/June2023/PCA_plots/PCA_by_passage_metabolites_PC12.pdf',w=8.45,h=6.8)
pca_output(sel_metab = sel_metab,4);dev.off()
pdf('../../output/June2023/PCA_plots/PCA_by_passage_metabolites_PC1.pdf',w=7,h=1.5)
pca_output(sel_metab = sel_metab,2)+theme(legend.position = 'none');dev.off()
pdf('../../output/June2023/PCA_plots/PCA_by_passage_metabolites_PC1_facet.pdf',w=8.3,h=12)
pca_output(sel_metab = sel_metab,3)+theme(legend.position = 'none');dev.off()

pdf('../../output/June2023/heatmaps/metabolites_changed_by_passage.pdf',w=12.9,h=8.7)
draw(make_hm(dat.imputed[,sel_metab],setNA = F,include_p0 = T,fs = 10,col.clust = 'ward.D2',
             hm.title="metabolites correlated with passage with padj < .05,\nsamples ordered by PC1",
             supplied.order = order(-prcomp(dat.imputed[,sel_metab],center = T,scale. = T)$x[,1])),merge_legend=T)
dev.off()

# Example that heatmap may not be the best visualization
plot.m(dat_select = dat.missing,m = 'Cystathionine',scatter = T)
ggscatter(m_df.i,x='Passage',y='Cystathionine',color='navy',alpha=.3,add = 'reg.line')+stat_cor()
################################################################################################################################
# What metabolic features or pathways are conserved over the course of passaging? Is this the same for all tumors? Some tumors?
################################################################################################################################
by_origin_pv<-function(m_df_select,use.cor=T){
  stat.df<-do.call(rbind,lapply(unique(m_df_select$Tumor.Name),function(tum){
    data.table(tum,colnames(dat),
               t(m_df_select[Tumor.Name==tum,lapply(.SD,function(x) 
                 tryCatch({
                   if(use.cor){
                     unlist(cor.test(x,Passage)[c('estimate','p.value')])
                   }else{
                     coefficients(summary(lm(x~Passage)))['Passage',c('t value','Pr(>|t|)')]
                   }
                 },error=function(e) as.numeric(rep(NA,2)))),.SDcols=-(1:3)]))
  }))
  colnames(stat.df)<-c('origin','metabolite','r','pv');stat.df
}

pv_df.by_origin<-rbind(data.table(input='with_zeros',by_origin_pv(m_df[use.ind])),
                       data.table(input='with_imputation',by_origin_pv(m_df.i[use.ind])),
                       data.table(input='with_NAs',by_origin_pv(m_df.m[use.ind])))
pv_df.by_origin[,padj:=p.adjust(pv,'BH'),by=.(input)]
sel_input<-'with_imputation'#'with NAs'
tmp<-pv_df.by_origin[input==sel_input][abs(r)>.3,.N,by=origin][order(origin)]
tmp[,sample_size:=m_df[,.N,by=Tumor.Name][match(tmp$origin,Tumor.Name)][['N']]]
tmp[,`abs(r)>.3`:=N]
tmp[,`abs(r)>.4`:=pv_df.by_origin[input==sel_input][abs(r)>.4,.N,by=origin][match(tmp$origin,origin)][['N']]]
tmp[,`pv<.05`:=pv_df.by_origin[input==sel_input][pv<.05,.N,by=origin][match(tmp$origin,origin)][['N']]]
tmp[,`padj<.05`:=pv_df.by_origin[input==sel_input][padj<.05,.N,by=origin][match(tmp$origin,origin)][['N']]]
tmp[,N:=NULL];tmp[order(origin)]
pv_df.by_origin[,direction:=factor(ifelse(r>0,'up','down'),levels=c('up','down'))]
pv_df.by_origin[,change:=ifelse(r>0,1,-1)]
pv_df.by_origin[input==sel_input&abs(r)>.4&metabolite=='Cystathionine']
library(reshape2);library(ggrepel)
# r.mat<-acast(pv_df.by_origin[input==sel_input],origin~metabolite,value.var = 'r',fill = 0)
# plot(hclust(dist(r.mat),method = 'ward.D2'))
pc<-prcomp(acast(pv_df.by_origin[input==sel_input],origin~metabolite,value.var = 'r',fill = 0))$x[,1:4]
g1<-ggplot(data.table(patient=gsub("[A-Z]$",'',rownames(pc)),origin=rownames(pc),pc),aes(x=PC1,y=PC2,label=origin))+
  geom_point(alpha=.3,size=5,color='pink')+geom_text_repel()+geom_line(aes(group=patient),alpha=.3)+theme_bw()+
  ggtitle('PCA by correlation between passage and metabolites')
pc<-prcomp(acast(pv_df.by_origin[input==sel_input],origin~metabolite,value.var = 'change',fill = 0))$x[,1:4]
g2<-ggplot(data.table(patient=gsub("[A-Z]$",'',rownames(pc)),origin=rownames(pc),pc),aes(x=PC1,y=PC2,label=origin))+
  geom_point(alpha=.3,size=5,color='pink')+geom_text_repel()+geom_line(aes(group=patient),alpha=.3)+theme_bw()+
  ggtitle('PCA by direction of change between passage and metabolites')
g1|g2

sig.metab<-pv_df.mxLinear_passage[input==sel_input&interaction=='no'][passage_padj<.05][order(passage_padj)][['metabolite']]
r.mat<-acast(pv_df.by_origin[input==sel_input],origin~metabolite,value.var = 'r',fill = NA)
pdf('../../output/June2023/heatmaps/correlation_between_metabolite_and_passage_by_PDX_origin.pdf',w=13.5,h=4.7)
Heatmap(r.mat[,sig.metab],clustering_method_rows = 'ward.D2',name='r',column_names_gp = gpar(fontsize=9),
        column_title = "Correlation between metabolite and passage")
dev.off()


N_df<-pv_df.by_origin[input==sel_input&pv<.05,list(N=.N*change),by=.(metabolite,change,direction)]
sum_df<-N_df[,sum(abs(N)),by=metabolite][order(V1)]
N_df[,metabolite:=factor(metabolite,levels=sum_df$metabolite)]
up_metab<-pv_df.mxLinear_passage[input==sel_input&interaction=='no'&beta>0][passage_padj<.05][order(passage_padj)][['metabolite']]
N_df.up<-N_df[metabolite%in%up_metab]
up.up<-N_df.up[direction=='up'][order(N)];up.dn<-N_df.up[direction=='down'][order(N)];up.up[match(up.dn$metabolite,metabolite),N:=N-up.dn$N*.1]
up_metab<-up.up[order(N)][['metabolite']];N_df.up[,metabolite:=factor(metabolite,levels=up_metab)]
g1<-ggplot(N_df.up,aes(x=N,y=metabolite,fill=direction))+geom_bar(stat='identity')+theme_bw()+ggtitle('Increased by passaging')

dn_metab<-pv_df.mxLinear_passage[input==sel_input&interaction=='no'&beta<0][passage_padj<.05][order(passage_padj)][['metabolite']]
N_df.dn<-N_df[metabolite%in%dn_metab]
dn.up<-N_df.dn[direction=='up'][order(N)];dn.dn<-N_df.dn[direction=='down'][order(N)];dn.dn[match(up.dn$metabolite,metabolite),N:=N-up.dn$N*.1]
dn_metab<-dn.dn[order(-N)][['metabolite']];N_df.dn[,metabolite:=factor(metabolite,levels=dn_metab)]
g2<-ggplot(N_df.dn,aes(x=N,y=metabolite,fill=direction))+geom_bar(stat='identity')+theme_bw()+ggtitle('Decreased by passaging')
(g1+theme(legend.position = 'none'))|g2


stat.wb<-createWorkbook()
addWorksheet(stat.wb,'lm_patient_vs_PDX');writeDataTable(stat.wb,'lm_patient_vs_PDX',pv_df.simpleLinear)
addWorksheet(stat.wb,'mx_lm_patient_vs_PDX');writeDataTable(stat.wb,'mx_lm_patient_vs_PDX',pv_df.mixedLinear)
addWorksheet(stat.wb,'mx_lm_wint_patient_vs_PDX');writeDataTable(stat.wb,'mx_lm_wint_patient_vs_PDX',pv_df.mixedLinearWithInt)
addWorksheet(stat.wb,'mx_lm_PDX_Passage');writeDataTable(stat.wb,'mx_lm_PDX_Passage',pv_df.mxLinear_passage)
addWorksheet(stat.wb,'individual_PDX_Passage');writeDataTable(stat.wb,'individual_PDX_Passage',pv_df.by_origin)
saveWorkbook(stat.wb,'../../output/June2023/stas.xlsx',overwrite = T)

# We have previously identified many differences between human samples and PDX samples. Are there metabolic features that are conserved or maintained?
prev.stats<-as.data.table(read.xlsx('../../output/Feb2023/patient_vs_PDX_stats.xlsx'))
prev.stats[,metab:=gsub("_pos|_neg",'',prev.stats$metabolite)]
keep.ind<-prev.stats[,.I[!metab%in%prev.stats.sub[,.N,by=metab][N==2][['metab']]]]
prev.stats.sub<-merge(prev.stats[grepl('_neg',metabolite),c('metab','P.Value'),with=F],
                      prev.stats[grepl('_pos',metabolite),c('metab','P.Value'),with=F],by='metab')
prev.stats<-rbind(prev.stats[keep.ind],prev.stats[metab%in%prev.stats.sub$metab&grepl('neg',metabolite)])
ggplot(prev.stats.sub,aes(x=-log10(P.Value.x),y=-log10(P.Value.y)))+geom_point()+geom_abline(intercept = 0,slope = 1)

prev.stats<-prev.stats[metab%in%pv_df.mixedLinear$metabolite]
m_stat<-merge(pv_df.mixedLinear[input=='with_imputation',c('metabolite','beta','padj')],
      prev.stats[,c('metabolite','logFC','adj.P.Val')],by='metabolite')
file.edit('/project/CRI/')

freq.df<-m_stat[,.N,by=.(logFC<0,beta<0,padj<=.05,adj.P.Val<=.05)]
colnames(freq.df)[1:4]<-c('prev_change','current_change','prev_significant','current_significant')
freq.df[,prev_change:=ifelse(prev_change,'down','up')][,current_change:=ifelse(current_change,'down','up')]
freq.df[,prev_significant:=ifelse(prev_significant,'padj<.05','padj>=.05')]
freq.df[,current_significant:=ifelse(current_significant,'padj<.05','padj>=.05')]
freq.df[,prev_significant:=factor(prev_significant,levels=c('padj>=.05','padj<.05'))]
fisher.test(with(m_stat[padj<=.05&adj.P.Val<=.05],table(sign(logFC),sign(beta))))
pdf('../../output/June2023/count_by_change_and_significance.pdf',w=4.4,h=3.52)
ggplot(freq.df,aes(x=prev_change,y=current_change,size=N,fill=N,label=N))+
  geom_point(pch=21)+geom_text(nudge_x = .3,nudge_y = .3)+scale_fill_viridis_c()+
  facet_grid(prev_significant~current_significant)+theme_bw()+
  xlab('previous PDX - patient tumor')+ylab('current PDX - patient tumor')+
  ggtitle('count by change and significance')
dev.off()
pdf('../../output/June2023/compare_changes.pdf',w=3.75,h=3.52)
ggplot(m_stat,aes(x=-log10(padj)*sign(beta),y=-log10(adj.P.Val)*sign(logFC)))+
  geom_point(alpha=.5)+
  geom_vline(xintercept = c(-log10(.05),log10(.05)),linetype='dashed',alpha=.7)+
  geom_hline(yintercept = c(-log10(.05),log10(.05)),linetype='dashed',alpha=.7)+
  xlab("This dataset\n-log10(padj)*sign(PDX-patient tumor)")+
  ylab("Previous dataset\n-log10(padj)*sign(PDX-patient tumor)")+
  geom_abline(intercept = 0,slope = 1)+
  theme_bw()+ggtitle('Compare PDX-patient tumor changes')
dev.off()
