rm(list=ls());gc()
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper/')
setwd('GitHub_repo/batch2/')
library(data.table);library(openxlsx);library(nlme)
dat.list<-readRDS('cleaned_data_refined.rds')
metab<-dat.list$data
s.info<-dat.list$s.info
wb<-createWorkbook()
addWorksheet(wb,'sample_info');writeDataTable(wb,'sample_info',s.info)
addWorksheet(wb,'metabolomics');writeDataTable(wb,'metabolomics',data.table(sample=colnames(metab),t(metab)))
saveWorkbook(wb,'processed_data_batch2.xlsx',overwrite = T)
# we use a random-effects linear model without an interaction term
pv_df.mxLinear_passage<-data.table(metabolite=colnames(metab),pv=t(data.table(s.info[,1:3,with=F],metab)[Passage!=0,lapply(.SD,function(x){
  tmp<-data.table(x,Passage,Tumor.Name)[!is.na(x)]
  tryCatch({
    fit<-lme(x ~ Passage, random=~1|Tumor.Name,data = tmp)
    c(fit$coefficients$fixed[['Passage']],anova(fit)[c('Passage'),'p-value'])
  },error=function(e) as.numeric(NA))
}
),.SDcols=-(1:3)]))
colnames(pv_df.mxLinear_passage)[2:3]<-c('beta','passage_pv')
pv_df.mxLinear_passage<-pv_df.mxLinear_passage[order(passage_pv)]
pv_df.mxLinear_passage[,passage_padj:=p.adjust(passage_pv,'BH')]
source('../common/beautify_df.R')
pv_df.mxLinear_passage<-beautify.dt(pv_df.mxLinear_passage)
fwrite(pv_df.mxLinear_passage,'TableSx_stats_compare_PDX_by_passage_mixed_effect_linear_with_imputed_data.csv')
metab.ref<-fread('metab_ref.csv')
pv_df.mxLinear_passage[,HMDB:=metab.ref[match(pv_df.mxLinear_passage$metabolite,Metabolite)][['official_HMDB']]]

source('../common/HyperGeometricTest.R')
hmdb_pw<-readRDS('../common/HMDB_pathway_list.rds')
stat.1<-pv_df.mxLinear_passage

stat.df.1<-do.call(rbind,lapply(names(hmdb_pw)[c(5,7,9,11,12,14,15)],function(pw){
  bg<-stat.1[['HMDB']] #unique(unlist(hmdb_pw[[pw]])) 
  rbind(data.table(change='increase',lib=pw,LC.hyper2(genesets = hmdb_pw[[pw]],genes = bg,testset = stat.1[passage_padj<.05&beta>0][['HMDB']],include.intersection = T)),
        data.table(change='decrease',lib=pw,LC.hyper2(genesets = hmdb_pw[[pw]],genes = bg,testset = stat.1[passage_padj<.05&beta<0][['HMDB']],include.intersection = T)))
}))
stat.df.1[,overlapping_metabolites:=unlist(sapply(overlapping.genes,function(x) paste(metab.ref[match(unlist(strsplit(x,split = ',')),official_HMDB)][['Metabolite']],collapse = ',')))]
stat.df.1[,overlapping.genes:=NULL]
stat.df.1[lib=='main_class'][order(pv)][pv<.1]
stat.df.1[lib=='kegg_pathway'][order(pv)][pv<.1]
stat.df.1[lib=='smpdb_pathway'][order(pv)][pv<.1]
pdf('Fig3A_pw_anal_pssg_main_class.pdf',w=4,h=2.8)
ggplot(stat.df.1[lib=='main_class'&!is.na(pv)&!is.na(enrichment.ratio)],aes(x=enrichment.ratio,y=-log10(pv),size=enrichment.ratio,color=change,label=feature))+
  geom_point(alpha=.8)+
  geom_text_repel(stat.df.1[lib=='main_class'&!is.na(pv)&!is.na(enrichment.ratio)][order(pv)][pv<.1],force = 10,mapping=aes(x=enrichment.ratio,y=-log10(pv),label=feature))+
  theme_bw()+scale_color_manual(values = structure(rev(brewer.pal(7,'Reds'))[c(5,1)],names=c('decrease','increase')))+
  theme(legend.position = 'bottom',legend.direction = 'horizontal',legend.background = element_blank())+labs(color='change with passage')+
  guides(size = "none",label='none')+ggtitle("Pathway analaysis by main chemical classes")
dev.off()

pdf('output/batch2/Fig3B_pw_anal_pssg_KEGG.pdf',w=4,h=2.8)
ggplot(stat.df.1[lib=='kegg_pathway'&!is.na(pv)&!is.na(enrichment.ratio)],aes(x=enrichment.ratio,y=-log10(pv),size=enrichment.ratio,color=change,label=feature))+
  geom_point(alpha=.8)+
  geom_text_repel(stat.df.1[lib=='kegg_pathway'&!is.na(pv)&!is.na(enrichment.ratio)][pv<.1],mapping=aes(x=enrichment.ratio,y=-log10(pv),label=feature),show.legend = FALSE)+
  theme_bw()+
  scale_color_manual(values = structure(rev(brewer.pal(7,'Reds'))[c(5,1)],names=c('decrease','increase')))+
  theme(legend.position = 'bottom',legend.direction = 'horizontal',legend.background = element_blank())+labs(color='change with passage')+
  guides(size = "none",color = guide_legend(override.aes = list(size = 4)))+ggtitle("Pathway analaysis by KEGG pathway")
dev.off()
##### Make Heatmaps ##### 
library(ComplexHeatmap)

process.metab<-function(selected.metabolites){
  dat.sub<-metab[,selected.metabolites]
  dat.sub[]<-apply(dat.sub,2,scale)
  up.sub<-stat.1[metabolite%in%selected.metabolites][beta>0][['metabolite']]
  dn.sub<-stat.1[metabolite%in%selected.metabolites][beta<0][['metabolite']]
  if(length(up.sub)>0)  avg.up<-apply(dat.sub[,up.sub,drop=F],1,function(x) mean(x,na.rm=T)) else avg.up<-rep(0,nrow(dat.sub))
  if(length(dn.sub)>0)  avg.dn<-apply(dat.sub[,dn.sub,drop=F],1,function(x) mean(x,na.rm=T)) else avg.dn<-rep(0,nrow(dat.sub))
  list(dat.sub=dat.sub,avg=avg.up-avg.dn)
}
library(scales)
all.ori<-sort(unique(s.info$Tumor.Name))
all.ori<-all.ori[order(as.numeric(gsub("[A-Z]",'',all.ori)))]
ori.map<-do.call(cbind,lapply(all.ori,function(o) ifelse(s.info$Tumor.Name==o,o,'other')))
ori.col<-structure(c(hue_pal()(13),'#e8eff4'),names=c(all.ori,'other'))
mk.hm<-function(processed.list,hm.title=''){
  no<-order(s.info$Passage,processed.list$avg)
  no<-no[s.info$Passage[no]!=0]
  ha<-HeatmapAnnotation(which = 'column',df=data.frame(Passage=paste0('P',s.info$Passage[no])),
                        # origin=ori.map[no,],
                        show_annotation_name = F,#annotation_legend_param = list(origin = list(at = c(all.ori,'other'))),
                        col = list(Passage=structure(brewer.pal(7,'Reds'),names=paste0('P',6:0))
                                   # ,origin=ori.col
                                   ),simple_anno_size = unit(.2, "cm"))
  Heatmap(t(processed.list$dat.sub[no,]),cluster_rows = F,cluster_columns = F,show_column_names = F,name='z score',row_names_gp = gpar(fontsize=9),
          row_title = hm.title,top_annotation = ha,row_title_gp = gpar(fontsize=9))
}

####### Pyrimidines ####### up
pyrimidines<-unlist(strsplit(stat.df.1[grep('pyrimidine',feature,ignore.case = T)][pv<.1][['overlapping_metabolites']],split = ','))
pyrimidines<-stat.1[metabolite%in%pyrimidines][order(beta)][['metabolite']]
pyrmidine.dat<-process.metab(pyrimidines)

####### Glycerophosphocholines ####### down
Glycerophosphocholines<-unlist(strsplit(stat.df.1[feature=='Glycerophosphocholines'][pv<.1][['overlapping_metabolites']],split = ','))
Glycerophosphocholines<-setdiff(stat.1[metabolite%in%Glycerophosphocholines][order(metabolite)][['metabolite']],'Platelet-activating factor')
Glycerophosphocholines.dat<-process.metab(Glycerophosphocholines)

####### Glycerophosphoethanolamines ####### down
Glycerophosphoethanolamines<-unlist(strsplit(stat.df.1[feature=='Glycerophosphoethanolamines'][pv<.1][['overlapping_metabolites']],split = ','))
Glycerophosphoethanolamines<-stat.1[metabolite%in%Glycerophosphoethanolamines][order(metabolite)][['metabolite']]
Glycerophosphoethanolamines.dat<-process.metab(Glycerophosphoethanolamines)


pdf('Fig3C_heatmap_Pssg_composite.pdf',w=4.6,h=4.4)
draw(mk.hm(Glycerophosphocholines.dat,"Glycerophospho-\ncholines")%v%
       mk.hm(Glycerophosphoethanolamines.dat,"Glycerophospho-\nethanolamines")%v%
       mk.hm(pyrmidine.dat,'Pyrimidines'),
     merge_legend=T,show_heatmap_legend=T,show_annotation_legend=T,padding = unit(c(2, 2, 2, 15), "mm"))
dev.off()

# PCA by all metabolites
pca<-prcomp(metab,scale. = T,center = T)
pc.df<-data.table(pca$x[,1:2],s.info)
pc.df[,Passage:=paste0('P',Passage)]
pc.df[,Tumor.Name:=factor(Tumor.Name,levels=all.ori)]
pc.lab<-paste0('PC',1:2,' (',100*round(pca$sdev[1:2]^2/sum(pca$sdev^2),3),'%)')
pdf('Fig3E_PCA_all.pdf',w=4,h=3)
ggplot(pc.df,aes(x=PC1,y=PC2,color=Tumor.Name,alpha=Passage,size=Passage))+geom_point()+theme_bw()+
  scale_color_manual(values=ori.col[-length(ori.col)])+
  scale_alpha_manual(values=structure(seq(from = 0.1, to = 1, length.out = 7),names=paste0('P',0:6)))+
  scale_size_manual(values=structure(seq(from = 5, to = 1, length.out = 7),names=paste0('P',0:6)))+
  xlab(pc.lab[1])+ylab(pc.lab[2])+guides(color = "none")+ggtitle('PCA with all metabolites')
dev.off()
pdf('Fig3F_PCA_by_pt.pdf',w=4,h=3)
ggplot(pc.df,aes(x=PC1,y=PC2,color=Tumor.Name,alpha=Passage,size=Passage))+geom_point()+theme_bw()+
  scale_color_manual(values=ori.col[-length(ori.col)])+
  scale_alpha_manual(values=structure(seq(from = 0.1, to = 1, length.out = 7),names=paste0('P',0:6)))+
  scale_size_manual(values=structure(seq(from = 5, to = 1, length.out = 7),names=paste0('P',0:6)))+
  xlab(pc.lab[1])+ylab(pc.lab[2])+ggtitle('PCA with all metabolites')+facet_wrap(~Tumor.Name,nrow = 3)+
  theme(panel.grid = element_blank(),strip.background = element_rect(fill = NA),legend.position = 'none',panel.spacing = unit(.1,'lines'),
        axis.ticks = element_blank(),axis.text = element_blank())
dev.off()
