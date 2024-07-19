rm(list=ls());gc()
library(data.table);library(ggplot2);library(limma);library(variancePartition);library(openxlsx);source('/project/CRI/DeBerardinis_lab/lcai/CommonFunctions/beautify_df.R')
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/')
dat.list<-readRDS('2024_paper/work/batch1/cleaned_data.rds')
metab<-dat.list$data
metab<-metab[grep('TIC',rownames(metab),invert = T),]
# write.csv(metab,'output/Feb2023/data.csv',row.names = T)
setDT(dat.list$s.info)
s.info<-dat.list$s.info
################################################################################################
################################ ADDED TO RESPOND ##############################################
################################################################################################
# However, we would now like to limit this to the samples we are discussing in this paper – 
# namely on the PDX samples and their matched original patient sample 
# (i.e. only samples in "File 2_112923 Human PDX Pairs">Sheet "Sample Annotations">Column B="yes"). 
# Are you able to generate these figures again using this smaller dataset?
s.info2<-read.xlsx('2024_paper/data/batch1/File 2_112923 Human PDX Pairs.xlsx',sheet=3)
setDT(s.info2)
identical(colnames(metab),s.info2$Sample)
s.info2[,.N,by=PAIR.FOR.ANALYSIS]
#    PAIR.FOR.ANALYSIS  N
# 1:                no 79
# 2:               yes 36
keep.ind<-s.info2[,.I[PAIR.FOR.ANALYSIS=='yes']]
metab<-metab[,keep.ind];s.info<-s.info[match(s.info2$Sample[keep.ind],Sample)]
stat.1<-readRDS('2024_paper/work/batch1/multi_sample_stat.rds')
stat.1<-data.table(metabolite=rownames(stat.1),stat.1)

########################################
hmdb_pw<-readRDS('/project/CRI/DeBerardinis_lab/lcai/CommonData/metaboanalyst/2023/pathway_qs/HMDB_pathway_list.rds')
metab.ref<-fread('2024_paper/work/batch1/processed_metab_ref.csv')
stat.1[,HMDB:=metab.ref[match(stat.1$metabolite,metab.ref$Metabolite)][['official_HMDB']]]
stat.1[nchar(stat.1$HMDB)==0,HMDB:='']
source('/project/CRI/DeBerardinis_lab/lcai/CommonFunctions/HyperGeometricTest.R')

stat.df.1<-do.call(rbind,lapply(names(hmdb_pw)[c(5,7,9,11,12,14,15)],function(pw){
  bg<-stat.1[['HMDB']] #unique(unlist(hmdb_pw[[pw]])) 
  rbind(data.table(change='high in PDX',lib=pw,LC.hyper2(genesets = hmdb_pw[[pw]],genes = bg,testset = stat.1[P.Value<.05&logFC>0][['HMDB']],include.intersection = T)),
        data.table(change='high in Human',lib=pw,LC.hyper2(genesets = hmdb_pw[[pw]],genes = bg,testset = stat.1[P.Value<.05&logFC<0][['HMDB']],include.intersection = T)))
}))
stat.df.2<-do.call(rbind,lapply(names(hmdb_pw)[c(5,7,9,11,12,14,15)],function(pw){
  bg<-stat.1[['HMDB']] #unique(unlist(hmdb_pw[[pw]])) 
  data.table(lib=pw,LC.hyper2(genesets = hmdb_pw[[pw]],genes = bg,testset = stat.1[P.Value<.05][['HMDB']],include.intersection = T))
}))
stat.df.1[,padj:=p.adjust(pv,'BH'),by=.(lib,change)];stat.df.2[,padj:=p.adjust(pv,'BH'),by=lib]

source('/project/CRI/DeBerardinis_lab/lcai/CommonFunctions/beautify_df.R')
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper/')
stat.df.2[,overlapping_metabolites:=unlist(sapply(overlapping.genes,function(x) paste(metab.ref[match(unlist(strsplit(x,split = ',')),official_HMDB)][['Metabolite']],collapse = ',')))];stat.df.2[,overlapping.genes:=NULL]
stat.df.2[lib=='kegg_pathway'][order(pv)][1:10]
stat.df.2[lib=='smpdb_pathway'][order(pv)][1:10] #
stat.df.2[lib=='main_class'][order(pv)][1:10] # 

pdf('output/batch1/Fig2H_pw_anal_KEGG.pdf',w=3.6,h=2.5)
ggplot(stat.df.2[lib=='kegg_pathway'],aes(x=enrichment.ratio,y=-log10(pv),label=feature))+geom_point(alpha=.3,color='magenta4',mapping=aes(size=enrichment.ratio))+
  geom_text_repel(stat.df.2[lib=='kegg_pathway'][order(pv)][pv<.05],mapping=aes(x=enrichment.ratio,y=-log10(pv),label=feature),size=4.5)+theme_bw()+
  theme(legend.position = 'bottom',legend.direction = 'horizontal',legend.background = element_blank())+
  guides(size = "none")+ggtitle("Pathway analaysis by KEGG pathways")
dev.off()




##### Make Heatmaps ##### 
library(ComplexHeatmap)

pyrimidines<-unlist(strsplit(stat.df.2[lib=='kegg_pathway'&feature=='Pyrimidine metabolism'&pv<.05][['overlapping_metabolites']],split = ','))
pyrimidines<-setdiff(stat.1[metabolite%in%pyrimidines][order(logFC)][['metabolite']],'Glutamine')
purines<-unlist(strsplit(stat.df.2[lib=='kegg_pathway'&feature=='Purine metabolism'][['overlapping_metabolites']],split = ','))
purines<-setdiff(stat.1[metabolite%in%purines][order(logFC)][['metabolite']],'Glutamine')
BA<-unlist(strsplit(stat.df.2[lib=='kegg_pathway'&feature=='Primary bile acid biosynthesis'&pv<.05][['overlapping_metabolites']],split = ','))
pyrimidines<-stat.1[metabolite%in%pyrimidines][order(logFC)][['metabolite']];BA<-stat.1[metabolite%in%BA][order(logFC)][['metabolite']]
purines<-stat.1[metabolite%in%purines][order(logFC)][['metabolite']]

# BA<-purines
test.pw<-'Citrate cycle (TCA cycle)';test<-unlist(strsplit(stat.df.2[lib=='kegg_pathway'&feature==test.pw][['overlapping_metabolites']],split = ','));test<-stat.1[metabolite%in%test][order(logFC)][['metabolite']]
test<-c('Citrate/Isocitrate','Fumarate','Malate','Pyruvate')
dat.sub<-t(metab[c(pyrimidines,BA),])
dat.sub[]<-apply(dat.sub,2,scale)
avg<-apply(dat.sub[,stat.1[metabolite%in%c(pyrimidines,BA)][logFC>0][['metabolite']]],1,function(x) mean(x,na.rm=T))-
  apply(dat.sub[,stat.1[metabolite%in%c(pyrimidines,BA)][logFC<0][['metabolite']]],1,function(x) mean(x,na.rm=T))

no<-order(s.info$`Human/PDX`,avg)
hm1<-Heatmap(t(dat.sub[no,pyrimidines]),cluster_rows = F,cluster_columns = F,show_column_names = F,name='z score',row_title = 'pyrimidines',
             top_annotation = HeatmapAnnotation(which = 'column',df=data.frame(sample=s.info$`Human/PDX`[no]),
                                                col = list(sample=structure(brewer.pal(6,'Reds')[c(6,1)],names=c('Human','PDX')))))
hm2<-Heatmap(t(dat.sub[no,BA]),cluster_rows = F,cluster_columns = F,show_column_names = F,name='z score',row_title = "Purine") #"Primary\nBile\nAcids"
pdf('output/batch1/Fig2_I_heatmaps_pyrimidine_BA.pdf',w=9,h=2.6);draw(hm1%v%hm2,merge_legend=T);dev.off()
