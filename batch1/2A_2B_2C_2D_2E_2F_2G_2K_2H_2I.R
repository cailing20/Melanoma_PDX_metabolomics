rm(list=ls());gc()
library(data.table);library(ggplot2);library(limma);library(variancePartition);library(openxlsx);library(reshape2);library(patchwork)
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper')
setwd('GitHub_repo/batch1/')
source('../common/beautify_df.R')
dat.list<-readRDS('cleaned_data.rds')
metab<-dat.list$data
metab<-metab[grep('TIC',rownames(metab),invert = T),]
setDT(dat.list$s.info)
s.info<-dat.list$s.info
################################################################################################
################################ ADDED TO RESPOND ##############################################
################################################################################################
# However, we would now like to limit this to the samples we are discussing in this paper – 
# namely on the PDX samples and their matched original patient sample 
# (i.e. only samples in "File 2_112923 Human PDX Pairs">Sheet "Sample Annotations">Column B="yes"). 
# Are you able to generate these figures again using this smaller dataset?
s.info2<-read.xlsx('File 2_112923 Human PDX Pairs.xlsx',sheet=3)
setDT(s.info2)
identical(colnames(metab),s.info2$Sample)
s.info2[,.N,by=PAIR.FOR.ANALYSIS]
#    PAIR.FOR.ANALYSIS  N
# 1:                no 79
# 2:               yes 36
keep.ind<-s.info2[,.I[PAIR.FOR.ANALYSIS=='yes']]
metab<-metab[,keep.ind];s.info<-s.info[match(s.info2$Sample[keep.ind],Sample)]
# saving the processed data with refined samples
wb<-createWorkbook()
addWorksheet(wb,'sample_info');writeDataTable(wb,'sample_info',s.info)
addWorksheet(wb,'metabolomics');writeDataTable(wb,'metabolomics',data.table(sample=colnames(metab),t(metab)))
saveWorkbook(wb,'processed_data.xlsx',overwrite = T)
################################################################################################
################################################################################################
################################################################################################
# univariate analysis to determine the % variance contribution
table(s.info$`Primary/LN/Met.for.ANALYSIS`)
table(s.info$Pigment.Status)
table(s.info$Sex)
table(s.info$Bin.BRAF)
s.info[Sex==0,Sex:=NA]
selected.features<-c('Human/PDX','Pigment.Status','Bin.BRAF','Sex','Primary/LN/Met.for.ANALYSIS')
univ.stat.list<-structure(lapply(selected.features,function(cl){
  keep.ind<-which(!is.na(s.info[[cl]]))
  form<-as.formula(paste0('~`',cl,'`'))
  design = model.matrix(form,s.info[keep.ind]) 
  fit <- lmFit(metab[,keep.ind], design);eBayesfit <- eBayes(fit,trend = T)
  co.df<-topTable(eBayesfit, sort.by = ifelse(length(unique(setdiff(s.info[[cl]],NA)))>2,'F','P'), n = Inf, coef=grep(cl,colnames(coef(fit)),value = T))
  sel.ind<-which(co.df$adj.P.Val<.05);if(length(sel.ind)>0) rownames(co.df)[sel.ind] # P.Value # adj.P.Val
}),names=selected.features)
univ.count.df<-data.table(feature=names(univ.stat.list),count=unlist(lapply(univ.stat.list,length)))
univ.var.list<-lapply(selected.features,function(cl){
  keep.ind<-which(!is.na(s.info[[cl]]))
  form<-as.formula(paste0('~`',cl,'`'))
  design = model.matrix(form,s.info[keep.ind]) 
  fit <- lmFit(metab[,keep.ind], design);eBayesfit <- eBayes(fit,trend = T)
  varPart <- fitExtractVarPartModel(metab, form, s.info)
  varPart[,1]
})
univ.var.df<-do.call(cbind,univ.var.list)
rownames(univ.var.df)<-rownames(metab)
colnames(univ.var.df)<-selected.features
univ.var.df<-as.data.table(reshape2::melt(univ.var.df))
univ.var.df<-univ.var.df[Var2!='Primary/LN/Met.for.ANALYSIS']
univ.var.df[,Var2:=factor(Var2,levels=selected.features)]
univ.count.df[feature=='Bin.BRAF',feature:='Binary BRAF mutation status'][feature=='Primary/LN/Met.for.ANALYSIS',feature:='Primary/LN/Met'][feature=='Pigment.Status',feature:='Pigment Status']
univ.count.df[,feature:=factor(feature,levels=c("Human/PDX","Pigment Status","Binary BRAF mutation status",'Primary/LN/Met',"Sex"))]
g1<-ggplot(univ.count.df,aes(x=count,y=feature,label=count))+geom_bar(stat='identity',fill='gray')+theme_bw()+theme(axis.title.y = element_blank())+geom_text(hjust=0)+
  xlab('significant metabolites by univariate test')+xlim(0,110)+ggtitle('Univariate association')


unique(s.info[is.na(Sex)][['Patient']])
colSums(is.na(as.matrix(s.info[,match(selected.features[1:4],colnames(s.info)),with=F])))
# multivariate analysis to determine the % variance contribution
form<-as.formula(paste0('~`',paste(selected.features[1:2],collapse = '`+`'),'`'))
multi.var.df<-fitExtractVarPartModel(metab, form, s.info)
tmp<-multi.var.df[,-ncol(multi.var.df)]
multi.var.df<-reshape2::melt(data.table(metabolite=rownames(metab)[which(rowSums(is.na(metab))==0)],tmp))
setDT(multi.var.df)
multi.var.df[,metabolite:=factor(metabolite,levels=multi.var.df[,sum(value),by=metabolite][order(-V1)]$metabolite)]
multi.var.df[variable=='`Human/PDX`',variable:='Human/PDX'][variable=='Pigment.Status',variable:='Pigment Status']

# sample type stats
design = model.matrix(form,s.info) 
fit <- lmFit(metab, design);eBayesfit <- eBayes(fit,trend = T)

sample.stat <- topTable(eBayesfit, sort.by = "P", n = Inf, coef='`Human/PDX`PDX')
pigment.stat <- topTable(eBayesfit, sort.by = "P", n = Inf, coef='Pigment.StatusPigmented')
multi.var.df[variable=='Human/PDX'&metabolite%in%rownames(sample.stat)[sample.stat$adj.P.Val<.05],significance:='yes']
multi.var.df[variable=='Pigment Status'&metabolite%in%rownames(pigment.stat)[pigment.stat$adj.P.Val<.05],significance:='yes']
multi.var.df[is.na(significance),significance:='no']
colnames(multi.var.df)[colnames(multi.var.df)=='variable']<-'feature'
g3<-ggplot(multi.var.df,aes(y=value,x=metabolite,fill=feature,alpha=significance))+geom_bar(stat = 'identity', width=1)+theme_classic()+
  ggtitle('Multivariate association (Human/PDX + Pigmentation Status)')+ylab('variance explained')+
  theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),panel.border = element_rect(fill=NA,size = 1),legend.position = 'none')+
  scale_fill_manual(values=c('Human/PDX'='#DE2D26','Pigment Status'='#553d3a'))+scale_alpha_manual(values=c('yes'=1,'no'=.5))
g1<-g1+ggtitle("Univariate association")+xlab('Number of statistically significant metabolites')
pdf('Fig2AB_univ_multi_side_by_side.pdf',w=11.6,h=2);g1+g3+plot_layout(design = "ABB");dev.off()
lab.df<-data.table(feature=rep(c('Human/PDX','Pigmentation Status'),each=2),significance=rep(c('significant','insignificant'),times=2))
lab.df[,significance:=factor(significance,levels = c('significant','insignificant'))]
lab.df[,feature:=factor(feature,levels = c('Human/PDX','Pigmentation Status'))]

pdf('Fig2B_multivar_legend.pdf',w=3.3,h=1)
ggplot(lab.df,aes(x=significance,y=feature,fill=feature,alpha=significance))+geom_tile(color='black')+theme_minimal()+
  scale_fill_manual(values=c('Human/PDX'='#DE2D26','Pigmentation Status'='#553d3a'))+scale_alpha_manual(values=c('significant'=1,'insignificant'=.5))+
  theme(legend.position = 'none',panel.grid = element_blank(),axis.title= element_blank(),axis.text.y = element_blank())+xlab("Statistical significance")+facet_wrap(~feature,scales = 'free')
dev.off()
wb<-createWorkbook()
sample.stat<-data.table(metabolite=rownames(sample.stat),beautify.dt(sample.stat,pv.key = c('P.Value','adj.P.Val')))
pigment.stat<-data.table(metabolite=rownames(pigment.stat),beautify.dt(pigment.stat,pv.key = c('P.Value','adj.P.Val')))
addWorksheet(wb,'patient_vs_PDX');writeDataTable(wb,'patient_vs_PDX',sample.stat)
addWorksheet(wb,'pigmentation');writeDataTable(wb,'pigmentation',pigment.stat)
saveWorkbook(wb,'multivar_stat.xlsx',overwrite = T)

########## Volcano Plot ###########
sample.stat[,log2FC:=logFC*log(10)/log(2)]
library(RColorBrewer);library(ggrepel)
p.col<-structure(rev(brewer.pal(7,'Reds')),names=as.character(0:6))
pdf('Fig2C_volcano_sample_type.pdf',w=4.2,h=4.2)
ggplot(sample.stat,aes(x=log2FC,y=-log10(adj.P.Val)))+geom_point(color='lightgray')+
  geom_text_repel(sample.stat[adj.P.Val<.05&log2FC<(-1)],mapping=aes(x=log2FC,y=-log10(adj.P.Val),label=metabolite),nudge_y = -.5,nudge_x = -.2,color=p.col[1],size=3,max.overlaps = 10)+
  geom_text_repel(sample.stat[adj.P.Val<.05&log2FC>1],mapping=aes(x=log2FC,y=-log10(adj.P.Val),label=metabolite),nudge_y = -.5,nudge_x=.2,color=p.col[5],size=3,max.overlaps = 10)+
  geom_point(sample.stat[adj.P.Val<.05&log2FC<(-1)],mapping=aes(x=log2FC,y=-log10(adj.P.Val)),color=p.col[1])+
  geom_point(sample.stat[adj.P.Val<.05&log2FC>1],mapping=aes(x=log2FC,y=-log10(adj.P.Val)),color=p.col[5])+
  theme_bw()+
  geom_vline(xintercept = c(-1,1),linetype='dashed')+geom_hline(yintercept = -log10(.05),linetype='dashed')
dev.off()
selected.metab.list<-list(c("1-Hydroxy-2-naphthoate","N-Methylglutamate","4-Guanidinobutanoate"),
                          c("Theophylline","7-Methylxanthine"),
                          'Creatinine',
                          c("Urate", "Allantoin"))

plot.select<-function(i){
  plot.tmp<-reshape2::melt(data.table(source=s.info$`Human/PDX`,t(metab[selected.metab.list[[i]],,drop=F])))
  ggplot(plot.tmp,aes(x=source,y=value,color=source))+
    geom_boxplot(outlier.shape = NA)+geom_jitter(width = .2,alpha=.3,size=1)+facet_wrap(~variable,scales = 'free_y')+theme_classic()+
    scale_color_manual(values=structure(p.col[c(1,5)],names=c('Human','PDX')))+stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))))+ylab('log10 value')+scale_y_continuous(expand = expansion(mult = .1)) +
    theme(plot.background = element_rect(fill = NA),strip.background = element_blank(),strip.text = element_text(size=11),legend.position = 'none',axis.title.x = element_blank())  
}
library(ggpubr)
pdf('Fig2D_boxplot_PDX_vs_Human_1.pdf',w=6.5,h=2.1);plot.select(1);dev.off()
pdf('Fig2F_boxplot_PDX_vs_Human_2.pdf',w=4.4,h=2.1);plot.select(2);dev.off()
pdf('Fig2E_boxplot_PDX_vs_Human_3.pdf',w=2.3,h=2.1);plot.select(3);dev.off()
pdf('Fig2G_boxplot_PDX_vs_Human_4.pdf',w=4.4,h=2.1);plot.select(4);dev.off()

################# visualize distance #################
comp.dist<-function(dat.sub){
  dat.sub[]<-apply(dat.sub,2,scale)
  # dist1<-as.matrix(as.dist(1-cor(t(dat.sub))))
  dist1<-as.matrix(dist(dat.sub))
  dist1[upper.tri(dist1,diag = T)]<-NA
  dist.df<-data.table(melt(dist1))
  setnames(dist.df,c('s1','s2','d'))
  dist.df<-dist.df[!is.na(d)]
  dist.df[,p1:=s.info[match(s1,Sample)][['Patient']]][,t1:=s.info[match(s1,Sample)][['Human/PDX']]][,p2:=s.info[match(s2,Sample)][['Patient']]][,t2:=s.info[match(s2,Sample)][['Human/PDX']]]
  dist.df[,`same patient`:=ifelse(p1==p2,'yes','no')]
  dist.df[t1=='Human' & t2=='Human',comparison:='both human'][t1=='PDX' & t2=='PDX',comparison:='both PDX'][is.na(comparison),comparison:='human vs. PDX']
  dist.df<-dist.df[comparison=='human vs. PDX']
  
  ggplot(dist.df,aes(x=`same patient`,color=`same patient`,y=d))+geom_boxplot(outlier.shape = NA)+ geom_jitter(width = .2)+
    stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))))+scale_y_continuous(expand = expansion(mult = .1))+
    scale_color_manual(values=c(yes='black',no='gray'))+theme_bw()+xlab('Same patient')+ylab("human vs. PDX\npairwise distance")+theme(legend.position = 'none')
}

dev.off()
pdf('Fig2K_distance_compare_boxplot_original.pdf',h=2.1,w=1.5);comp.dist(t(metab));dev.off()



#######################################
hmdb_pw<-readRDS('../common/HMDB_pathway_list.rds')
metab.ref<-fread('processed_metab_ref.csv')
sample.stat[,HMDB:=metab.ref[match(sample.stat$metabolite,metab.ref$Metabolite)][['official_HMDB']]]
sample.stat[nchar(sample.stat$HMDB)==0,HMDB:='']
source('../common/HyperGeometricTest.R')
pw<-'kegg_pathway'
stat.df.2<-data.table(lib=pw,LC.hyper2(genesets = hmdb_pw[[pw]],genes = sample.stat[['HMDB']] ,testset = sample.stat[P.Value<.05][['HMDB']],include.intersection = T))
stat.df.2[,padj:=p.adjust(pv,'BH'),by=lib]
stat.df.2[,overlapping_metabolites:=unlist(sapply(overlapping.genes,function(x) paste(metab.ref[match(unlist(strsplit(x,split = ',')),official_HMDB)][['Metabolite']],collapse = ',')))];stat.df.2[,overlapping.genes:=NULL]
stat.df.2[lib=='kegg_pathway'][order(pv)][1:10]
pdf('Fig2H_pw_anal_KEGG.pdf',w=3.6,h=2.5)
ggplot(stat.df.2[lib=='kegg_pathway'],aes(x=enrichment.ratio,y=-log10(pv),label=feature))+geom_point(alpha=.3,color='magenta4',mapping=aes(size=enrichment.ratio))+
  geom_text_repel(stat.df.2[lib=='kegg_pathway'][order(pv)][pv<.05],mapping=aes(x=enrichment.ratio,y=-log10(pv),label=feature),size=4.5)+theme_bw()+
  theme(legend.position = 'bottom',legend.direction = 'horizontal',legend.background = element_blank())+
  guides(size = "none")+ggtitle("Pathway analaysis by KEGG pathways")
dev.off()




##### Make Heatmaps ##### 
library(ComplexHeatmap)

pyrimidines<-unlist(strsplit(stat.df.2[lib=='kegg_pathway'&feature=='Pyrimidine metabolism'&pv<.05][['overlapping_metabolites']],split = ','))
pyrimidines<-setdiff(sample.stat[metabolite%in%pyrimidines][order(logFC)][['metabolite']],'Glutamine')
purines<-unlist(strsplit(stat.df.2[lib=='kegg_pathway'&feature=='Purine metabolism'][['overlapping_metabolites']],split = ','))
purines<-setdiff(sample.stat[metabolite%in%purines][order(logFC)][['metabolite']],'Glutamine')
BA<-unlist(strsplit(stat.df.2[lib=='kegg_pathway'&feature=='Primary bile acid biosynthesis'&pv<.05][['overlapping_metabolites']],split = ','))
pyrimidines<-sample.stat[metabolite%in%pyrimidines][order(logFC)][['metabolite']];BA<-sample.stat[metabolite%in%BA][order(logFC)][['metabolite']]
purines<-sample.stat[metabolite%in%purines][order(logFC)][['metabolite']]

# BA<-purines
test.pw<-'Citrate cycle (TCA cycle)';test<-unlist(strsplit(stat.df.2[lib=='kegg_pathway'&feature==test.pw][['overlapping_metabolites']],split = ','));test<-sample.stat[metabolite%in%test][order(logFC)][['metabolite']]
test<-c('Citrate/Isocitrate','Fumarate','Malate','Pyruvate')
dat.sub<-t(metab[c(pyrimidines,BA),])
dat.sub[]<-apply(dat.sub,2,scale)
avg<-apply(dat.sub[,sample.stat[metabolite%in%c(pyrimidines,BA)][logFC>0][['metabolite']]],1,function(x) mean(x,na.rm=T))-
  apply(dat.sub[,sample.stat[metabolite%in%c(pyrimidines,BA)][logFC<0][['metabolite']]],1,function(x) mean(x,na.rm=T))

no<-order(s.info$`Human/PDX`,avg)
hm1<-Heatmap(t(dat.sub[no,pyrimidines]),cluster_rows = F,cluster_columns = F,show_column_names = F,name='z score',row_title = 'pyrimidines',
             top_annotation = HeatmapAnnotation(which = 'column',df=data.frame(sample=s.info$`Human/PDX`[no]),
                                                col = list(sample=structure(brewer.pal(6,'Reds')[c(6,1)],names=c('Human','PDX')))))
hm2<-Heatmap(t(dat.sub[no,BA]),cluster_rows = F,cluster_columns = F,show_column_names = F,name='z score',row_title = "Purine") #"Primary\nBile\nAcids"
pdf('Fig2I_heatmaps_pyrimidine_BA.pdf',w=9,h=2.6);draw(hm1%v%hm2,merge_legend=T);dev.off()
