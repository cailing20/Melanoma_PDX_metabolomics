rm(list=ls());gc()
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/JenniferGill/2024_paper/')
setwd('GitHub_repo/batch2/')
library(data.table);library(openxlsx);library(nlme);library(ggpubr)
dat.list<-readRDS('cleaned_data_refined.rds')
sheets(loadWorkbook('processed_data.xlsx'))
# load data without imputed missing value (with NAs)
dat.wNA<-read.xlsx('processed_data.xlsx',sheet=4,check.names = F)
colnames(dat.wNA)<-gsub('.',' ',colnames(dat.wNA),fixed = T)
metab<-dat.list$data
s.info<-dat.list$s.info
identical(rownames(metab),dat.wNA$row.names)
identical(rownames(metab),s.info$Mouse.ID)
dat.wNA<-dat.wNA[,-1];rownames(dat.wNA)<-rownames(metab)

all.ori<-sort(unique(s.info$Tumor.Name))
all.ori<-all.ori[order(as.numeric(gsub("[A-Z]",'',all.ori)))]
ori.col<-structure(hue_pal()(13),names=all.ori)
saveRDS(ori.col,'ori_col.rds')
ori.map<-do.call(cbind,lapply(all.ori,function(o) ifelse(s.info$Tumor.Name==o,o,'other')))
ori.col2<-structure(c(hue_pal()(13),'#e8eff4'),names=c(all.ori,'other'))
s.info[,Tumor.Name:=factor(Tumor.Name,levels=all.ori)]
select.ind<-s.info[,.I[Passage!=0]]
multi.var.df<-fitExtractVarPartModel(t(metab[select.ind,]), ~Passage, s.info[select.ind,])
summary(multi.var.df[,1])
psg.stat<-topTable(eBayes(lmFit(t(metab[select.ind,]), model.matrix(~Passage,s.info[select.ind,])),trend = T), sort.by = "P", n = Inf, coef='Passage')
multi.var.df<-fitExtractVarPartModel(t(metab[select.ind,]), ~Tumor.Name, s.info[select.ind,])
summary(multi.var.df[,1])
origin.stat<-topTable(eBayes(lmFit(t(metab[select.ind,]), model.matrix(~Tumor.Name,s.info[select.ind,])),trend = T), sort.by = "F", n = Inf)
table(origin.stat$adj.P.Val<.05)
multi.var.df<-fitExtractVarPartModel(t(metab[select.ind,]), ~Passage+Tumor.Name, s.info[select.ind,])
multi.var.df<-as.data.frame(multi.var.df)
multi.var.df<-melt(data.table(metabolite=rownames(multi.var.df),multi.var.df))
setDT(multi.var.df)
multi.var.df[,metabolite:=factor(metabolite,levels=multi.var.df[variable=='Tumor.Name'][order(-value)][['metabolite']])]
multi.var.df[variable=='Tumor.Name',variable:='Origin']
multi.var.df[,variable:=factor(variable,levels = rev(c('Origin','Passage','Residuals')))]
pdf('Fig3D_var_exp.pdf',w=5,h=2)
ggplot(multi.var.df,aes(x=metabolite,y=value,fill=variable))+geom_bar(stat='identity',width=1)+
  scale_fill_manual(values=c('Origin'='#C80002', 'Passage'='#FFCA48', 'Residuals'="#444444"))+theme_minimal()+
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),panel.grid = element_blank(),legend.box.spacing = unit(0,'lines'),
        # axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = -10, r = 0, b = -10, l = 0)),
        legend.position = 'bottom')+ylab('Variance')+
  labs(fill="variance contribution")+ggtitle('Variance partition in PDX')

dev.off()


metab.ref<-fread('metab_ref.csv')
source('../common/HyperGeometricTest.R')
hmdb_pw<-readRDS('../common/HMDB_pathway_list.rds')
# compute concordance between patient and PDX average
median.0.df<-melt(data.table(s.info,metab)[Passage==0,lapply(.SD,median),by=Tumor.Name,.SDcols=-(1:ncol(s.info))])
median.6.df<-melt(data.table(s.info,metab)[Passage==6,lapply(.SD,median),by=Tumor.Name,.SDcols=-(1:ncol(s.info))])
setDT(median.0.df);setDT(median.6.df)
median.df<-data.table(median.0.df,median.6.df[match(with(median.0.df,paste(Tumor.Name,variable)),with(median.6.df,paste(Tumor.Name,variable)))][['value']])
colnames(median.df)[2:4]<-c('metabolite','P0','P6')
cor.stat<-median.df[,cor.test(P0,P6)[c('estimate','p.value')],by=metabolite]
cor.stat[,variance:=multi.var.df[variable=='Origin'][match(cor.stat$metabolite,as.character(metabolite))][['value']]]
library(ggrepel);library(ggside)

length(as.character(cor.stat[estimate>.3][['metabolite']]))
char.metab<-as.character(multi.var.df[variable=='Origin'][value>.3][['metabolite']])
length(char.metab)
char.metab<-intersect(char.metab,as.character(cor.stat[estimate>.3][['metabolite']]))
length(char.metab)
char.metab<-setdiff(char.metab,'Platelet-activating factor')
length(char.metab)
char.hmdb<-setdiff(metab.ref[match(char.metab,Metabolite)][['official_HMDB']],'')
stat.df.1<-do.call(rbind,lapply(names(hmdb_pw)[c(5,7,9,11,12,14,15)],function(pw){
  bg<-setdiff(metab.ref[Metabolite%in%multi.var.df$metabolite][['official_HMDB']],'')
  data.table(lib=pw,LC.hyper2(genesets = hmdb_pw[[pw]],genes = bg,testset = char.hmdb,include.intersection = T))
}))
stat.df.1[,overlapping_metabolites:=unlist(sapply(overlapping.genes,function(x) paste(metab.ref[match(unlist(strsplit(x,split = ',')),official_HMDB)][['Metabolite']],collapse = ',')))]
stat.df.1[,overlapping.genes:=NULL]
stat.df.1[lib=='kegg_pathway'][order(pv)][1:10] # Glutathione metabolism: Cys-Gly,Glutathione disulfide,NADP+,Ornithine
stat.df.1[lib=='smpdb_pathway'][order(pv)][1:10] # Lysine Degradation, Carnitine Synthesis, Phospholipid Biosynthesis
stat.df.1[lib=='super_class'][order(pv)][1:10] # too broad
stat.df.1[lib=='main_class'][order(pv)][1:10] # Fatty esters (carnitines)
stat.df.1[lib=='sub_class'][order(pv)][1:10] # Hydroxy Fatty Acids
stat.df.1[order(pv)][1:10]
nonspec.metab<-as.character(multi.var.df[variable=='Origin'][value<.3][['metabolite']])
nonspec.metab<-setdiff(nonspec.metab,'Platelet-activating factor')
nonspec.metab<-setdiff(metab.ref[match(nonspec.metab,Metabolite)][['official_HMDB']],'')
stat.df.2<-do.call(rbind,lapply(names(hmdb_pw)[c(5,7,9,11,12,14,15)],function(pw){
  bg<-setdiff(metab.ref[Metabolite%in%multi.var.df$metabolite][['official_HMDB']],'')
  data.table(lib=pw,LC.hyper2(genesets = hmdb_pw[[pw]],genes = bg,testset = nonspec.metab,include.intersection = T))
}))
stat.df.2[,overlapping_metabolites:=unlist(sapply(overlapping.genes,function(x) paste(metab.ref[match(unlist(strsplit(x,split = ',')),official_HMDB)][['Metabolite']],collapse = ',')))]
stat.df.2[,overlapping.genes:=NULL]
stat.df.2[lib=='kegg_pathway'][order(pv)][1:10] # pyrimidines
stat.df.2[lib=='smpdb_pathway'][order(pv)][1:10] # 
stat.df.2[lib=='super_class'][order(pv)][1:10] # 
stat.df.2[lib=='main_class'][order(pv)][1:10] # nucleic acids
stat.df.2[lib=='sub_class'][order(pv)][1:10]
stat.df.2[order(pv)][1:10]
cor.stat[,HMDB:=metab.ref[match(cor.stat$metabolite,Metabolite)][['official_HMDB']]]
cor.stat[,class:='other']
pyrimidines<-intersect(hmdb_pw$main_class$Pyrimidines,setdiff(intersect(hmdb_pw$kegg_pathway$`Pyrimidine metabolism`,
                                                           hmdb_pw$smpdb_pathway$`Pyrimidine Metabolism`),''))
cor.stat[HMDB%in%pyrimidines,class:='Pyrimidines']
cor.stat[class=='Pyrimidines']
cor.stat[grep('carnitine',metabolite,ignore.case = T),class:='Carnitines']
cor.stat[grep('hydroxycarnitine',metabolite,ignore.case = T),class:='Hydroxycarnitines']
cor.stat[,highlight:=class!='other']
cor.stat[,label:=F]
#'Carnitine (3:0)','3-Hydroxymethylglutarate','Cytosine','Deoxyuridine',
select.m<-c('2-Hydroxyglutarate','Myoinositol','Guanosine','Cystathionine'#,'Lyso-PC (16:0)','Hypotaurine','Mevalonate'
            )
for(m in select.m) {cor.stat[metabolite==m,label:=T][metabolite==m,highlight:=T]}
plot.m2<-function(m){
  ggplot(median.df[metabolite==m],aes(x=P0,y=P6))+geom_point(mapping=aes(color=Tumor.Name),size=3,alpha=.8)+
    theme_classic()+scale_color_manual(values=ori.col)+ylab('P6 median')+
    theme(legend.position = 'none')+geom_abline(slope = 1,intercept = 0)+
    geom_smooth(method='lm',alpha=.1)+stat_cor()+
    ggtitle(m)
}
cor.stat[,class:=factor(class,levels=c('Carnitines','Hydroxycarnitines','Pyrimidines','other'))]
pdf('Fig3G_cor_vs_var_exp.pdf',w=5.1,h=2.5)
ggplot(cor.stat,aes(x=variance,y=estimate))+
  geom_point(mapping=aes(color=class,alpha=highlight,size=highlight))+theme_bw()+
  geom_text_repel(cor.stat[which(label)],mapping=aes(label=metabolite))+
  xlab('Variance explained by origin')+ylab('Correlation between P0 and P6')+
  scale_alpha_manual(values=c('TRUE'=1,'FALSE'=.5))+scale_size_manual(values=c('TRUE'=2,'FALSE'=.5))+
  scale_color_manual(values = structure(c(brewer.pal(3,'Set1'),'#444444'),names=c('Pyrimidines','Carnitines','Hydroxycarnitines','other')))+
  geom_ysidedensity(aes(x= after_stat(density)),alpha= 0.4,size= 1,position = "stack",color='deeppink4')+
  geom_ysideabline(slope = 0,intercept = 0,linetype='dashed',color='#444444')+
  theme(ggside.panel.grid=element_blank(),ggside.axis.text=element_blank(),ggside.axis.ticks=element_blank())+
  guides(size = "none",alpha='none')+labs(color='metabolic class')
dev.off()
pdf('Ext7_cor_vs_var_exp_full.pdf',w=12,h=15)
ggplot(cor.stat,aes(x=variance,y=estimate))+
  geom_point(alpha=.4)+theme_bw()+
  geom_text_repel(cor.stat,mapping=aes(label=metabolite))+
  xlab('Variance explained by origin')+ylab('Correlation between P0 and P6')
dev.off()

colnames(cor.stat)[2:4]<-c('Pearson_corr','p-value','origin_variance')
source('/project/CRI/DeBerardinis_lab/lcai/CommonFunctions/beautify_df.R')
cor.stat<-cor.stat[,1:4,with=F]
cor.stat<-beautify.dt(cor.stat,pv.key = 'p-value')
fwrite(cor.stat,'TableSx_passage_cor_and_origin_var.csv')
plot.m<-function(m){
  tmp<-data.table(s.info,val=dat.wNA[,m])
  ggplot(tmp,aes(y=val,color=Tumor.Name))+geom_point(mapping=aes(x=paste0('P',Passage)),alpha=.8)+
    geom_smooth(tmp[Passage!=0],mapping=aes(x=1+Passage),se = F,method='lm')+theme_bw()+scale_color_manual(values=ori.col)+
    theme(legend.position = 'none')+
    ggtitle(m)+xlab('Passage')+ylab('log10 value')+ #,', Residuals: ',round(value[3],2)
    labs(subtitle = with(multi.var.df[metabolite==m],paste0("(Origin: ",round(value[2],2),', Passage: ',round(value[1],2),')')))
}
pdf('Fig3HIJK_select_passage_fidelity.pdf',w=10,h=4.8)
wrap_plots(c(lapply(select.m,plot.m),lapply(select.m,plot.m2)),nrow = 2)
dev.off()

# heatmap

p.col<-structure(rev(brewer.pal(7,'Reds')),names=paste0('P',0:6))
char.metab<-as.character(multi.var.df[variable=='Origin'][value>.6][['metabolite']])
metab.sub<-metab[,char.metab]
metab.sub[]<-apply(metab.sub,2,scale)
# order cluster
dd <- as.dendrogram(hclust(dist(metab.sub),method = 'ward.D2'))
dd.reorder <- reorder(dd, ifelse(s.info$Passage==0,-20,match(s.info$Tumor.Name,all.ori)))
pdf('Ext6_All_passage_characteristic_metabolite_heatmap.pdf',w=11,h=8.2)
draw(Heatmap(t(metab.sub),show_column_names = F,cluster_columns = dd.reorder,
             clustering_method_rows = 'ward.D2',
             name='z score',row_names_gp = gpar(fontsize=8),
             right_annotation = HeatmapAnnotation(which = 'row',`P0 P6 correlation`=cor.stat[match(colnames(metab.sub),metabolite)][['Pearson_corr']],simple_anno_size = unit(.2, "cm"),
                                                  col=list(`P0 P6 correlation`=circlize::colorRamp2(c(-1,0,1), colors = brewer.pal(7,'RdBu')[c(7,4,1)])),show_annotation_name = F),
             top_annotation = HeatmapAnnotation(which='column',passage=paste0('P',s.info$Passage),annotation_legend_param = list(origin = list(at = c(all.ori,'other'))),
                                                origin=ori.map,col = list(passage=p.col,origin=ori.col2),simple_anno_size = unit(.1, "cm"))),merge_legend=T)
dev.off()
