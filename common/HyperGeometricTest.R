#####This is a function for gene set enrichment analysis#####
#default setting
options(stringsAsFactors = F)

#load gene sets
try(load('/project/CRI/DeBerardinis_lab/lcai/CommonData/c2cp.RData'),silent=T)
try(load('Z:/CommonData/c2cp.RData'),silent=T)

#load('/project/CRI/DeBerardinis_lab/lcai/CommonData/c3mir.RData')
#load('/project/CRI/DeBerardinis_lab/lcai/CommonData/c2cgp.RData')
try(load('/archive/CRI/DeBerardinis_lab/lcai/TCGA_NeighborProject/Figures/ProteinCorrelation/RData/GeneFamily.RData'),silent = T)
try(load('Z:/TCGA_NeighborProject/Figures/ProteinCorrelation/RData/GeneFamily.RData'),silent = T)
rm(GeneFamilyByID)

##Original hypergeometric method
#N: Number of all genes
#B: Size of geneset (within all genes)
#b: overlapping genes between selected genes and gene set
#select: number of selected genes
#phyper(b,B,N-B,select,lower.tail=F)

#problem: it is always easier for larger gene sets to get more
#significant p-value

#To caluclate null distribution for the test
#geneset<-sample(select,N)

###tuning function, can be deleted later
LC.hyper<-function(genesets,genes,testset){
  N<-length(genes)
  p.set<-vector('numeric',length(genesets))
  names(p.set)<-names(genesets)
  for(i in 1:length(genesets)){
    #print(i)
    geneset<-genesets[[i]]
    geneset<-intersect(genes,geneset)
    if(length(geneset)<5) {
      p.set[i]<-NA
      next
    }
    B<-length(geneset)
    b<-length(intersect(geneset,testset))
    select<-length(testset)
    p.set[i]<-phyper(b-1,B,N-B,select,lower.tail=F)
  }
  return(p.set)
}
# with enrichment ratio
LC.hyper2<-function(genesets,genes,testset,include.intersection=F){
  N<-length(genes)
  overlapping.genes<-set.size<-overlap<-e.ratio<-p.set<-vector('numeric',length(genesets))
  names(p.set)<-names(genesets)
  for(i in 1:length(genesets)){
    #print(i)
    geneset<-genesets[[i]]
    geneset<-intersect(genes,geneset)
    if(length(geneset)<5) {
      p.set[i]<-NA
      next
    }
    B<-length(geneset)
    b<-length(intersect(geneset,testset))
    select<-length(testset)
    p.set[i]<-phyper(b-1,B,N-B,select,lower.tail=F)
    overlapping.genes[i]<-paste(intersect(geneset,testset),collapse = ',')
    set.size[i]<-B;overlap[i]<-b;e.ratio[i]<-b/(B*select/N)
  }
  if(include.intersection){
    data.table(feature=names(p.set),pv=p.set,enrichment.ratio=e.ratio,set.size,overlap,overlapping.genes)
  }else{
    data.table(feature=names(p.set),pv=p.set,enrichment.ratio=e.ratio,set.size,overlap)
  }
}

###Examples:
#p.RNA<-LC.hyper(genesets=c2cp,
#                genes=names(r.RNA),
#                testset=names(sort(r.RNA,decreasing=T))[1:50])
#p.RES<-LC.hyper(genesets=c2cp,
#                genes=names(r.RNA),
#                testset=names(sort(r.RES,decreasing=T))[1:50])
#head(sort(p.RNA))
#temp<-data.frame(p.RNA,p.RES)
#head(temp[order(p.RES,decreasing=F),])
#head(temp[order(p.RNA,decreasing=F),])
