# works with data.table
beautify.dt<-function(dt,pv.key=NULL){
  if(is.data.frame(dt) & !is.data.table(dt)) {df<-T;dt.rname<-rownames(dt);setDT(dt)} else df<-F
  num.col<-which(apply(as.matrix(dt),2,FUN=function(x) identical(is.na(x),is.na(as.numeric(x)))))
  int.col<-which(apply(as.matrix(dt),2,FUN=function(x) identical(as.numeric(x),as.numeric(as.integer(x)))))
  pv.cols<-num.col[grep(paste(c('pv','padj',pv.key),collapse ='|'),colnames(dt)[num.col])]
  for(pv.col in pv.cols) dt[[pv.col]]<-signif(dt[[pv.col]],2)
  for(other.col in setdiff(num.col,c(pv.cols,int.col))) dt[[other.col]]<-round(dt[[other.col]],2)
  if(df){dt<-as.data.frame(dt);rownames(dt)<-dt.rname}
  dt
}
