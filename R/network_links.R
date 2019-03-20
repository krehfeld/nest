network_links<-function(WBTobj,conflevel=0.1,detr=FALSE,dt.lag=NULL){
# use this wrapper function to get similarities
f1 <- function(a,b, fun){
  outer(a, b, function(x,y) vapply(seq_along(x), function(i) fun(x[[i]], y[[i]]), numeric(1)))
}


if (!class(WBTobj)=="window-by-time-object") stop("Check arguments. Function needs a window-by-time-object created with window_by_time()")

WBTobj.tme<-WBTobj[[1]]$tmid
WBTobj.tls<-lapply(WBTobj,function(x){x$tsplit}) #list of time series only

CMAT<-PMAT<-array(NA,dim=c(length(WBTobj.tme),length(WBTobj),length(WBTobj)))
 for (i in 1:length(WBTobj.tme)){
 TSlist<-lapply(WBTobj.tls,function(x){x[[i]]}) # list with 
 CMAT[i,,]<-f1(TSlist,TSlist,function(x,y){nexcf_ci(x,y,detr=detr,dt=dt.lag)$rxy}) #lag1 cross-correlation
 PMAT[i,,]<-f1(TSlist,TSlist,function(x,y,dt=dt.lag){nexcf_ci(x,y,conflevel=conflevel,detr=detr)$pval}) #lag1 cross-correlation
 }
 results<-list(tmid=WBTobj.tme,CMAT=CMAT,PMAT=PMAT)
 class(results)<-"palnet.links"
 return(results)
 }
 
