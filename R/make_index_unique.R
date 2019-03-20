make_index_unique<-function(zoobj){
if (length(zoobj)==0) return(zoobj)
tx<-index(zoobj)
dind<-which(diff(index(zoobj))==0)
if (length(dind)>0){
zoobj<-zoobj[-dind]
}
if (any(is.na(zoobj))){zoobj<-na.omit(zoobj)}
if (mean(diff(index(zoobj)))<0) {warning("decreasing index - check time definition")}
return(zoobj)
}