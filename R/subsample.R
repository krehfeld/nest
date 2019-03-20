subsample<-function(zoohres,dtx=diff(tx),tx,avg=TRUE,fast=FALSE){
# subsampling a time series by taking means or by taking point samples
#subsample(zoo(rnorm(100),generate_t()),tx=seq(5,75,by=5))
#subsample(zoo(rnorm(100),generate_t()),dtx=rep(3,length.out=length(diff(tx))),tx=seq(5,95,by=10))
#? How can I increase the speed? by splitting the start/end cutting from the mean operation?
xnew<-zoo(rep(NA,length(tx)),order.by=tx)
# sampling integration of a high resolution time series
if (avg){
startvals<-tx[1:(length(tx))]-c(dtx[1],dtx)/2
endvals<-tx[1:(length(tx))]+c(dtx,dtx[length(dtx)])/2

if (fast){

trep<-rep(list(zoohres),length(startvals))# list of time series

tmp<-mapply(function(startval,endval,temp){tims<-unlist(temp);
return(window(tims,start=startval,end=endval))},startvals,endvals,rep(list(zoohres),length(startvals))) #cut them

xnew<-zoo(sapply(tmp,mean,na.omit=TRUE),order.by=tx) # take means


#zoo(mapply(function(startval,endval,temp){tims<-unlist(temp);
#return(mean(window(tims,start=startval,end=endval),na.omit=TRUE))},startvals,endvals,rep(list(zoohres),length(startvals))),order.by=tx)
}
else {
# correct but very slow
for (n in 1:length(xnew)){xnew[n]<-mean(window(zoohres,start=startvals[n],end=endvals[n]),na.omit=TRUE)}
}
}

if (!avg){# point sampling/aliasing
xnew<-zoo(approx(index(zoohres),coredata(zoohres),tx)$y,tx)
}

if(any(is.na(xnew))) {

warning("NA values in results - fixing by interpolation")
ind<-which(is.na(xnew))
for (i in 1:length(ind)){
xnew[ind[i]]<-mean(c(coredata(xnew)[ind[i]-1],coredata(xnew)[ind[i]+1]),na.rm=TRUE)
}
}

return(xnew)
}