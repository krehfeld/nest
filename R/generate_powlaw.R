generate_powlaw <-
function(TX,beta=1,SNR=1,avg=TRUE,sd.in=1){
#Generate power law time series, with a log-log-slope of the spectrum equal to beta

X<-list()


if (is.list(TX)){
nsur=length(TX)
 tl<-TX
} else {nsur=1 
tl<-list(TX)}

#browser()
for (n in 1:nsur){
####<----------------------------------------------------------- beginning of loop

tx<-as.numeric(tl[[n]])
dtfactor=10
dtsur=min(mean(diff(tx)))/dtfactor;
t1=min(tx); # subtract lag  to compensate transient effects;
t2=max(tx); #add lag to avoid extrapolation at the end
#make interpolation timescale

t=seq(t1,t2+2*dtsur,dtsur);


#persistence time rescaled for higher time series resolution
#fi=exp(-dtsur/(-1/log(phi)));

# First step: Create driving process at high resolution
xsur=sim.powerlaw(beta,length(t))

# interpolate linearly to the original timescale desired
#filter.pTs1 <- function(data,filter,method=1,...)
# pad the time series with long-term mean (see filter.R)
#filter.pTs1(xsur,c(1,1,1),method=1)
#zoohres<-zoo(xsur,order.by=t)

xsur<-(xsur-mean(xsur))/sd(xsur)

## Add white noise according to SNR
#print(SNR)


if (!SNR==1){
wn<-rnorm(length(t),mean=0,sd=1)
wn<-(wn-mean(wn))/sd(wn)
xsur<-SNR*xsur+(1-SNR)*wn
#x<-(x-mean(x))/sd(x) # normalize to 0 mean, unit variance; again [to have overall var(x)==1]
}
xsur<-(xsur-mean(xsur))/sd(xsur)*sd.in



zoohres<-zoo(xsur,order.by=t)

#xsursmoo=gausssmooth(tx=t,x=xsur,h=mean(diff(tx))/6) #<--- reduces the overall variance
if (!avg){
x=approx(index(zoohres),coredata(zoohres),tx,method="linear",rule=1)$y

} else {

x<-zoo(subsample(zoohres,dtx=diff(tx),tx),order.by=tx)
}





xts<-(zoo(x,order.by=tx))
X[[n]]<-xts

rm(x,xts)
### <------------------------------------------- end of loop
}


return(X)

}
