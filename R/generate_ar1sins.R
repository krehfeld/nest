generate_ar1sins <-
function(TX,phi,ps=NULL,vars=1,SNR=1){
#Generate AR1 data, with temporal sampling as in tx and an autocorrelation coefficient of phi. Then add sinusoids.
 

# set the resolution at which the AR1process is simulated at high
# resolution as 10 times the actual mean resolution
#tl<-list()

#ps-> vector with n_p periods
#vars-> vector with n_p+1 variances (that should add up to one), where the first one is for the ar1
#snr -> signal-to-noise ratio: equal to one = no additional white noise
if (sum(vars)>1) stop("Variances don't add up to one")

X<-list()


if (is.list(TX)) {nsur=length(TX); tl<-TX;} else {nsur=1; tl<-list(TX)}


for (n in 1:nsur){

tx<-as.numeric(tl[[n]])

dtsur=min(mean(diff(tx)))/10;
t1=min(tx); # subtract lag  to compensate transient effects;
t2=max(tx); #add lag to avoid extrapolation at the end
#make interpolation timescale

t=seq(t1,t2+2*dtsur,dtsur);


#persistence time rescaled for higher time series resolution
fi=exp(-dtsur/(-1/log(phi)));


####<----------------------------------------------------------- beginning of loop

# First step: Create driving process at high resolution
eps1=rnorm(length(t),0,1);
xsur=eps1;
for(i in 2:length(t)){
xsur[i]=fi*xsur[i-1] +sqrt(1-fi^2)*eps1[i];
# do something
}

xsur=(xsur-mean(xsur))/sd(xsur); #zero mean unit variance for ar


# interpolate linearly to the original timescale desired
intoutx=approx(t,xsur,tx,method="linear",rule=1)


# assign output of ar1
x<-intoutx$y;




Sx<-matrix(NA,nrow=length(tx),ncol=(length(ps)+1))
Sx[,1]<-cbind(x) # ar1 in first col

# simulate sinusoids on the time axis of time series X[[n]] 
if (is.vector(ps)){ # check if periodicities are requested
# generate random phases for the sinusoids
randphases<-matrix(runif(n=length(vars),min=0,max=pi),nrow=length(vars))

for (m in 2:(length(ps)+1)){

Sx[,m]<-scale(sin(2*pi*tx/ps[m-1]+randphases[m-1]))


}
# do variance mixing



}

Sx<-scale(Sx) # standardize all subsignals to unit mean and zero variance
Xin<-apply(sqrt(vars)*t(Sx),2,sum,na.omit=TRUE) # mix the signals


### observational noise

obsnoise<-scale(matrix(rnorm(length(tx)),nrow=length(tx),ncol=1))

Xout<-sqrt(SNR)*Xin+sqrt(1-SNR)*obsnoise



### end of loop

xts<-zoo(Xout,order.by=tx)

## normalize
xts<-(xts-mean(xts))/sd(xts)

X[[n]]<-xts

rm(x,xts,Xout,obsnoise,Xin,Sx)

}


return(X)

}
