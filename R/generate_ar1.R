#' Simulate irregular time series
#' 
#' Simulate ensembles of irregular time series with pre-defined properties.
#' 
#' 
#' @aliases generate_ar1 generate_powlaw generate_ar1sins
#' @usage generate_ar1(TX, phi,avg=TRUE,sd.in=1,osamp=10) generate_powlaw(TX,
#' beta = 1, SNR = 1, avg=TRUE, sd.in = 1) generate_ar1sins(TX, phi, ps = NULL,
#' vars = 1, SNR = 1)
#' @param TX Sampling times given as numerical vector or list of numerical
#' vectors
#' @param phi Number between 0 and 1, giving the lag-1 autocorrelation for AR-1
#' processes.
#' @param avg If true, sampling by averaging, else point-sampling of the
#' high-resolution process.
#' @param sd.in Standard deviation for the output.
#' @param osamp Oversampling factor for generation.
#' @param beta Number giving the power law exponent. 0 for white noise. For
#' beta>1 the time series is not stationary.
#' @param SNR Signal-to-noise-ratio: Number between 0 and 1. For SNR=0 the
#' signal is white noise.
#' @param ps Vector with N periods for sinusoids
#' @param vars Vector with N+1 variances (that should add up to one), where the
#' first one is for the AR1 component
#' @return Returns a numerical vector (if TX is a vector) or a list of
#' numerical vectors (if TX is a list of sampling times).
#' @author Kira Rehfeld
#' @seealso \code{\link{generate_t}}, \code{\link{car}}
#' @references Rehfeld, K. and Kurths, J.: Similarity estimators for irregular
#' and age-uncertain time series, Clim. Past, 10, 107-122,
#' doi:10.5194/cp-10-107-2014, 2014. Rehfeld, K., Marwan, N., Heitzig, J., and
#' Kurths, J.: Comparison of correlation analysis techniques for irregularly
#' sampled time series, Nonlin. Processes Geophys., 18, 389-404,
#' doi:10.5194/npg-18-389-2011, 2011.
#' @keywords irregular time series
#' @examples
#' 
#' #TX<-replicate(1,generate_t(dt=1,tmin=0,tmax=100,method="gamma"),simplify=FALSE)
#' #TY<-replicate(1,generate_t(dt=1,tmin=0,tmax=100,method="gamma"),simplify=FALSE)
#' #TZ<-replicate(1,generate_t(dt=1,tmin=0,tmax=100,method="gamma"),simplify=FALSE)
#' TX<-generate_t(dt=1,tmin=0,tmax=100,method="gamma");
#' TY<-generate_t()
#' TZ<-generate_t(method="linear")
#' #print(generate_t())
#' X<-generate_ar1(TX=list(generate_t()),phi=0.5)
#' 
#' 
#' Y<-generate_powlaw(list(generate_t()),beta=1,SNR=0.5)
#' #print(str(Y))
#' #Z<-generate_ar1sins(TX=list(TX), phi=0.5, ps = c(10,30), vars = c(.5,.25,.25), SNR = 1)
#' #print((Y))
#' ## Get their spectra with significance
#' #nestobj.X<-nestspec.sig(data=X,W=10)
#' #nestobj.Y<-nestspec.sig(data=Y[[1]],W=10)
#' #nestobj.Z<-nestspec.sig(data=Z[[1]],W=10)
#' 
#' ## Plot the ts/ spectra
#' #x11(width=10,height=4)
#' #par(mfrow=c(2,3))
#' #plot(X[[1]],xlab="time")
#' #plot(Y[[1]],xlab="time")
#' #plot(Z[[1]],xlab="time")
#' #nestspec.plot(nestobj.X,addtitle=FALSE,speconly=TRUE); title(main="AR1")
#' #nestspec.plot(nestobj.Y,addtitle=FALSE,speconly=TRUE); title(main="Power law")
#' #nestspec.plot(nestobj.Z,addtitle=FALSE,speconly=TRUE); title(main="AR1+Sinusoid+White noise")
#' #abline(v=1/c(10,30))# marking the frequencies
#' 
#' 
#' @export generate_ar1
generate_ar1 <-
function(TX,phi,avg=TRUE,sd.in=1,osamp=10){
#Generate AR1 data, with temporal sampling as in tx and an autocorrelation coefficient of phi. 

# set the resolution at which the AR1process is simulated at high
# resolution as 10 times the actual mean resolution
#tl<-list()
X<-list()
tl<-list()
nsur=matrix()

if (is.list(TX)){tl<-TX; nsur<-length(TX); returnlist=TRUE} 
if (is.matrix(TX)) {nsur=dim(TX)[1]; tl<-split(TX,row(TX)); returnlist=TRUE}
if (length(TX)==1){returnlist=FALSE;nsur=1}#vector
#if (is.vector(TX)) {nsur=1; tl[[1]]<-TX; print("vector")} # das funktioniert so nicht - lists are also vectors

for (n in 1:nsur){



tx<-as.numeric(tl[[n]])

dtsur=median(diff(tx))/osamp;
t1=min(tx); # subtract lag  to compensate transient effects;ge
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
}

xsur=sd.in*(xsur-mean(xsur))/sd(xsur);

# interpolate linearly to the original timescale desired || sampling point-by-point
if (!avg){
x=approx(t,xsur,tx,method="linear",rule=2)$y
#x=(x-mean(x))/sd(x);
xts<-zoo(x,order.by=tx)
}

else {zoohres<-zoo(xsur,order.by=t); xts=subsample(zoohres,tx=tx,fast=FALSE);

# x=approx(t,xsur,tx,method="linear",rule=2)$y
# #x=(x-mean(x))/sd(x);
# xts2<-zoo(x,order.by=tx)
# 
# plot(zoohres)
# points(xts,col="blue",pch=8)
# points(xts2,col="magenta",pch=24)
# 

}
# Problem: sehr schiefe Zeitverteilungen führen zu leeren bins im Sampling -> NAs
# Lösungsansatz: Oversampling hoch setzen. 
# Nächstes Problem: Memory/ Zeit geht viel zu hoch...
if (any(is.na(xts))) print(xts)

# assign output
#x<-intoutx$y;




### end of loop

X[[n]]<-xts
rm(xts)

}



if (n==1&&returnlist==FALSE){res<-X[[1]]}
if (returnlist){res<-X}
return(res)

}
