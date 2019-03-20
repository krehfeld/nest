#' Generate irregular time axes with controlled properties
#' 
#' Generates regular and irregular time axes with defined start, end, and
#' sampling interval distribution.
#' 
#' 
#' @usage generate_t(dt = 1, tmin = 0, tmax = 100, method = "gamma", skew = 1)
#' @param dt Number giving the average time step
#' @param tmin Number giving the minimum time step
#' @param tmax Number giving the maximum time step
#' @param method Sampling \code{method}: "linear" for irregular or "gamma" for
#' regular
#' @param skew Number larger than zero giving the skewness of the
#' gamma-distribution.
#' @return Numeric vector containing monotonically increasing sampling times.
#' @author Kira Rehfeld
#' @seealso
#' \code{\link{generate_ar1}},\code{\link{generate_powlaw}},\code{\link{generate_ar1sins}}
#' @references Rehfeld, K., Marwan, N., Heitzig, J., and Kurths, J.: Comparison
#' of correlation analysis techniques for irregularly sampled time series,
#' Nonlin. Processes Geophys., 18, 389-404, doi:10.5194/npg-18-389-2011, 2011.
#' @keywords irregular time series
#' @examples
#' 
#' 
#' ## Generate one gamma-distributed and one regular time axis
#' tx<-generate_t(dt=1,tmin=0,tmax=100,method="gamma")
#' ty<-generate_t(dt=1,tmin=0,tmax=100,method="linear")
#' ## Simulate one coupled AR1 process (see reference for details)
#' Proc<-car(tx,ty,coupl_strength=0.7,phi=0.5,lag=0,nsur=1)
#' ## Bind the results to zoo time series
#' x<-zoo(Proc$x,order.by=tx)
#' y<-zoo(Proc$y,order.by=ty)
#' ## Estimate the autocorrelation
#' phi.est=nexcf(x,lag=1,h=0.25)
#' ## Estimate the cross-correlation
#' couplstr.est=nexcf(x=x,y=y,lag=0,h=0.05)
#' 
#' @export generate_t
generate_t <-
function(dt=1,tmin=0,tmax=100,method="gamma",skew=1){
# generate arbitrary time axes for time series
# function call:
# tx<-generate_t(dt,tmin,tmax,method) 
# methods: linear/gamma sampling - the parameter skew determines the irregularity of the sampling (cf. Rehfeld et al.,NPG, 2011)
# 20.12.13 krehfeld@awi.de
#
# linear
if((method=="linear")||(skew==0)){
temp=seq(tmin,tmax,dt)
} 
else{ #method=="gamma"
#Gamma-distributed sampling intervals are generated. The irregularity is determined by the skew variable.
shape=1/skew;
scale=dt/shape;
# generate more random variables than necessary
L=2*ceiling((tmax-tmin)/dt)
R=rgamma(L,shape,scale=scale)
# kick out timesteps that are smaller than an arbitrary precision (some interpolation routines complain if the timesteps are too small)
R=R[R>0.0000001];
# cumulative sum
temp=cumsum(R);

dtin<-mean(diff(temp));
#print(dtin)

#temp<-temp/dtin*dt
#select times in window
#tx=temp[(tmin*0.95<=temp)&(temp<tmax*1.05)];
# if(any(temp)>tmax){
# 
# }



}


tx=temp[1:(floor((tmax-tmin)/dt)+1)] # fix the length of the records
#dtin<-mean(diff(tx));
#tx<-tx/dtin*dt
#if (any(tx>tmax)||any(tx<tmin)){
#tx<-rescale(tx,to=c(tmin,tmax))
#tx<-scales:::rescale(tx,to=c(tmin,tmax))

#}
#if (){tx<-tx-min(tx)+tmin}

#tx=tx[(tmin<=tx)&(tx<tmax)];

if (any(is.na(tx))){warning("NA values");}

return(tx)

}
