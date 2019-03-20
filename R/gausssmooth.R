#' Gaussian kernel smoother
#' 
#' Smoothing and detrending of arbitrarily sampled time series with a Gaussian
#' smoother.
#' 
#' 
#' @usage gausssmooth(tx = 1:length(x), x, h = 0.5)
#' @param tx Observation time vector
#' @param x Observations
#' @param h Kernel width in time units.
#' @return A list with \item{Xsmooth }{Smoothed time series}
#' \item{detr}{Detrended time series}
#' @note Legacy. Use gaussdetr or gaussbandpass now. Endpoints are currently
#' not treated.
#' @author Kira Rehfeld
#' @seealso \code{\link{gaussdetr}},\code{\link{gaussbandpass}}
#' @references Rehfeld, K., Marwan, N., Heitzig, J., and Kurths, J.: Comparison
#' of correlation analysis techniques for irregularly sampled time series,
#' Nonlin. Processes Geophys., 18, 389-404, doi:10.5194/npg-18-389-2011, 2011.
#' @examples
#' 
#' tx<-generate_t()
#' x<-sin(2*pi*tx/100)+sin(2*pi*tx/10)+rnorm(length(tx))
#' plot(tx,x)
#' lines(tx,gausssmooth(tx=tx,x=x,h=100/6)$Xsmooth,col="red")
#' lines(tx,gausssmooth(tx=tx,x=x,h=100/6)$detr,col="blue")
#' legend("topright",legend=c("original data","trend","residual"),col=c("black","red","blue"),lty=1)
#' 
#' @export gausssmooth
gausssmooth <-
function(tx=1:length(x),x,h=0.5){
#GAUSSSMOOTH smoothes irregularly sampled time series with gaussian weights.
#   [Xsmooth] = GAUSSSMOOTH(tx,x,h) computes a smoothed version of
#   the time series observations x, where the amount of smoothing is
#   determined by the kernel width h. The first points should be taken with
#   caution as they are only smoothed onesidedly.
#
# Function calls:
# est<-gausssmooth(ty,y+trend,10) # use kernel width of 10 (high smoothing)
# est<-gausssmooth(ty,y+trend,) # use default kernel width of 0.5 (low smoothing)
# est<-gausssmooth(,y+trend,) # use only if ty is regular, and the default kernel 
#
# width of 0.5 should be # applied to the index of the time series
#    REFERENCE:
#       K. Rehfeld, N. Marwan, J. Heitzig, J. Kurths: Comparison of
#       correlation analysis techniques for irregularly sampled time series,
#       Nonlinear Processes in Geophysics, 18(3), 389-404 (2011).
# (c) Alfred-Wegner Institute for Polar and Marine Research, 2013
# Kira Rehfeld
# ver: 0.1  last rev: 2013-12-29
#
# Note: The function was ported from the Octave function in The NESToolbox. To cater for irregular time series, the observation time is considered in the vector, the proxy observations in vector x, here. In many applications in R, regular sampled time series are stored handily as time series objects, but this is not a feature that can be used with this function - yet.


if(is.ts(x)){
#print("use time series properties of the original timeseries!")
 Xsmooth=detr=x;
} else {
Xsmooth=rnorm(length(tx));detr<-Xsmooth;
}

Xsmooth[1]=mean(x[1:floor(h/2)]);

for (k in 2:length(tx)) {
WL=1/sqrt(2*pi*h^2)*(exp(-((tx-tx[k])^2)/(2*h^2)));
Xsmooth[k]=sum(WL*x)/sum(WL); 
}
   
Xsmooth[1]=Xsmooth[2];

#detr<-c(x)-c(Xsmooth);
#plot(ts(x)-Xsmooth)
detr<-x-Xsmooth;
result<-list(Xsmooth=Xsmooth,detr=detr)   
class(result)<-"smoothingobject" 
return(result)
#return(list(Xsmooth=Xsmooth,detr=detr))
}
