#' Gaussian kernel bandpass and smoothing
#' 
#' Bandpassing, smoothing and detrending of irregular time series using a
#' gaussian kernel smoother.
#' 
#' 
#' @aliases gaussbandpass gaussdetr
#' @usage gaussbandpass(X, per1, per2,prune=FALSE) gaussdetr(X, tsc.in =
#' mean(diff(index(X)))*10,prune=FALSE)
#' @param X Input time series (\code{zoo}-object)
#' @param per1 Timescale 1 for lowpass
#' @param per2 Timescale 2 for highpass
#' @param prune Logical; prune \code{(per1/mean(diff(index(X))))} points at the
#' edges of the time series
#' @param tsc.in Timescale for detrending in \code{gaussdetr}
#' @return A list consisting of \item{trend}{Trend (used for highpass)}
#' \item{smoothed}{Smoothed time series, used for lowpass} \item{filt}{Filtered
#' time series (filtered=smoothed-lowpass}
#' @section Warning: No elaborate edge-treatment yet
#' @author Kira Rehfeld krehfeld@@awi.de
#' @seealso \code{\link{tsc_dep_var}}
#' @references Rehfeld, K., Marwan, N., Heitzig, J. and Kurths, J. (2011)
#' Comparison of correlation analysis techniques for irregularly sampled time
#' series, Nonlinear Processes in Geophysics, 18 (3), pp. 389-404.
#' doi:10.5194/npg-18-389-2011
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' 
#' ## Generate one gamma-distributed and one regular time axis
#' tx<-generate_t(dt=1,tmin=0,tmax=100,method="gamma")
#' ty<-generate_t(dt=1,tmin=0,tmax=100,method="linear")
#' ## Simulate one coupled AR1 process (see reference for details)
#' Proc<-car(tx,ty,coupl_strength=0.5,phi=0.5,lag=0,nsur=1)
#' ## Bind the results to zoo time series
#' x<-zoo(Proc$x,order.by=tx)
#' y<-zoo(Proc$y,order.by=ty)
#' 
#' 
#' plot(x)
#' lines(gaussbandpass(x,10,50)$trend,lwd=2)
#' lines(gaussbandpass(x,10,50)$smoothed,lwd=2,col="limegreen")
#' lines(gaussbandpass(x,10,50)$filt,col="red2",lwd=2)
#' colrs<-c("black","black","limegreen","red")
#' legend("bottom",c("original ts","trend","smoothed","filtered"),col=colrs,lty=1,lwd=c(1,2,2,2))
#' 
#' @export gaussbandpass
gaussbandpass<-function(X,per1,per2,prune=FALSE){
if (!(per1<per2)){stop("smoothing timescale per1 must be smaller than detrending timescale per2")}

tsx.trend<-gaussdetr(X,tsc.in=per2)$Xsmooth


tsx.smoo<-gaussdetr(X,tsc.in=per1)$Xsmooth

DT<-mean(diff(index(X)))

if (prune){
tsx.smoo[c(1:ceiling(per1/(DT)),length(tsx.smoo):(length(tsx.smoo)-ceiling(per1/(DT))))]<-NA
}

tsx.filt<-tsx.smoo-tsx.trend
return(list(trend=na.omit(tsx.trend),smoothed=na.omit(tsx.smoo),filt=na.omit(tsx.filt)))
}

