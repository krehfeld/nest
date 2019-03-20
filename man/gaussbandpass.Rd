\name{gaussbandpass}
\alias{gaussbandpass}
\alias{gaussdetr}
\title{Gaussian kernel bandpass and smoothing}
\description{Bandpassing, smoothing and detrending of irregular time series using a gaussian kernel smoother.
}
\usage{
gaussbandpass(X, per1, per2,prune=FALSE)
gaussdetr(X, tsc.in = mean(diff(index(X)))*10,prune=FALSE) 
}
\arguments{
  \item{X}{Input time series (\code{zoo}-object)}
  \item{per1}{Timescale 1 for lowpass}
  \item{per2}{Timescale 2 for highpass}
  \item{prune}{Logical; prune \code{(per1/mean(diff(index(X))))} points at the edges of the time series}	
  \item{tsc.in}{Timescale for detrending in \code{gaussdetr}}
	
}
\value{
A list consisting of
\item{trend}{Trend (used for highpass)}
\item{smoothed}{Smoothed time series, used for lowpass}
\item{filt}{Filtered time series (filtered=smoothed-lowpass}
}
\references{Rehfeld, K., Marwan, N., Heitzig, J. and Kurths, J. (2011)
Comparison of correlation analysis techniques for irregularly sampled time series,
Nonlinear Processes in Geophysics, 18 (3), pp. 389-404. doi:10.5194/npg-18-389-2011}

\author{Kira Rehfeld krehfeld@awi.de}



\section{Warning}{No elaborate edge-treatment yet}
\seealso{\code{\link{tsc_dep_var}}}
\examples{

## Generate one gamma-distributed and one regular time axis
tx<-generate_t(dt=1,tmin=0,tmax=100,method="gamma")
ty<-generate_t(dt=1,tmin=0,tmax=100,method="linear")
## Simulate one coupled AR1 process (see reference for details)
Proc<-car(tx,ty,coupl_strength=0.5,phi=0.5,lag=0,nsur=1)
## Bind the results to zoo time series
x<-zoo(Proc$x,order.by=tx)
y<-zoo(Proc$y,order.by=ty)


plot(x)
lines(gaussbandpass(x,10,50)$trend,lwd=2)
lines(gaussbandpass(x,10,50)$smoothed,lwd=2,col="limegreen")
lines(gaussbandpass(x,10,50)$filt,col="red2",lwd=2)
colrs<-c("black","black","limegreen","red")
legend("bottom",c("original ts","trend","smoothed","filtered"),col=colrs,lty=1,lwd=c(1,2,2,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
