\name{generate_t}
\alias{generate_t}
\title{Generate irregular time axes with controlled properties}

\description{Generates regular and irregular time axes with defined start, end, and sampling interval distribution.}
\usage{generate_t(dt = 1, tmin = 0, tmax = 100, method = "gamma", skew = 1)}

\arguments{
  \item{dt}{Number giving the average time step}
  \item{tmin}{Number giving the minimum time step}
  \item{tmax}{Number giving the maximum time step}
  \item{method}{Sampling \code{method}: "linear" for irregular or "gamma" for regular}
  \item{skew}{Number larger than zero giving the skewness of the gamma-distribution.}
}

\value{Numeric vector containing monotonically increasing sampling times.}
\references{Rehfeld, K., Marwan, N., Heitzig, J., and Kurths, J.: Comparison of correlation analysis techniques for irregularly sampled time series, Nonlin. Processes Geophys., 18, 389-404, doi:10.5194/npg-18-389-2011, 2011.}
\author{Kira Rehfeld}


\seealso{\code{\link{generate_ar1}},\code{\link{generate_powlaw}},\code{\link{generate_ar1sins}} }

\examples{

## Generate one gamma-distributed and one regular time axis
tx<-generate_t(dt=1,tmin=0,tmax=100,method="gamma")
ty<-generate_t(dt=1,tmin=0,tmax=100,method="linear")
## Simulate one coupled AR1 process (see reference for details)
Proc<-car(tx,ty,coupl_strength=0.7,phi=0.5,lag=0,nsur=1)
## Bind the results to zoo time series
x<-zoo(Proc$x,order.by=tx)
y<-zoo(Proc$y,order.by=ty)
## Estimate the autocorrelation
phi.est=nexcf(x,lag=1,h=0.25)
## Estimate the cross-correlation
couplstr.est=nexcf(x=x,y=y,lag=0,h=0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{irregular time series} 	
