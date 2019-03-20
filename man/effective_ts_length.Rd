\name{effective_ts_length}
\alias{effective_ts_length}
\title{degrees of freedom of an irregular autocorrelated time series.}
\description{Get the effective length, or degrees of freedom, of an irregular autocorrelated time series, using different methods.}
\usage{effective_ts_length(X, method = "tau")}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A zoo-time-series object.}
  \item{method}{One of "tau", "mudelsee", or "simple"}
}

\value{Approximated number of degrees of freedom, estimated using the autocorrelation of the time series.}
\references{
von Storch, H., and Zwiers, F.W., 1999. Statistical Analysis in Climate Research. Cambridge University Press [Page number to be checked]
Dawdy, D.R., and Matalas, N.C., 1964, Statistical and probability analysis of hydrologic data,part III: Analysis of variance, covariance and time series, in Ven Te Chow, ed., Handbook of applied hydrology, a compendium of water-resources technology: New York, McGraw-Hill Book Company, p. 8.68-8.90.
Mudelsee,M., 2010. Climate Time Series Analysis - Classical Statistical and Bootstrap Methods. Springer ISBN: 9789048194827, 474 p [Page number for reference to be checked] 
}
\author{Kira Rehfeld, with contribution from Thom Laepple.}


\seealso{\code{\link{nexcf_ci}}}
\examples{

## Generate one gamma-distributed and one regular time axis
tx<-generate_t(dt=1,tmin=0,tmax=100,method="gamma")
ty<-generate_t(dt=1,tmin=0,tmax=100,method="linear")
## Simulate one coupled AR1 process (see reference for details)
Proc<-car(tx,ty,coupl_strength=0.5,phi=0.5,lag=0,nsur=1)
## Bind the results to zoo time series
x<-zoo(Proc$x,order.by=tx)
y<-zoo(Proc$y,order.by=ty)
## Estimate the autocorrelation
phi.est=nexcf(x,lag=1,h=0.25)

plot(x)
lines(y,col="red")

rxy=nexcf(x,y,lag=0,h=0.25)
rx1=nexcf(x,lag=1)
ry1=nexcf(y,lag=1)

nexcf_ci(x,y)


effective_ts_length(x,method="tau") # adapted from von Storch/Zwiers 1999
effective_ts_length(x,method="mudelsee") # after Mudelsee, 2010
effective_ts_length(x,method="simple") # using both time directions
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
