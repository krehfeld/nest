\name{tsc_dep_var}
\alias{tsc_dep_var}
\title{Timescale-dependent variance estimation}
\description{Estimate the timescale-dependent variance of a time series.}
\usage{tsc_dep_var(timser, tsc.in, min.res =
                 min(mean(diff(index(timser))), floor(tsc.in[1]/2) -
                 1), start.val = start(timser), end.val = end(timser),
                 pval = 0.1, detrend = FALSE)
}

\arguments{
  \item{timser}{Input time series. Can be a \code{zoo}, \code{ts} or vector object.}
  \item{tsc.in}{Vector of two timescales, \code{c(tsc1,tsc2)}.}
  \item{min.res}{Minimal resolution of the time series (if this is not met, NAs are returned).}
  \item{start.val}{Starting value of the time series window.}
  \item{end.val}{End value of the time series window.}
  \item{pval}{\code{pval} gives the p-value for the confidence intervals of the variance estimate.}
  \item{detrend}{logical argument, if set to \code{TRUE} the timesereis is linearly detrended prior to spectral analysis}
}

\value{
A list of
\item{std}{Standard deviation}
\item{var.tsc}{timescale-dependent variance estimated}
\item{var.tot}{total variance of the time series in the window}
\item{dof}{estimated degrees of freedom of the variance estimate}
\item{ts.used}{used time series}
\item{var.ci}{list of \code{lo} and \code{up}, lower and upper confidence level for the variance estimate}
\item{spec.win}{spectrum object}
}

\author{
Kira Rehfeld, with contributions from Thom Laepple
}




\seealso{\code{\link{var_ci}}}

\examples{
library(scales)
## Generate one gamma-distributed and one regular time axis
tx<-generate_t(dt=1,tmin=0,tmax=250,method="gamma")
ty<-generate_t(dt=1,tmin=0,tmax=250,method="linear")
## Simulate one coupled AR1 process (see reference for details)
Proc<-car(tx,ty,coupl_strength=0.5,phi=0.5,lag=0,nsur=1)
## Bind the results to zoo time series
x<-zoo(Proc$x,order.by=tx)
y<-zoo(Proc$y,order.by=ty)


tsc_dep_var(y,tsc.in=c(30,100)) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
