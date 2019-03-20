\name{car}
\alias{car}

\title{Coupled AR-1 processes}
\description{Generate coupled AR-1 processes on arbitrary time axes}
\usage{
car(tx, ty, coupl_strength, phi, lag, nsur = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tx}{time axis for time series 1 (vector)}
  \item{ty}{time axis for time series 2}
  \item{coupl_strength}{Double value between \code{-1} and \code{1}}
  \item{phi}{Autocorrelation of the time series 1}
  \item{lag}{Time lag for the coupling}
  \item{nsur}{Number of realizations to be generated}
}

\value{List of time series}
\references{Rehfeld, K., Marwan, N., Heitzig, J. and Kurths, J. (2011)
Comparison of correlation analysis techniques for irregularly sampled time series,
Nonlinear Processes in Geophysics, 18 (3), pp. 389-404. doi:10.5194/npg-18-389-2011}
\author{Kira Rehfeld}


\seealso{
\code{\link{generate_t}},\code{\link{ar1sur}},
}
\examples{
tx<-generate_t(dt=1,tmin=0,tmax=100,method="gamma")
ty<-generate_t(dt=1,tmin=0,tmax=100,method="gamma")
Proc<-car(tx,ty,coupl_strength=0.7,phi=0.5,lag=0,nsur=1)
x<-zoo(Proc$x,order.by=tx)
y<-zoo(Proc$y,order.by=ty)
phi.est=nexcf(x,lag=1,h=0.25)
couplstr.est=nexcf(x=x,y=y,lag=0,h=0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
