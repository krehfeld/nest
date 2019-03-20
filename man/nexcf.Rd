\name{nexcf}
\alias{nexcf}
\alias{ixcf}
\alias{nexcf_ci}

\title{Pearson correlation for irregular time series}
\description{Compute the correlation for irregular time series with or without interpolation.}
\usage{nexcf(x, y = x, lag = 0, h = 0.25, enforce = FALSE)
       ixcf(x, y = x, lag = 0)
       nexcf_ci(x, y, conflevel = 0.05, dt = NULL, smoo = FALSE,detr=FALSE)}

\arguments{
  \item{x}{time series one to be analyzed. Note that time increments should be positive (i.e. time is reported as age before present).}
  \item{y}{time series two to be analyzed (if omitted, autocorrelation is computed)}
  \item{lag}{numeric vector giving lags for which correlation is to be estimated}
  \item{h}{number giving the Gaussian kernel width}
  \item{enforce}{Coerce the result to the interval (-1,1) in case of numerical problems}
  \item{conflevel}{Confidence-level for the two-sided confidence intervals around zero correlation}
  \item{dt}{Average time step for the correlation function evaluation (by default set to max(mean(dtx),mean(dty)}
  \item{smoo}{logical; set to TRUE for smoothing of variability below the time step}
\item{detr}{logical; set to TRUE for linear detrending of the time series}
}

\value{Numerical vector of the length of the lag vector.}
\references{Rehfeld, K. and Kurths, J.: Similarity estimators for irregular and age-uncertain time series, Clim. Past, 10, 107-122, doi:10.5194/cp-10-107-2014, 2014.
Rehfeld, K., Marwan, N., Heitzig, J., and Kurths, J.: Comparison of correlation analysis techniques for irregularly sampled time series, Nonlin. Processes Geophys., 18, 389-404, doi:10.5194/npg-18-389-2011, 2011.
}
\author{Kira Rehfeld}


\seealso{
\code{\link{generate_t}}
}
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
## Estimate the cross-correlation with interpolation
couplstr.est.interp=ixcf(x=x,y=y,lag=0)

## Estimate the effective degrees of freedom, and the significance of the correlation estimate,
## from the t-distribution.
temp<-nexcf_ci(x,y,conflevel=0.1)
# degrees of freedom used for the significance test
print(temp$neff)
print(temp$rxy)
# if rxy is outside 
temp$ci
# the confidence intervals, the correlation estimate is significant at the (1-conflevel)*100\% level
}

\keyword{correlation}
\keyword{irregular time series}
