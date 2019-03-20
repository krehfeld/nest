\name{generate_ar1}
\alias{generate_ar1}
\alias{generate_powlaw}
\alias{generate_ar1sins}

\title{Simulate irregular time series}
\description{Simulate ensembles of irregular time series with pre-defined properties.}
\usage{
generate_ar1(TX, phi,avg=TRUE,sd.in=1,osamp=10)
generate_powlaw(TX, beta = 1, SNR = 1, avg=TRUE, sd.in = 1)
generate_ar1sins(TX, phi, ps = NULL, vars = 1, SNR = 1)
}

\arguments{
  \item{TX}{Sampling times given as numerical vector or list of numerical vectors}
  \item{phi}{Number between 0 and 1, giving the lag-1 autocorrelation for AR-1 processes.}
  \item{avg}{If true, sampling by averaging, else point-sampling of the high-resolution process.}
  \item{sd.in}{Standard deviation for the output.}
  \item{osamp}{Oversampling factor for generation.}
  \item{beta}{Number giving the power law exponent. 0 for white noise. For beta>1 the time series is not stationary.}

  \item{SNR}{Signal-to-noise-ratio: Number between 0 and 1. For SNR=0 the signal is white noise.}
  \item{ps}{Vector with N periods for sinusoids}
  \item{vars}{Vector with N+1 variances (that should add up to one), where the first one is for the AR1 component}
}

\value{Returns a numerical vector (if TX is a vector) or a list of numerical vectors (if TX is a list of sampling times).}
\references{Rehfeld, K. and Kurths, J.: Similarity estimators for irregular and age-uncertain time series, Clim. Past, 10, 107-122, doi:10.5194/cp-10-107-2014, 2014.
Rehfeld, K., Marwan, N., Heitzig, J., and Kurths, J.: Comparison of correlation analysis techniques for irregularly sampled time series, Nonlin. Processes Geophys., 18, 389-404, doi:10.5194/npg-18-389-2011, 2011.
}

\author{Kira Rehfeld}


\seealso{
\code{\link{generate_t}}, \code{\link{car}}
}
\examples{
#TX<-replicate(1,generate_t(dt=1,tmin=0,tmax=100,method="gamma"),simplify=FALSE)
#TY<-replicate(1,generate_t(dt=1,tmin=0,tmax=100,method="gamma"),simplify=FALSE)
#TZ<-replicate(1,generate_t(dt=1,tmin=0,tmax=100,method="gamma"),simplify=FALSE)
TX<-generate_t(dt=1,tmin=0,tmax=100,method="gamma");
TY<-generate_t()
TZ<-generate_t(method="linear")
#print(generate_t())
X<-generate_ar1(TX=list(generate_t()),phi=0.5)


Y<-generate_powlaw(list(generate_t()),beta=1,SNR=0.5)
#print(str(Y))
#Z<-generate_ar1sins(TX=list(TX), phi=0.5, ps = c(10,30), vars = c(.5,.25,.25), SNR = 1)
#print((Y))
## Get their spectra with significance
#nestobj.X<-nestspec.sig(data=X,W=10)
#nestobj.Y<-nestspec.sig(data=Y[[1]],W=10)
#nestobj.Z<-nestspec.sig(data=Z[[1]],W=10)

## Plot the ts/ spectra
#x11(width=10,height=4)
#par(mfrow=c(2,3))
#plot(X[[1]],xlab="time")
#plot(Y[[1]],xlab="time")
#plot(Z[[1]],xlab="time")
#nestspec.plot(nestobj.X,addtitle=FALSE,speconly=TRUE); title(main="AR1")
#nestspec.plot(nestobj.Y,addtitle=FALSE,speconly=TRUE); title(main="Power law")
#nestspec.plot(nestobj.Z,addtitle=FALSE,speconly=TRUE); title(main="AR1+Sinusoid+White noise")
#abline(v=1/c(10,30))# marking the frequencies

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{irregular time series}
