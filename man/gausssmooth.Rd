\name{gausssmooth}
\alias{gausssmooth}

\title{Gaussian kernel smoother}
\description{Smoothing and detrending of arbitrarily sampled time series with a Gaussian smoother.}
\usage{gausssmooth(tx = 1:length(x), x, h = 0.5)}

\arguments{
  \item{tx}{Observation time vector}
  \item{x}{Observations}
  \item{h}{Kernel width in time units.}
}

\value{
A list with
\item{Xsmooth }{Smoothed time series}
\item{detr}{Detrended time series}
}
\references{
Rehfeld, K., Marwan, N., Heitzig, J., and Kurths, J.: Comparison of correlation analysis techniques for irregularly sampled time series, Nonlin. Processes Geophys., 18, 389-404, doi:10.5194/npg-18-389-2011, 2011.
}
\author{Kira Rehfeld}
\note{Legacy. Use gaussdetr or gaussbandpass now. Endpoints are currently not treated.}

\seealso{\code{\link{gaussdetr}},\code{\link{gaussbandpass}}}
\examples{
tx<-generate_t()
x<-sin(2*pi*tx/100)+sin(2*pi*tx/10)+rnorm(length(tx))
plot(tx,x)
lines(tx,gausssmooth(tx=tx,x=x,h=100/6)$Xsmooth,col="red")
lines(tx,gausssmooth(tx=tx,x=x,h=100/6)$detr,col="blue")
legend("topright",legend=c("original data","trend","residual"),col=c("black","red","blue"),lty=1)
}

