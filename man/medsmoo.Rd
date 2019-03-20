\name{medsmoo}
\alias{medsmoo}

\title{Median Smoother}
\description{Mean and median smoothing of numerical vectors with optional edge preservation.}
\usage{medsmoo(x, m, m2, edge = c("both"))}

\arguments{
  \item{x}{numeric vector}
  \item{m}{half width of the smoother}
  \item{m2}{half width of the smoother at the edges}
  \item{edge}{"left", "right" or "both" for application of smaller windows towards the edges}
}

\value{

A list of 
\item{rmedian }{running median}
\item{rmean }{running mean}
\item{rmad }{running median absolute deviation}
\item{smootherwidth}{smoother width for each point}
\item{edge}{which edges have been treated}
}
\references{Mann, M. E., & Lees, J. M. (1996). Robust estimation of background noise and signal detection in climatic time series. Climatic change, 33(3), 409-445.}

\seealso{
code{\link{gausssmooth}}
}
\examples{
X<-rnorm(50,mean=0,sd=1)+seq(0,10,length.out=50)
plot(X,type="l")
lines(seq(0,10,length.out=50))
lines(medsmoo(X,10,5)$rmean,col="red")
points(medsmoo(X,10,5,"right")$rmean,col="orange")
points(medsmoo(X,10,5,"left")$rmean,col="blue2")
}

