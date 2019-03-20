\name{gaussfun}
\alias{gaussfun}

\title{Simple Gaussian kernel function}
\description{Simple Gaussian kernel function}
\usage{gaussfun(dt, h = 0.25)}

\arguments{
  \item{dt}{time difference}
  \item{h}{kernel width}
}

\value{Numeric}


\examples{
#Compare the weights given in nexcf to observations at a time lag of 0.1 and 1.
nest:::gaussfun(0.1,h=0.25)
nest:::gaussfun(1,h=0.25)
# If the sampling of two time series is very skewed, it may pay off to increase the kernel width h.
# h is set to 0.25 (lag units) by default, and increasing it moderately, e.g. to 0.5 is acceptable.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
