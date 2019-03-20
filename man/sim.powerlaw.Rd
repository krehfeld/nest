\name{sim.powerlaw}
\alias{sim.powerlaw}

\title{Simulate a power law time series}
\description{Simulate a power-law time series}
\usage{
sim.powerlaw(beta, N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{power-law exponent}
  \item{N}{length of the time series}
}
\note{The Filter used in the power law generation is currently unscaled. To force Variance(X)=1, a scaling is performed on the final timeseries. This is problematic for short time series and if the distribution of variances is of interest. A workaround is to generate a very long powerlaw series (with overall variance=1) and cut it into many shorter ones. x<-matrix(sim.powerlaw(betas[b],Nsur*Nlength),ncol=Nsur,nrow=Nlength)}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
