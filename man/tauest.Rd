\name{tauest}
\alias{tauest}
\title{
Estimate the AR1 persistence time
}
\description{Estimate the AR1 persistence time}
\usage{
tauest(X,method="acf",deltas=mean(diff(index(X))))
}

\arguments{
  \item{X}{\code{zoo} object}
  \item{method}{"acf" or "ls"}
  \item{deltas}{reference time step}
}

\value{
An estimate of the persistence time
}
\references{
Mann, M. E., & Lees, J. M. (1996). Robust estimation of background noise and signal detection in climatic time series. Climatic change, 33(3), 409-445.
Rehfeld, K., Marwan, N., Heitzig, J., & Kurths, J. (2011). Comparison of correlation analysis techniques for irregularly sampled time series. Nonlinear Processes in Geophysics.
}
\author{Kira Rehfeld}

\examples{
dtx=1
X<-generate_ar1(list(generate_t()),phi=0.5)
tau=tauest(X,deltas=1,method="ls")
phi.est=exp(-1/tau);
}
