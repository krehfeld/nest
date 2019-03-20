\name{ar1sur}
\alias{ar1sur}

\title{AR1 surrogates for irregularly sampled time series.}
\description{This function generates time series with a lag-1 autocorrelation mimicking that of the supplied data.}
\usage{
ar1sur(X, N = 3)
}
\arguments{
  \item{X}{\code{zoo-object}, time series for which surrogates are to be generated.}
  \item{N}{Integer number of surrogates desired}
}

\value{List of surrogate time series}

\references{Rehfeld, K., Marwan, N., Heitzig, J. and Kurths, J. (2011)
Comparison of correlation analysis techniques for irregularly sampled time series,
Nonlinear Processes in Geophysics, 18 (3), pp. 389-404. doi:10.5194/npg-18-389-2011}
\author{Kira Rehfeld}
\note{Time series with linear trends should be detrended prior to modeling. The estimation of the autocorrelation strength will otherwise likely fail.}

\seealso{
\code{\link{car}}
}
\examples{
ar1sur(zoo(rnorm(100)))
}

