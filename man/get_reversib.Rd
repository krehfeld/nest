\name{get_reversib}
\alias{get_reversib}

\title{Get the reversibility of time series}
\description{Get the horizontal visibility graph of time series, and from that, test the reversibility}
\usage{get_reversib(xin, plt = TRUE)}

\arguments{
  \item{xin}{vector or \code{zoo}-object}
  \item{plt}{logical; set to \code{TRUE} if plot is desired}
}

\value{A list of
\item{res}{List with sub-items for \code{clustering} and \code{degree}. Large p-values indicate reversible, small p-values irreversible dynamics}
\item{A}{Adjacency matrix}
\item{tvec}{Sampling times}
\item{xin}{Input data}
}
\references{TBD}
\author{Kira Rehfeld}

\seealso{\code{\link{degclust}}}
\examples{

# nonlinear process
x<-zoo(rnorm(100)^2,order.by=generate_t())
out<-get_reversib(x,plt=)
out[["res"]][["clustering"]][["pval"]]
out[["res"]][["degree"]][["pval"]]

# linear/reversible process
y<-zoo(rnorm(100),order.by=generate_t())
out2<-get_reversib(y,plt=TRUE)
out2[["res"]][["clustering"]][["pval"]]
out2[["res"]][["degree"]][["pval"]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
