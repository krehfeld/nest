\name{degclust}
\alias{degclust}

\title{Compute degree and clustering for small complex networks}
\description{Compute degree and clustering for small complex networks. For time-ordered nodes, the network measures are computed irrespectively of the time (total), in forward, and backward sense.}
\usage{
degclust(A)
}

\arguments{
  \item{A}{Adjacency matrix. diag(A) should be zero (no self-loops).}
}

\value{
%%  ~Describe the value returned
A list of
\item{degree }{A list of \code{forward},\code{backward}, and \code{tot}al degree}
\item{clustering }{A list of \code{forward},\code{backward}, and \code{tot}al local clustering}
}
\references{
Watts, D.J. & Strogatz S. (1998) Collective dynamics of small-world networks, Nature 393 (1998), doi:10.1038/30918
Rehfeld, K. , Molkenthin, N. and Kurths, J. (2014) Testing the detectability of spatio-temporal climate transitions from paleoclimate networks with the START model, Nonlinear Processes in Geophysics, 21 (3), pp. 691-703. doi:10.5194/npg-21-691-2014}
\author{Kira Rehfeld}


\examples{
A<-matrix(1,nrow=4,ncol=4); diag(A)=0;
degclust(A)$degree$tot
degclust(A)$clustering$tot
A[3,4]=A[4,3]=0
degclust(A)$clustering$tot
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
