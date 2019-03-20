\name{quality_check}
\alias{quality_check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quality checks and convenient splitting of time series}
\description{
Quality checks and convenient splitting of time series
}
\usage{
quality_check(prxlist, meta, T0, T1, maxhiat = 800, length.min = 15, 
min.res = 500, min.range=(T1-T0)/3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prxlist}{List of zoo objects}
  \item{meta}{metadata}
  \item{T0}{Vector of start points}
  \item{T1}{Vector of end points}
  \item{maxhiat}{Maximum hiatus length in time series}
  \item{length.min}{Minimum acceptable length of the time series}
  \item{min.res}{Minimum acceptable resolution}
  \item{min.range}{Minimum time covered (in time units)}
  
}

\value{
Returns a list of
\item{WBT.out}{Window-by-time-object}
\item{metafilt}{Filtered metadata object}
\item{ind.ok}{Indices of acceptable objects}
}

\note{TBD: Floating window definition (e.g. Eemian by thermal max. in window)}


\seealso{
\code{\link{window_by_time}}
}
\examples{

# Generate a list of long time series time series, one for each network node
Xlist<-lapply(replicate(3,zoo(x=rnorm(1000)*seq(1,1000)),simplify=FALSE),sample,250)
metafake<-data.frame(c(-75,30,75),c(120,120,120),c("a","b","c"))
colnames(metafake)<-c("Lat","Lon","Name")

# window each time series
#OUT<-lapply(Xlist,window_by_time,T0=0,T1=1000,shift=0.75)

OUT<-window_by_time(Xlist,T0=c(0,500),T1=c(500,1000))

OUT<-quality_check(Xlist,metafake,T0=0,T1=1000,maxhiat=10,length.min=100,min.res=10)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
