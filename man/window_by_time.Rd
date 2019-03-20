\name{window_by_time}
\alias{window_by_time}
\title{Cut an irregular zoo time series into windows using the value of the index}
\description{
This function subsets a time series using the value of the index, not the index number. It respects the irregular sampling of the zoo-timeseries, and returns a list of sub-time series.
}
\usage{
window_by_time(X, T0 = start(X)[1], T1 = end(X)[1], wwidth = (T1 - T0)/10, shift = 0.5)
}

\arguments{\item{X}{zoo time series object, or list of zoo time series.}
  \item{T0}{Start cutting value of the time series. Set, for example, to \code{T0}=0 and \code{T1}=1000 if multiple series are to be cut.}
  \item{T1}{Last value of the time series for cutting.}
  \item{wwidth}{Width of the sub-time-series window (in time units)}
  \item{shift}{Relative shift of the time series window. Should be above 0 and below/equal to 1. A value of 0.5 results in windows that are 50 \% overlapping. A value of 0.9 in 90\% overlaps.}
}

\value{Returns a "window-by-time-object" containing a list of:
\item{tsplit}{List of sub-time series}
\item{tmid}{Mid-points of the windows}
}

\author{Kira Rehfeld}

\seealso{
\code{\link{network_stats}},\code{\link{network_links}},\code{\link{quality_check}}
}

\examples{
## Example for window_by_time (univariate windowing)

## time series with increasing variance
X.full<-zoo(x=rnorm(1000)*seq(1,1000))
## Make time series gappy/irregular by omitting a percentage of points
X<-sample(X.full,250)
## To see the impact this subsampling has on the time series, plot the difference of the time index 
plot(X.full)
lines(X,col="red")
plot(diff(index(X)))
hist(diff(index(X)))

## Obtain windowed time series - by time, not by index
X.split<-window_by_time(X,shift=0.5,wwidth=15)
X.win<-X.split$tsplit


# Get time series of the variance and the median
X.var<-zoo(x=sapply(X.win,var),order.by=X.split$tmid)
X.median<-zoo(x=sapply(X.win,median),order.by=X.split$tmid)

plot(X.var,main="Checking stationarity of the variance")
plot(X.median,main="Checking stationarity of the mean")
abline(h=0,lty=2)

######################### Multivariate use

# Generate a list of long time series time series, one for each network node

Xlist<-lapply(replicate(3,zoo(x=rnorm(1000)*seq(1,1000)),simplify=FALSE),sample,250)

# window each time series

OUT<-window_by_time(Xlist,T0=0,T1=1000,shift=0.75)

# strip out only time series
OUT.tls<-lapply(OUT,function(x){x$tsplit}) #list of time series only
# get the number of points in each window
OUT.length<-sapply(OUT.tls,function(sublist){sapply(sublist,length)}) 

# get the standard deviations in each time window at each node
OUT.var<-sapply(OUT.tls,function(sublist){sapply(sublist,sd)}) 

# get the time vector (the mid point of each time window)
OUT.tme<-OUT[[1]]$tmid

# plot variance against time for each node
matplot(OUT.tme,OUT.var,type="l")
matpoints(OUT.tme,OUT.var)

# Another simple example: Get the covariances between all time series (no windowing)
f1 <- function(a,b, fun){
  outer(a, b, function(x,y) vapply(seq_along(x), function(i) fun(x[[i]], y[[i]]), numeric(1)))
}
f1(Xlist,Xlist,nexcf)
f1(Xlist,Xlist,function(x,y){nexcf(x,y,lag=1)}) #lag1 cross-correlation for all

}

