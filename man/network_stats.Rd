\name{network_stats}
\alias{network_stats}
\alias{network_links}
\alias{networkmap_simple}
\title{Paleoclimate network statistics}

\description{Functions to compute paleoclimate network node and link statistics, and to visualize them on a simple map.}
\usage{
network_stats(WBTobj,fcn="length",metanames=1:length(WBTobj),optargs=NULL)
network_links(WBTobj,conflevel=0.1,detr=FALSE,dt.lag=NULL)
networkmap_simple(CMAT,lat,lon,weights=rep(1,dim(CMAT)[1]),thresh= NULL)
}

\arguments{
  \item{WBTobj}{"window-by-time"-object created from window_by_time}
  \item{fcn}{Function for node statistics. One of \code{c("length","var","mean","median","tau","quantile","dof","
","own")}}
  \item{metanames}{Optional: Names for the time series objects}
  \item{optargs}{List of optional additional arguments for network_stats, supplied to the chosen function. For example \code{optargs=list(probs=0.9))}, which is a double value between 0 and 1 used when \code{fcn}=="quantile", and \code{method}, which gives the method to compute the effective time series length/ degrees of freedom for \code{fcn} equal to "dof". For \code{fcn}=="var.tsc" (timescale-dependent variance) - timescales in time units of the supplied time series. See \code{\link{tsc_dep_var}} for more information. For \code{fcn}=="own", a function named "myfunc" can be supplied with additional input arguments. This is then applied to all sub-time-windows. See example below.}

For \code{network_links}:
\item{conflevel}{Confidence level for the correlation estimates.}
\item{detr}{Logical; detrending option for \code{nexcf_ci}}
\item{dt.lag}{Lag at which \code{nexcf} is to be evaluated}

For \code{networkmap_simple}:
\item{lat}{Latitude (in \code{c(-90,+90)})}
\item{lon}{Longitude (in \code{c(-180,+180)})}
\item{weights}{Weights for the node sizes. Symbols will be scaled by the square-root of the weight vector.}
\item{CMAT}{Correlation matrix}
\item{thresh}{Numeric; Correlation threshold}
}

\value{
For \item{network_stats}{an object of the class "palnet.nodes"}
For \item{network_links}{an object of the class "palnet.links"}
}
\references{Rehfeld, K. , Molkenthin, N. and Kurths, J. (2014) Testing the detectability of spatio-temporal climate transitions from paleoclimate networks with the START model, Nonlinear Processes in Geophysics, 21 (3), pp. 691-703. doi:10.5194/npg-21-691-2014
Rehfeld, K., Marwan, N., Breitenbach, S. F. M. and Kurths, J. (2013) Late Holocene Asian summer monsoon dynamics from small but complex networks of paleoclimate data,Climate Dynamics, 41 (1), pp. 3-19. doi:10.1007/s00382-012-1448-3}
\author{Kira Rehfeld krehfeld@awi.de}


\seealso{\code{\link{nexcf}}, \code{\link{window_by_time}},\code{\link{effective_ts_length}}}
\examples{
# for networks
library(maps)
library(mapdata)
#library(PaleoSpec)
#library(zoo)
# Generate a list of long time series time series, one for each network node
Xlist<-lapply(replicate(3,zoo(x=rnorm(1000)*seq(1,1000)),simplify=FALSE),sample,500)
metafake<-data.frame(c(-75,30,75),c(120,120,120),c("a","b","c"))
colnames(metafake)<-c("Lat","Lon","Name")

# window each time series
#OUT<-lapply(Xlist,window_by_time,T0=0,T1=1000,shift=0.75)

OUT<-window_by_time(Xlist,T0=0,T1=1000,shift=0.75)



# Several functions can be computed for each time window and each node
statlist<-c("length","var","mean","median","tau","quantile","dof")

for (i in 1:(length(statlist))){
 #ST<-network_stats(OUT,fcn=statlist[i],optargs=list(probs=0.9,tsc.in=c(1,10)))
ST<-network_stats(OUT,fcn=statlist[i])
 
 matplot(ST$tvec,ST[["MAT"]],type="l",main=statlist[i])
}

# "quantile" and "own" are to be used with additional input arguments:
ST<-network_stats(OUT,fcn="quantile",optargs=list(probs=0.5))
matplot(ST$tvec,ST[["MAT"]],type="l")
#matplot(ST$tvec,network_stats(OUT,fcn="quantile",optargs=list(probs=0.9))$MAT,type="l")

tsc.in=c(10,40)
ST<-network_stats(OUT,fcn="var.tsc",optargs=list(pval=0.1,tsc.in=c(10,40)))
matplot(ST$tvec,ST[["MAT"]],type="p")
matlines(ST$tvec,ST[["misc"]][["var.ci"]]$lo,type="l")
matlines(ST$tvec,ST[["misc"]][["var.ci"]]$up,type="l")

matplot(ST$tvec,ST[["misc"]][["var.ci"]][["up"]],type="n")
x0<-ST$tvec;
lower.ci<-ST[["misc"]]$var.ci$lo;upper.ci<-ST[["misc"]]$var.ci$up;
arrows(x0=x0,y0=lower.ci[,1],y1=upper.ci[,1],angle=90,col=1,code=3,length=0.05)
arrows(x0=x0,y0=lower.ci[,2],y1=upper.ci[,2],angle=90,col=2,code=3,length=0.05)
arrows(x0=x0,y0=lower.ci[,3],y1=upper.ci[,3],angle=90,col=3,code=3,length=0.05)
matpoints(ST$tvec,ST$MAT,pch=c(16:18),lwd=2)




# Example with user-supplied function (median absolute deviation of the time series) in each window
#library(psych)
#myfunc<-function(x){skew(x)}
myfunc<-function(x){median(abs(x-median(x)))}
ST<-network_stats(OUT,fcn="own",optargs=list(myfunc=myfunc))

# myfunc<-function(x){effective_ts_length(x)}
# ST<-network_stats(OUT,fcn="own",optargs=list(myfunc=myfunc))

matplot(ST$tvec,ST[["MAT"]],type="l")


##################### GET THE CORRELATIONS BETWEEN ALL NETWORK TIME SERIES
## network_links: Get the cross-correlation matrix for each time window

#OUT<-lapply(Xlist,window_by_time,T0=0,T1=1000,shift=0.75)
Node.sd<- sqrt(network_stats(OUT,fcn="var")$MAT)

NW<-network_links(OUT,conflevel=0.95)

# Adjacency matrix for the network
A<-1*(NW[["PMAT"]]<0.1)
image(A[1,,])

# node weights
par(mfcol=c(4,4),mar=c(1,1,1,1))
for (i in 1:length(NW$tmid)){
timept<-i
CMAT<-NW$CMAT[timept,,]
networkmap_simple(CMAT=CMAT,lat=metafake$Lat,lon=metafake$Lon,weights=Node.sd[timept,])
rm(CMAT)
}

# make some fake connection matrix with many nodes
nnodes=10
CMAT<-matrix(runif(nnodes^2,min=-1,max=1),nnodes,nnodes)
weights<-runif(nnodes,1,3)/nnodes
networkmap_simple(CMAT=CMAT,lat=runif(nnodes,-85,85),lon=runif(nnodes,-180,180),weights=weights)
# automatic node size setting is to be fine-tuned

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
