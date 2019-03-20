#' Median Smoother
#' 
#' Mean and median smoothing of numerical vectors with optional edge
#' preservation.
#' 
#' 
#' @usage medsmoo(x, m, m2, edge = c("both"))
#' @param x numeric vector
#' @param m half width of the smoother
#' @param m2 half width of the smoother at the edges
#' @param edge "left", "right" or "both" for application of smaller windows
#' towards the edges
#' @return
#' 
#' A list of \item{rmedian }{running median} \item{rmean }{running mean}
#' \item{rmad }{running median absolute deviation}
#' \item{smootherwidth}{smoother width for each point} \item{edge}{which edges
#' have been treated}
#' @seealso code\link{gausssmooth}
#' @references Mann, M. E., & Lees, J. M. (1996). Robust estimation of
#' background noise and signal detection in climatic time series. Climatic
#' change, 33(3), 409-445.
#' @examples
#' 
#' X<-rnorm(50,mean=0,sd=1)+seq(0,10,length.out=50)
#' plot(X,type="l")
#' lines(seq(0,10,length.out=50))
#' lines(medsmoo(X,10,5)$rmean,col="red")
#' points(medsmoo(X,10,5,"right")$rmean,col="orange")
#' points(medsmoo(X,10,5,"left")$rmean,col="blue2")
#' 
#' @export medsmoo
medsmoo <-
function(x,m,m2,edge=c("both")){

if (m2>m|m2<0) { stop("medsmoo: m2 should be smaller than m and larger than 0")
}

e=pmatch(edge,c("left","right","both"))


N<-length(x)
rmedian<-rmad<-rmean<-vector()

M<-rep(m,length(x)) # set up width vector

for (i in 1:(m2+2)) { # up to the m2+1th data point, use the width of m2
if (e==1|e==3){M[i]=m2}
if (e==3||e==2){M[N-i+1]=m2}
}

for (i in (m2+2):m){ # up to the m'th data point, and from the m2th, use increasingly wider 
if (e==1|e==3){M[i]=i-1}
if (e==3||e==2){M[N-i+1]=i-1}
}

for (i in 1:N) # for each point

	{
		i0<-max(i-M[i],1)#from the first point
		i1<-min(i+M[i],N)#to maximum the last point
		part<-x[i0:i1]#select these points between the indices
		rmedian[i]<-median(part)#get their median/mad/mean
		rmad[i]<-mad(part)
		rmean[i]<-mean(part)
	}


result<-list(rmedian=rmedian,rmad=rmad,rmean=rmean,smootherwidth=M,edge=edge)
class(result)<-"msmooth"
return(result)	

}
