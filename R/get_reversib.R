#' Get the reversibility of time series
#' 
#' Get the horizontal visibility graph of time series, and from that, test the
#' reversibility
#' 
#' 
#' @usage get_reversib(xin, plt = TRUE)
#' @param xin vector or \code{zoo}-object
#' @param plt logical; set to \code{TRUE} if plot is desired
#' @return A list of \item{res}{List with sub-items for \code{clustering} and
#' \code{degree}. Large p-values indicate reversible, small p-values
#' irreversible dynamics} \item{A}{Adjacency matrix} \item{tvec}{Sampling
#' times} \item{xin}{Input data}
#' @author Kira Rehfeld
#' @seealso \code{\link{degclust}}
#' @references TBD
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' 
#' # nonlinear process
#' x<-zoo(rnorm(100)^2,order.by=generate_t())
#' out<-get_reversib(x,plt=)
#' out[["res"]][["clustering"]][["pval"]]
#' out[["res"]][["degree"]][["pval"]]
#' 
#' # linear/reversible process
#' y<-zoo(rnorm(100),order.by=generate_t())
#' out2<-get_reversib(y,plt=TRUE)
#' out2[["res"]][["clustering"]][["pval"]]
#' out2[["res"]][["degree"]][["pval"]]
#' 
#' @export get_reversib
get_reversib<-function(xin,plt=TRUE){
A<-D<-matrix(0,length(xin),length(xin))
x<-coredata(xin)

for (i in 1:(length(x)-1)){
for (j in (i+1):length(x)){
cmp<-x[(i+1):j]
if (length(cmp)==1) {D[i,j]<-1; next}
cmp<-cmp[-length(cmp)]
if (all(cmp<=x[i])&&all(cmp<=x[j])) {D[i,j]<-1}
}
}

A<-(D+t(D))

t.out<-index(xin)


if (plt){
plot(t.out,coredata(xin))

for (i in 1:(length(x)-1)){
for (j in (i+1):length(x)){

if (A[i,j]==1) {lines(c(t.out[i],t.out[j]),c(x[i],x[j]),col=i);
}

}
}
}

res<-degclust(A)
out<-ks.test(res$degree$back,res$degree$forw,exact=FALSE)
res$degree$pval<- out[["p.value"]]


out<-ks.test(res$clustering$back,res$clustering$forw,exact=FALSE)
res$clustering$pval<- out[["p.value"]]


return(list(A=A,tvec=t.out,xin=xin,res=res))

}


