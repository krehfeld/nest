#' Compute degree and clustering for small complex networks
#' 
#' Compute degree and clustering for small complex networks. For time-ordered
#' nodes, the network measures are computed irrespectively of the time (total),
#' in forward, and backward sense.
#' 
#' 
#' @usage degclust(A)
#' @param A Adjacency matrix. diag(A) should be zero (no self-loops).
#' @return %% ~Describe the value returned A list of \item{degree }{A list of
#' \code{forward},\code{backward}, and \code{tot}al degree} \item{clustering
#' }{A list of \code{forward},\code{backward}, and \code{tot}al local
#' clustering}
#' @author Kira Rehfeld
#' @references Watts, D.J. & Strogatz S. (1998) Collective dynamics of
#' small-world networks, Nature 393 (1998), doi:10.1038/30918 Rehfeld, K. ,
#' Molkenthin, N. and Kurths, J. (2014) Testing the detectability of
#' spatio-temporal climate transitions from paleoclimate networks with the
#' START model, Nonlinear Processes in Geophysics, 21 (3), pp. 691-703.
#' doi:10.5194/npg-21-691-2014
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' A<-matrix(1,nrow=4,ncol=4); diag(A)=0;
#' degclust(A)$degree$tot
#' degclust(A)$clustering$tot
#' A[3,4]=A[4,3]=0
#' degclust(A)$clustering$tot
#' 
#' @export degclust
degclust<-function(A){
is.full.adj<-all(A==t(A))


degree.tot<-apply(A,2,sum)

clust.tot<-clust.back<-clust.forw<-degree.forw<-degree.back<-matrix(NA,ncol=dim(A)[2])
for (k in 1:dim(A)[2]){
ind.forw<-ind.back<-ind.conn<-NA

if (degree.tot[k]<1) {

clust.forw[k]<-clust.back[k]<-clust.tot[k]<-0

} else {

#if (degree.tot[k]==0) clust[k]=0
ind.conn<-which(A[k,]>0)

if (any(ind.conn<k)){ # backward} 
ind.back<-ind.conn[ind.conn<k]
degree.back[k]<-length(ind.back)
}
if (any(ind.conn>k)){ # backward} 
ind.forw<-ind.conn[ind.conn>k]
degree.forw[k]<-length(ind.forw)

}


## now split into forward and backward
if (length(ind.conn)>1){

clust.tot[k]<-(sum(A[ind.conn,ind.conn]))/((degree.tot[k]*(degree.tot[k]-1)))
# this works for the total degree, but for the forward/backward degree I still have problems
# nach rechts
if (length(ind.forw)>1){
clust.forw[k]<-(sum(A[ind.forw,ind.forw]))/((degree.forw[k]*(degree.forw[k]-1)))
#clust.back[k]<-(sum(A[which(A[,k]>0),which(A[,k]>0)])/((degree.forw[k]*(degree.forw[k]-1))/2))/2
} else {clust.forw[k]=0}

if (length(ind.back)>1){
clust.back[k]<-(sum(A[ind.back,ind.back]))/((degree.back[k]*(degree.back[k]-1)))
#clust.back[k]<-(sum(A[which(A[,k]>0),which(A[,k]>0)])/((degree.forw[k]*(degree.forw[k]-1))/2))/2
} else {clust.back[k]=0}


} else {clust.tot[k]<-clust.back[k]<-clust.forw[k]<-0}
}


}
res.clust<-list(tot=clust.tot,forw=clust.forw,back=clust.back)
res.degree<-list(tot=degree.tot,forw=degree.forw,back=degree.back)

return(list(degree=res.degree,clustering=res.clust))
}
