#' Paleoclimate network statistics
#' 
#' Functions to compute paleoclimate network node and link statistics, and to
#' visualize them on a simple map.
#' 
#' 
#' @aliases network_stats network_links networkmap_simple
#' @usage
#' network_stats(WBTobj,fcn="length",metanames=1:length(WBTobj),optargs=NULL)
#' network_links(WBTobj,conflevel=0.1,detr=FALSE,dt.lag=NULL)
#' networkmap_simple(CMAT,lat,lon,weights=rep(1,dim(CMAT)[1]),thresh= NULL)
#' @param WBTobj "window-by-time"-object created from window_by_time
#' @param fcn Function for node statistics. One of
#' \code{c("length","var","mean","median","tau","quantile","dof"," ","own")}
#' @param metanames Optional: Names for the time series objects
#' @param optargs List of optional additional arguments for network_stats,
#' supplied to the chosen function. For example
#' \code{optargs=list(probs=0.9))}, which is a double value between 0 and 1
#' used when \code{fcn}=="quantile", and \code{method}, which gives the method
#' to compute the effective time series length/ degrees of freedom for
#' \code{fcn} equal to "dof". For \code{fcn}=="var.tsc" (timescale-dependent
#' variance) - timescales in time units of the supplied time series. See
#' \code{\link{tsc_dep_var}} for more information. For \code{fcn}=="own", a
#' function named "myfunc" can be supplied with additional input arguments.
#' This is then applied to all sub-time-windows. See example below.
#' @param conflevel Confidence level for the correlation estimates.
#' @param detr Logical; detrending option for \code{nexcf_ci}
#' @param dt.lag Lag at which \code{nexcf} is to be evaluated
#' @param lat Latitude (in \code{c(-90,+90)})
#' @param lon Longitude (in \code{c(-180,+180)})
#' @param weights Weights for the node sizes. Symbols will be scaled by the
#' square-root of the weight vector.
#' @param CMAT Correlation matrix
#' @param thresh Numeric; Correlation threshold
#' @return For \item{network_stats}{an object of the class "palnet.nodes"} For
#' \item{network_links}{an object of the class "palnet.links"}
#' @author Kira Rehfeld krehfeld@@awi.de
#' @seealso \code{\link{nexcf}},
#' \code{\link{window_by_time}},\code{\link{effective_ts_length}}
#' @references Rehfeld, K. , Molkenthin, N. and Kurths, J. (2014) Testing the
#' detectability of spatio-temporal climate transitions from paleoclimate
#' networks with the START model, Nonlinear Processes in Geophysics, 21 (3),
#' pp. 691-703. doi:10.5194/npg-21-691-2014 Rehfeld, K., Marwan, N.,
#' Breitenbach, S. F. M. and Kurths, J. (2013) Late Holocene Asian summer
#' monsoon dynamics from small but complex networks of paleoclimate
#' data,Climate Dynamics, 41 (1), pp. 3-19. doi:10.1007/s00382-012-1448-3
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' # for networks
#' library(maps)
#' library(mapdata)
#' #library(PaleoSpec)
#' #library(zoo)
#' # Generate a list of long time series time series, one for each network node
#' Xlist<-lapply(replicate(3,zoo(x=rnorm(1000)*seq(1,1000)),simplify=FALSE),sample,500)
#' metafake<-data.frame(c(-75,30,75),c(120,120,120),c("a","b","c"))
#' colnames(metafake)<-c("Lat","Lon","Name")
#' 
#' # window each time series
#' #OUT<-lapply(Xlist,window_by_time,T0=0,T1=1000,shift=0.75)
#' 
#' OUT<-window_by_time(Xlist,T0=0,T1=1000,shift=0.75)
#' 
#' 
#' 
#' # Several functions can be computed for each time window and each node
#' statlist<-c("length","var","mean","median","tau","quantile","dof")
#' 
#' for (i in 1:(length(statlist))){
#'  #ST<-network_stats(OUT,fcn=statlist[i],optargs=list(probs=0.9,tsc.in=c(1,10)))
#' ST<-network_stats(OUT,fcn=statlist[i])
#'  
#'  matplot(ST$tvec,ST[["MAT"]],type="l",main=statlist[i])
#' }
#' 
#' # "quantile" and "own" are to be used with additional input arguments:
#' ST<-network_stats(OUT,fcn="quantile",optargs=list(probs=0.5))
#' matplot(ST$tvec,ST[["MAT"]],type="l")
#' #matplot(ST$tvec,network_stats(OUT,fcn="quantile",optargs=list(probs=0.9))$MAT,type="l")
#' 
#' tsc.in=c(10,40)
#' ST<-network_stats(OUT,fcn="var.tsc",optargs=list(pval=0.1,tsc.in=c(10,40)))
#' matplot(ST$tvec,ST[["MAT"]],type="p")
#' matlines(ST$tvec,ST[["misc"]][["var.ci"]]$lo,type="l")
#' matlines(ST$tvec,ST[["misc"]][["var.ci"]]$up,type="l")
#' 
#' matplot(ST$tvec,ST[["misc"]][["var.ci"]][["up"]],type="n")
#' x0<-ST$tvec;
#' lower.ci<-ST[["misc"]]$var.ci$lo;upper.ci<-ST[["misc"]]$var.ci$up;
#' arrows(x0=x0,y0=lower.ci[,1],y1=upper.ci[,1],angle=90,col=1,code=3,length=0.05)
#' arrows(x0=x0,y0=lower.ci[,2],y1=upper.ci[,2],angle=90,col=2,code=3,length=0.05)
#' arrows(x0=x0,y0=lower.ci[,3],y1=upper.ci[,3],angle=90,col=3,code=3,length=0.05)
#' matpoints(ST$tvec,ST$MAT,pch=c(16:18),lwd=2)
#' 
#' 
#' 
#' 
#' # Example with user-supplied function (median absolute deviation of the time series) in each window
#' #library(psych)
#' #myfunc<-function(x){skew(x)}
#' myfunc<-function(x){median(abs(x-median(x)))}
#' ST<-network_stats(OUT,fcn="own",optargs=list(myfunc=myfunc))
#' 
#' # myfunc<-function(x){effective_ts_length(x)}
#' # ST<-network_stats(OUT,fcn="own",optargs=list(myfunc=myfunc))
#' 
#' matplot(ST$tvec,ST[["MAT"]],type="l")
#' 
#' 
#' ##################### GET THE CORRELATIONS BETWEEN ALL NETWORK TIME SERIES
#' ## network_links: Get the cross-correlation matrix for each time window
#' 
#' #OUT<-lapply(Xlist,window_by_time,T0=0,T1=1000,shift=0.75)
#' Node.sd<- sqrt(network_stats(OUT,fcn="var")$MAT)
#' 
#' NW<-network_links(OUT,conflevel=0.95)
#' 
#' # Adjacency matrix for the network
#' A<-1*(NW[["PMAT"]]<0.1)
#' image(A[1,,])
#' 
#' # node weights
#' par(mfcol=c(4,4),mar=c(1,1,1,1))
#' for (i in 1:length(NW$tmid)){
#' timept<-i
#' CMAT<-NW$CMAT[timept,,]
#' networkmap_simple(CMAT=CMAT,lat=metafake$Lat,lon=metafake$Lon,weights=Node.sd[timept,])
#' rm(CMAT)
#' }
#' 
#' # make some fake connection matrix with many nodes
#' nnodes=10
#' CMAT<-matrix(runif(nnodes^2,min=-1,max=1),nnodes,nnodes)
#' weights<-runif(nnodes,1,3)/nnodes
#' networkmap_simple(CMAT=CMAT,lat=runif(nnodes,-85,85),lon=runif(nnodes,-180,180),weights=weights)
#' # automatic node size setting is to be fine-tuned
#' 
#' 
#' @export network_stats
network_stats<-function(WBTobj,fcn="length",metanames=1:length(WBTobj),optargs=NULL){
statlist<-c("length","var","mean","median","tau","quantile","dof","var.tsc","own")
ind<-match(fcn,statlist)
# optargs is a list of optional arguments and replaces the ...
# WBTobj must be a list of lists of ts created with window_by_time

if (!class(WBTobj)=="window-by-time-object") stop("Check arguments. Function needs a window-by-time-object created with window_by_time()")

# The outer list of WBTobj has two properties, tsplit and tmid (the time points) - these have to be separated
stripped_obj<-lapply(WBTobj,function(x){x$tsplit}) #list of time series only
tvec<-WBTobj[[1]]$tmid #list of time series only

OUT.stat=matrix(NA,nrow=length(tvec),ncol=length(WBTobj))
rownames(OUT.stat)<-tvec
colnames(OUT.stat)<-metanames

misc=list()
if (ind==1) {
#length

mylength<-function(x){if(is.null(x)){Q<-NA} else {Q<-length(x)}; return(Q)}
OUT.stat<-sapply(stripped_obj,function(sublist){sapply(sublist,mylength)}) # get matrix of subseries lengths
}

if (ind==2) {
# varian
myvar<-function(x){if(is.null(x)){Q<-NA} else {Q<-var(x)}; return(Q)}

OUT.stat<-sapply(stripped_obj,function(sublist){sapply(sublist,myvar)}) # get matrix of subseries lengths
}
if (ind==3) {
# mean
mymean<-function(x){if(is.null(x)){Q<-NA} else {Q<-mean(x)}; return(Q)}
mysd<-function(x){if(is.null(x)){Q<-NA} else {Q<-sd(x)}; return(Q)}
mylength<-function(x){if(is.null(x)){Q<-NA} else {Q<-length(x)}; return(Q)}

OUT.stat<-sapply(stripped_obj,function(sublist){sapply(sublist,mymean)}) # get matrix of subseries lengths
OUT.sd<-sapply(stripped_obj,function(sublist){sapply(sublist,mysd)}) # get matrix of subseries lengths
OUT.l<-sapply(stripped_obj,function(sublist){sapply(sublist,mylength)}) # get matrix of subseries lengths
OUT.se<-OUT.sd/sqrt(OUT.l)

misc$SE<-OUT.se
misc$sd<-OUT.sd
misc$l<-OUT.l
}

if (ind==4) {
# median
mymedian<-function(x){if(is.null(x)){Q<-NA} else {Q<-median(x)}; return(Q)}

OUT.stat<-sapply(stripped_obj,function(sublist){sapply(sublist,mymedian)}) # get matrix of subseries lengths
}

if (ind==5) {
# tau
OUT.stat<-sapply(stripped_obj,function(sublist){sapply(sublist,function(x){tauest(x)})}) # get matrix of subseries lengths
}

if (ind==6){
# quantile
## what if optargs is NULL? -> set probs to 0.9
#argnames <- names(optargs)
#print(argnames)
#print(args$probs)
#args<-list(...)
        # check whether probs is an argument
        #if(!("probs" %in% argnames)){
        if (length(optargs)==0){
		optargs<-list(probs=0.9) 
                #probs <- 0.9 #90% quantile
        #        warning("setting 90% quantile as default. If this is not what you want, supply the probs argument")
        } 
myquantfunc<-function(x,optargs){if (length(x)<3){Q=rep(NA,length(optargs$probs))} else {optargs$x<-x;Q<-do.call(quantile,args=optargs)}; return(Q)}
OUT.stat<-sapply(stripped_obj,function(sublist,optargs){sapply(sublist,function(x,optargs){myquantfunc(x,optargs)},optargs)},optargs) # get matrix of subseries lengths
}

if (ind==7){
#print("dof not yet estimated")
#argnames <- names(optargs)
 #  if(!("method" %in% argnames)) {
 if (length(optargs)==0){
 optargs<-list(method="tau") #90% quantile
                #warning("setting 90% quantile as default. If this is not what you want, supply the probs argument")
        } 
# method can be "tau", "mudelsee", "simple" 
#OUT.stat<-sapply(stripped_obj,function(sublist,method){sapply(sublist,function(x,method){effective_ts_length(x,method=method)},method=method)},method=method) # get matrix of subseries lengths

OUT.stat<-sapply(stripped_obj,function(sublist,optargs){sapply(sublist,function(x,optargs){if (is.null(x)){Q<-NA} else {optargs$X<-x;Q<-do.call(effective_ts_length,args=optargs)}; return(Q)},optargs)},optargs) # get matrix of subseries lengths


# DOFs!!
# dof.effective.mudelsee -> Mudelsee 2010
# dof.effective.simple # where does this formula come from?
#neff<-range/taux
#100/(-1/log(0.7)) # why is this a factor of 2 of the other two estimators? -> because this is the 1/e point?
#dof.effective.simple(0.7,100)
#dof.effective.simple<-function(a1,n) n*(1-a1)/(1+a1)
#mudelsee 2010?
#dof.effective.mudelsee<-function(a1,n) n/(1+2/n*1/(1-a1)*(a1*(n-1/(1-a1))-a1^n*(1-1/(1-a1))))
# von storch/zwiers 1999
#decorrelation time tau0=-dt/log(a1)
#dof.effective.storch=n*tau/tau0=-n*log(a1) # same as mine!
#100/(-2/(log(0.7))) # -> full decorrelation instead of half decorrelation
# why is there the difference? matze's idea: past & future autocorrelation? or full decorrelation length to 1/2e
#dof.effective.mudelsee(0.7,100)
#-100*log(0.7)
}

if ((ind==8)&(match("var.tsc",fcn))){
# get tsc dependent variance
# three additional input arguments - timescale, conf level and estimator
# three outputs
tsc.in=pval=NULL
print(optargs)

argnames <- names(optargs)
#attach(optargs)
#print(optargs)

#print(tsc.in)
# print(methvar)
   if(!("pval" %in% argnames)) {
                pval <- 0.1 #90% quantile
        } 
#    if(!("methvar" %in% argnames)) {#spec/sd/gaussbapa+sd
#                 methvar <- "spec" #90% quantile
#         }
   if(!("tsc.in" %in% argnames)||!(length(optargs$tsc.in)==2)) {#spec/sd/gaussbapa+sd
                #tsc.in <- "spec" #90% quantile
                stop("Supply timescale vector [e.g. tsc.in=c(1,10)] for variance")
        }

misc$var.ci<-list(lo=c(),up=c())


# misc$var.ci$lo<-sapply(stripped_obj,sapply,function(x,pval,methvar,tsc.in){tsc_dep_var(x,pval=pval,methvar=methvar,tsc.in=tsc.in)$var.ci$lo},pval=pval,methvar=methvar,tsc.in=tsc.in)

misc$var.ci$lo<-sapply(stripped_obj,sapply,function(x,optargs){if (is.null(x)){return(NA)} else {optargs$timser<-x;print(x);
Q<-do.call(tsc_dep_var,args=optargs);return(Q$var.ci$lo)}},optargs)



# da sind abwechselnd die lo/up CIs vermixt -> müssen irgendwie auseinander klamüsert werden
# misc$var.ci$up<-sapply(stripped_obj,sapply,function(x,pval,methvar,tsc.in){tsc_dep_var(x,pval=pval,methvar=methvar,tsc.in=tsc.in)$var.ci$up},pval=pval,methvar=methvar,tsc.in=tsc.in)
misc$var.ci$up<-sapply(stripped_obj,sapply,function(x,optargs){if (is.null(x)){return(NA)} else { optargs$timser<-x;
Q<-do.call(tsc_dep_var,args=optargs);return(Q$var.ci$up)}},optargs)

# misc$dof<-sapply(stripped_obj,function(sublist,...){sapply(sublist,function(x,...){tsc_dep_var(x,...)$dof},...)},...) # get matrix 
misc$dof<-sapply(stripped_obj,sapply,function(x,optargs){if (is.null(x)){return(NA)} else { optargs$timser<-x;
Q<-do.call(tsc_dep_var,args=optargs);return(Q$dof)}},optargs)

# misc$ts.used<-lapply(stripped_obj,function(sublist,...){lapply(sublist,function(x,...){tsc_dep_var(x,...)$ts.used},...)},...) # get matrix 
misc$ts.used<-lapply(stripped_obj,lapply,function(x,optargs){if (is.null(x)){return(NA)} else { optargs$timser<-x;
Q<-do.call(tsc_dep_var,args=optargs);return(Q$ts.used)}},optargs)


# dof, ci
#$misc$var.ci<-var.ci
# # OUT.stat<-sapply(stripped_obj,sapply,function(x,pval,methvar,tsc.in){tsc_dep_var(x,pval=pval,methvar=methvar,tsc.in=tsc.in)$var.out},pval=pval,methvar=methvar,tsc.in=tsc.in)
OUT.stat<-sapply(stripped_obj,sapply,function(x,optargs){if (is.null(x)){return(NA)} else { optargs$timser<-x;
Q<-do.call(tsc_dep_var,args=optargs);return(Q$var.tsc)}},optargs)

#OUT.stat<-sapply(stripped_obj,function(sublist,...){sapply(sublist,function(x,...){tsc_dep_var(x,...)$var.out},...)},...) # get matrix of subseries lengths       
}

if (ind==9){# own function
argnames <- names(optargs)
myfunc=NULL
        if(!("myfunc" %in% argnames)) {
                stop("myfunc function is not supplied, or not named myfunc")
        } 
attach(optargs)
# supply your own function in the ... argument as myfunc
OUT.stat<-sapply(stripped_obj,function(sublist){sapply(sublist,function(x){myfunc(x)})}) # get matrix of subseries lengths

}

results<-list(MAT=OUT.stat,tvec=tvec,misc=misc)
class(results)<-"palnet.nodes"
return(results)
}
