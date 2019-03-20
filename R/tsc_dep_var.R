#' Timescale-dependent variance estimation
#' 
#' Estimate the timescale-dependent variance of a time series.
#' 
#' 
#' @usage tsc_dep_var(timser, tsc.in, min.res = min(mean(diff(index(timser))),
#' floor(tsc.in[1]/2) - 1), start.val = start(timser), end.val = end(timser),
#' pval = 0.1, detrend = FALSE)
#' @param timser Input time series. Can be a \code{zoo}, \code{ts} or vector
#' object.
#' @param tsc.in Vector of two timescales, \code{c(tsc1,tsc2)}.
#' @param min.res Minimal resolution of the time series (if this is not met,
#' NAs are returned).
#' @param start.val Starting value of the time series window.
#' @param end.val End value of the time series window.
#' @param pval \code{pval} gives the p-value for the confidence intervals of
#' the variance estimate.
#' @param detrend logical argument, if set to \code{TRUE} the timesereis is
#' linearly detrended prior to spectral analysis
#' @return A list of \item{std}{Standard deviation}
#' \item{var.tsc}{timescale-dependent variance estimated} \item{var.tot}{total
#' variance of the time series in the window} \item{dof}{estimated degrees of
#' freedom of the variance estimate} \item{ts.used}{used time series}
#' \item{var.ci}{list of \code{lo} and \code{up}, lower and upper confidence
#' level for the variance estimate} \item{spec.win}{spectrum object}
#' @author Kira Rehfeld, with contributions from Thom Laepple
#' @seealso \code{\link{var_ci}}
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' library(scales)
#' ## Generate one gamma-distributed and one regular time axis
#' tx<-generate_t(dt=1,tmin=0,tmax=250,method="gamma")
#' ty<-generate_t(dt=1,tmin=0,tmax=250,method="linear")
#' ## Simulate one coupled AR1 process (see reference for details)
#' Proc<-car(tx,ty,coupl_strength=0.5,phi=0.5,lag=0,nsur=1)
#' ## Bind the results to zoo time series
#' x<-zoo(Proc$x,order.by=tx)
#' y<-zoo(Proc$y,order.by=ty)
#' 
#' 
#' tsc_dep_var(y,tsc.in=c(30,100)) 
#' 
#' @export tsc_dep_var
tsc_dep_var<-function(timser,tsc.in,min.res=min(mean(diff(index(timser))),floor(tsc.in[1]/2)-1),start.val=start(timser),end.val=end(timser),pval=0.1,detrend=FALSE){


### Definitions: 
results<-list(std=NA,var.tsc=NA,var.tot<-NA,dof=NA,ts.used=timser,var.ci=list(up=NA,lo=NA))
std<-var<-dof<-var.tsc<-var.tot<-ts.used<-c(NA);

### Checks:
# print(any(is.na(timser)))

if (!(is.ts(timser)|is.zoo(timser))){warning("Vector, not time series (zoo or ts) supplied - assuming regular sampling at unit intervals")}
if (length(timser)<15){warning("Time series length is lower than 15 pts, returning NA"); return(results)}
if (mean(diff(index(timser)))>min.res){warning("Time series mean sampling step is larger than minimum resolution, returning NA"); return(results)}

# Checking the order of the supplied timescales
if (tsc.in[2]<tsc.in[1]) tsc.in<-tsc.in[2:1]
fcut=1/tsc.in; #omega.upper omega.lower

if (2*min.res>tsc.in[1]){ 
 print(min.res)
 print(tsc.in[1])
 print(tsc.in[2])
stop("Desired bandpass attempts to go beyond nyquist frequency - check min.res and tsc.in")
}



# Other possible cutoff for the : freq.min =1/4*max.dt similar to nyquist frequency 

# select a minimum timestep for the high-resolution interpolation
dt.hres<-0.99*min(min.res/10,min(diff(index(timser))))

### Computation


##? does it have to be min.res, or does it have to be constant? 
#if (!all(diff(index(timser))==min.res)){
if ((!is.regular(timser,strict=TRUE))|(any(is.na(timser)))){ # if the time series is not regular, make it regular
#tmptmp<-nest:::interp.antialias(timser,dt.out=min.res,dt.hres=NULL,k=5,kf=1.2,
#startTime=start.val)
#first<-function(x){x[1]};last<-function(x){x[length(x)]};
tmptmp<-na.omit(MakeEquidistant(index(timser),coredata(timser),dt=min.res,dt.hres=
NULL , k = 5 , kf = 1.2, time.target=seq(from=start.val,to=index(timser)[length(timser)],by=min.res)))

# # var.in<-var(window(timser,start=start.val,end=end.val),na.rm=TRUE)
# # var.out<-var(window(tmptmp,start=start.val,end=end.val),na.rm=TRUE)
# # print(paste("Overall Variance change from",var.in,"to" , var.out, "due to interpolation i.e. by",round((var.out/var.in-1)*100,3),"percent"))
#print(tmptmp)
} else {tmptmp<-timser} # if the time series is already regular (and at a resolution better than min.res) then continue with it

# cut to start and end times
tmp<-window(tmptmp,start=start.val,end=end.val)


spec.win<-spec.pgram(tmp,pad=0,taper=FALSE,plot=FALSE,detrend=detrend,fast=FALSE)
#contains p [power] f [frequency] df [degrees of freedom] and 

# find the indices of the spectral estimates in the tsc.in-window 
ind<-c(which.min(abs(fcut[1]-spec.win$f)):which.min(abs(fcut[2]-spec.win$f)))

if (length(ind)< 3) warning("Less than three spectral estimates in the frequency window")

# take the mean over these spectral estimates and multiply by 2* bandwidth to get a correct variance
var.tsc<-mean(spec.win$spec[ind])*(max(spec.win$f[ind])-min(spec.win$f[ind]))*2
var.tot<-var(tmp)

# get the degrees of freedom of this variance estimate
dof<-length(ind)*spec.win$df

std<-sqrt(var.tsc)
ts.temp<-tmp




if (pval>0.1) warning("p-value larger than .1 - is this on purpose or misspecified significance level?")
var.ci<-var_ci(var.tsc,dof,pval=pval)

results<-list(std=std,var.tsc=var.tsc,var.tot=var.tot,dof=dof,ts.used=ts.temp,var.ci=var.ci,spec.win=spec.win)
return(results)
}

