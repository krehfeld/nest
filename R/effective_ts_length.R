#' degrees of freedom of an irregular autocorrelated time series.
#' 
#' Get the effective length, or degrees of freedom, of an irregular
#' autocorrelated time series, using different methods.
#' 
#' 
#' @usage effective_ts_length(X, method = "tau")
#' @param X A zoo-time-series object.
#' @param method One of "tau", "mudelsee", or "simple"
#' @return Approximated number of degrees of freedom, estimated using the
#' autocorrelation of the time series.
#' @author Kira Rehfeld, with contribution from Thom Laepple.
#' @seealso \code{\link{nexcf_ci}}
#' @references von Storch, H., and Zwiers, F.W., 1999. Statistical Analysis in
#' Climate Research. Cambridge University Press [Page number to be checked]
#' Dawdy, D.R., and Matalas, N.C., 1964, Statistical and probability analysis
#' of hydrologic data,part III: Analysis of variance, covariance and time
#' series, in Ven Te Chow, ed., Handbook of applied hydrology, a compendium of
#' water-resources technology: New York, McGraw-Hill Book Company, p.
#' 8.68-8.90. Mudelsee,M., 2010. Climate Time Series Analysis - Classical
#' Statistical and Bootstrap Methods. Springer ISBN: 9789048194827, 474 p [Page
#' number for reference to be checked]
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' 
#' ## Generate one gamma-distributed and one regular time axis
#' tx<-generate_t(dt=1,tmin=0,tmax=100,method="gamma")
#' ty<-generate_t(dt=1,tmin=0,tmax=100,method="linear")
#' ## Simulate one coupled AR1 process (see reference for details)
#' Proc<-car(tx,ty,coupl_strength=0.5,phi=0.5,lag=0,nsur=1)
#' ## Bind the results to zoo time series
#' x<-zoo(Proc$x,order.by=tx)
#' y<-zoo(Proc$y,order.by=ty)
#' ## Estimate the autocorrelation
#' phi.est=nexcf(x,lag=1,h=0.25)
#' 
#' plot(x)
#' lines(y,col="red")
#' 
#' rxy=nexcf(x,y,lag=0,h=0.25)
#' rx1=nexcf(x,lag=1)
#' ry1=nexcf(y,lag=1)
#' 
#' nexcf_ci(x,y)
#' 
#' 
#' effective_ts_length(x,method="tau") # adapted from von Storch/Zwiers 1999
#' effective_ts_length(x,method="mudelsee") # after Mudelsee, 2010
#' effective_ts_length(x,method="simple") # using both time directions
#' 
#' @export effective_ts_length
effective_ts_length<-function(X,method="tau"){
# Get effective length of a time series, corrected for autocorrelation.
# Useful for significance tests

poss.methods<-c("tau","mudelsee","simple")

methind<-match(method,poss.methods,nomatch=1)
nx<-length(X)

if (methind==1){
# Von Storch/ Zwiers (check page number)
# Formula adapted from N_eff derived from Dawdy and Matalas, 1964.
# Dawdy, D.R., and Matalas, N.C., 1964, Statistical and probability analysis of hydrologic data,part III: Analysis of variance, covariance and time series, in Ven Te Chow, ed., Handbook of applied hydrology, a compendium of water-resources technology: New York, McGraw-Hill Book Company, p. 8.68-8.90.
# Von Storch, Zwiers 1999
Rx<-range(index(X))[2]-range(index(X))[1]
taux<-tauest(X,method="acf")
neff<-min(Rx/taux,nx)
}

if (methind==2){
# Mudelsee 2010 (check page number)
r1<-nexcf(X,lag=mean(diff(index(X))),enforce=TRUE)

eff.mudelsee<-function(r1,n) n/(1+2/n*1/(1-r1)*(r1*(n-1/(1-r1))-r1^n*(1-1/(1-r1))))

neff<-min(eff.mudelsee(r1,nx),nx)

}

if (methind==3){
r1<-nexcf(X,lag=mean(diff(index(X))),enforce=TRUE)

eff.simple<-function(r1,n) n*(1-r1)/(1+r1)
neff<-min(eff.simple(r1,nx),nx)


}


return(neff)

}
