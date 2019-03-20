#' Estimate the AR1 persistence time
#' 
#' Estimate the AR1 persistence time
#' 
#' 
#' @usage tauest(X,method="acf",deltas=mean(diff(index(X))))
#' @param X \code{zoo} object
#' @param method "acf" or "ls"
#' @param deltas reference time step
#' @return An estimate of the persistence time
#' @author Kira Rehfeld
#' @references Mann, M. E., & Lees, J. M. (1996). Robust estimation of
#' background noise and signal detection in climatic time series. Climatic
#' change, 33(3), 409-445. Rehfeld, K., Marwan, N., Heitzig, J., & Kurths, J.
#' (2011). Comparison of correlation analysis techniques for irregularly
#' sampled time series. Nonlinear Processes in Geophysics.
#' @examples
#' 
#' dtx=1
#' X<-generate_ar1(list(generate_t()),phi=0.5)
#' tau=tauest(X,deltas=1,method="ls")
#' phi.est=exp(-1/tau);
#' 
#' @export tauest
tauest <- function(X,method="acf",deltas=mean(diff(index(X)))){
  
  y<-coredata(X)
  yt<-index(X)
  
  t<-NA
  
  meths<-c("acf","ls")
  methind<-match(method,meths,nomatch=3)
  
  if (length(y)<=10) {warning("time series too short"); return(t)}

  r1=nexcf(X,lag=deltas)
  if (is.na(r1)){return(t)}
  

  
  if (methind==3){
  warning("Check spelling of method for tauest, reverting to ACF method")
  methind=1
  }
  

  if (methind==1){#acf
  
  if (r1<0) {warning("Persistence smaller than sampling time, set to mean time step"); t<-deltas} else { 

  t<--deltas/log(r1) }
  
  }

  if (methind==2){# least squares
  
  if (r1>0.01) {
  taustart=-deltas/log(r1)
  yi<-y[-1]
  yL<-y[-length(y)]
  delt<-diff(yt)
  #mod<-nls(yi ~ int+exp(-delt/tau)*yL, start=list(tau=taustart,int=0),lower=c(tau=0,int=0),upper=c(tau=max(yt)),algorithm="port")
  mod<-nls(yi ~ int+exp(-delt/tau)*yL, start=list(tau=taustart,int=0), control=nls.control(warnOnly = TRUE))
  t<-coef(mod)[1]
  
  if (t<=0){t<-NA; warning("tauest: tau smaller than 0. This could be due to numerical issues.")} # if tau is determined (by whatever numerical mistake) as below zero, set to a numerically acceptable number close to zero.
  if (t>(max(yt)-min(yt))) warning("tauest: tau is larger than the time series is long!")
  
  #if (!(mod$convergence==0)){t=taustart} 
  } else {
  #if autocorrelation is practically zero
  if (r1>0){ #autocorrelation could be there, but minute -- set it from lag-1 acf as the numerical routines will fail 
  t=max(min(diff(yt)),-deltas/log(abs(r1)))
  warning("tauest: estimated tau from autocorrelation function due to numerical issues in the nonlinear regression routine")
  } else { #lag-1 autocorrelation is already negative. Set the persistence time to a numerically small value (0 would cause divide-by-zero errors).
  t=NA}
  }
  }
 
  
  return(t)
}
