#' AR1 surrogates for irregularly sampled time series.
#' 
#' This function generates time series with a lag-1 autocorrelation mimicking
#' that of the supplied data.
#' 
#' 
#' @usage ar1sur(X, N = 3)
#' @param X \code{zoo-object}, time series for which surrogates are to be
#' generated.
#' @param N Integer number of surrogates desired
#' @return List of surrogate time series
#' @note Time series with linear trends should be detrended prior to modeling.
#' The estimation of the autocorrelation strength will otherwise likely fail.
#' @author Kira Rehfeld
#' @seealso \code{\link{car}}
#' @references Rehfeld, K., Marwan, N., Heitzig, J. and Kurths, J. (2011)
#' Comparison of correlation analysis techniques for irregularly sampled time
#' series, Nonlinear Processes in Geophysics, 18 (3), pp. 389-404.
#' doi:10.5194/npg-18-389-2011
#' @examples
#' 
#' ar1sur(zoo(rnorm(100)))
#' 
#' @export ar1sur
ar1sur <-
function(X,N=3){
#Generate AR1 data, with temporal sampling as in tx and an autocorrelation coefficient of phi. 
#browser()
tx<-index(X)
x<-coredata(X)
taux<-tauest(X)
phi<-exp(-mean(diff(tx))/taux)
sdx<-sd(X)

if (phi>1) stop("phi shouldn't be above 1")


#obs. time matrix
TX<-matrix(tx,nrow=N,ncol=length(tx),byrow=TRUE) # Nxlength(tx) matrix

#get ar1 surrogates
Xsur<-generate_ar1(TX,phi) # list with output series

#Xout<-t(do.call(cbind,Xsur))
# match the variance of the input data
normalize<-function(xs,sdout=1,mout=0){xs<-sdout*(xs-mean(xs))/sd(xs)+mout}
Xsur<-lapply(Xsur,normalize,sdout=sdx,mout=mean(x))

return(Xsur) # returns a list of AR1 surrogates for the passed time series 

}
