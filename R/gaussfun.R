#' Simple Gaussian kernel function
#' 
#' Simple Gaussian kernel function
#' 
#' 
#' @usage gaussfun(dt, h = 0.25)
#' @param dt time difference
#' @param h kernel width
#' @return Numeric
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' #Compare the weights given in nexcf to observations at a time lag of 0.1 and 1.
#' nest:::gaussfun(0.1,h=0.25)
#' nest:::gaussfun(1,h=0.25)
#' # If the sampling of two time series is very skewed, it may pay off to increase the kernel width h.
#' # h is set to 0.25 (lag units) by default, and increasing it moderately, e.g. to 0.5 is acceptable.
#' 
gaussfun <-
function(dt,h=0.25){
	result<-1/(sqrt(2*pi)*h) * exp(-dt^2/(2*h^2))}
