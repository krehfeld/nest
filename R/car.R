#' Coupled AR-1 processes
#' 
#' Generate coupled AR-1 processes on arbitrary time axes
#' 
#' 
#' @usage car(tx, ty, coupl_strength, phi, lag, nsur = 1)
#' @param tx time axis for time series 1 (vector)
#' @param ty time axis for time series 2
#' @param coupl_strength Double value between \code{-1} and \code{1}
#' @param phi Autocorrelation of the time series 1
#' @param lag Time lag for the coupling
#' @param nsur Number of realizations to be generated
#' @return List of time series
#' @author Kira Rehfeld
#' @seealso \code{\link{generate_t}},\code{\link{ar1sur}},
#' @references Rehfeld, K., Marwan, N., Heitzig, J. and Kurths, J. (2011)
#' Comparison of correlation analysis techniques for irregularly sampled time
#' series, Nonlinear Processes in Geophysics, 18 (3), pp. 389-404.
#' doi:10.5194/npg-18-389-2011
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' tx<-generate_t(dt=1,tmin=0,tmax=100,method="gamma")
#' ty<-generate_t(dt=1,tmin=0,tmax=100,method="gamma")
#' Proc<-car(tx,ty,coupl_strength=0.7,phi=0.5,lag=0,nsur=1)
#' x<-zoo(Proc$x,order.by=tx)
#' y<-zoo(Proc$y,order.by=ty)
#' phi.est=nexcf(x,lag=1,h=0.25)
#' couplstr.est=nexcf(x=x,y=y,lag=0,h=0.05)
#' 
#' @export car
car <-
function(tx,ty,coupl_strength,phi,lag,nsur=1){
#Coupled autoregressive processes 

# set the resolution at which the AR1process is simulated at high
# resolution as 10 times the actual mean resolution

dtsur=min(mean(diff(tx)),mean(diff(ty)))/10;
t1=min(min(tx), min(ty))-lag; # subtract lag  to compensate transient effects;
t2=max(max(tx),max(ty))+lag; #add lag to avoid extrapolation at the end
#make interpolation timescale

t=seq(t1,t2+2*dtsur,dtsur);


#persistence time rescaled for higher time series resolution
fi=exp(-dtsur/(-1/log(phi)));


####<----------------------------------------------------------- beginning of loop

# First step: Create driving process at high resolution
eps1=rnorm(length(t),0,1);
xsur=eps1;
for(i in 2:length(t)){
xsur[i]=fi*xsur[i-1] +sqrt(1-fi^2)*eps1[i];
# do something
}

xsur=(xsur-mean(xsur))/sd(xsur);

# Second step: Couple the second time series
eps2=rnorm(length(t),0,1);
ysur=eps2;

L=ceiling(lag/dtsur)

end=length(t);
ysur[(L+1):end]=coupl_strength*xsur[1:(length(t)-L)]+sqrt(1-coupl_strength^2)*eps2[(L+1):length(t)];


ysur=(ysur-mean(ysur))/sd(ysur);


# interpolate linearly to the original timescale desired
intoutx=approx(t,xsur,tx,method="linear",rule=1)
intouty=approx(t,ysur,ty,method="linear",rule=1)


# assign output
x<-intoutx$y;
y<-intouty$y;	

### end of loop


return(list(x=x,y=y))

}
