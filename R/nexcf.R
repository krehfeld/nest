#' Pearson correlation for irregular time series
#' 
#' Compute the correlation for irregular time series with or without
#' interpolation.
#' 
#' 
#' @aliases nexcf ixcf nexcf_ci
#' @usage nexcf(x, y = x, lag = 0, h = 0.25, enforce = FALSE) ixcf(x, y = x,
#' lag = 0) nexcf_ci(x, y, conflevel = 0.05, dt = NULL, smoo =
#' FALSE,detr=FALSE)
#' @param x time series one to be analyzed. Note that time increments should be
#' positive (i.e. time is reported as age before present).
#' @param y time series two to be analyzed (if omitted, autocorrelation is
#' computed)
#' @param lag numeric vector giving lags for which correlation is to be
#' estimated
#' @param h number giving the Gaussian kernel width
#' @param enforce Coerce the result to the interval (-1,1) in case of numerical
#' problems
#' @param conflevel Confidence-level for the two-sided confidence intervals
#' around zero correlation
#' @param dt Average time step for the correlation function evaluation (by
#' default set to max(mean(dtx),mean(dty)
#' @param smoo logical; set to TRUE for smoothing of variability below the time
#' step
#' @param detr logical; set to TRUE for linear detrending of the time series
#' @return Numerical vector of the length of the lag vector.
#' @author Kira Rehfeld
#' @seealso \code{\link{generate_t}}
#' @references Rehfeld, K. and Kurths, J.: Similarity estimators for irregular
#' and age-uncertain time series, Clim. Past, 10, 107-122,
#' doi:10.5194/cp-10-107-2014, 2014. Rehfeld, K., Marwan, N., Heitzig, J., and
#' Kurths, J.: Comparison of correlation analysis techniques for irregularly
#' sampled time series, Nonlin. Processes Geophys., 18, 389-404,
#' doi:10.5194/npg-18-389-2011, 2011.
#' @keywords correlation irregular time series
#' @examples
#' 
#' ## Generate one gamma-distributed and one regular time axis
#' tx<-generate_t(dt=1,tmin=0,tmax=100,method="gamma")
#' ty<-generate_t(dt=1,tmin=0,tmax=100,method="linear")
#' ## Simulate one coupled AR1 process (see reference for details)
#' Proc<-car(tx,ty,coupl_strength=0.7,phi=0.5,lag=0,nsur=1)
#' ## Bind the results to zoo time series
#' x<-zoo(Proc$x,order.by=tx)
#' y<-zoo(Proc$y,order.by=ty)
#' ## Estimate the autocorrelation
#' phi.est=nexcf(x,lag=1,h=0.25)
#' ## Estimate the cross-correlation
#' couplstr.est=nexcf(x=x,y=y,lag=0,h=0.05)
#' ## Estimate the cross-correlation with interpolation
#' couplstr.est.interp=ixcf(x=x,y=y,lag=0)
#' 
#' ## Estimate the effective degrees of freedom, and the significance of the correlation estimate,
#' ## from the t-distribution.
#' temp<-nexcf_ci(x,y,conflevel=0.1)
#' # degrees of freedom used for the significance test
#' print(temp$neff)
#' print(temp$rxy)
#' # if rxy is outside 
#' temp$ci
#' # the confidence intervals, the correlation estimate is significant at the (1-conflevel)*100% level
#' 
#' @export nexcf
nexcf <- function(x,y=x,lag=0,h=.25,enforce=FALSE) {
# Gaussian kernel correlation

gaussfun <-function(dt,h=0.25){	result<-1/(sqrt(2*pi)*h) * exp(-dt^2/(2*h^2))}

# initialize output
corr<-matrix(NA,length(lag))

### prechecks

if (class(x)!=class(y)) warning("nexcf: x and y don't have the same class")

if ((length(x)==0)||length(y)==0) { warning("nexcf: one or both time series have length zero, returning NA."); return(corr) }
if (all(diff(index(x))<0)|all(diff(index(y))<0)) {
if (all(diff(index(x))<0)&&all(diff(index(y))<0)) {
x<-x[length(x):1];y<-y[length(y):1]
warning("nexcf: time increments are negative for both time series. This renders the windowing ineffective, therefore the time series have been flipped around. The most likely cause is that your time series is in AD/CE, rather than BP")
} else {
stop("Time increments positive for one but not the other time series. Check time definitions")}

}

T0<-max(index(x)[1],index(y)[1]) #index(x)[1].... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T1<-min(index(x)[length(x)],index(y)[length(y)])
# Using this cutting definition assumes that your time series is ordered, with positive time increments, otherwise the windowing will be ineffective







# check overlap in time axes, if no overlap -> error
if ((T1-T0)<=0) { warning("nexcf: time series don't overlap, returning NA."); return(corr) }

## cut to overlap exactly!
x<-window(x,start=T0,end=T1)
y<-window(y,start=T0,end=T1)


if ((length(x)<10)||length(y)<10) { warning("nexcf: less than 10 points in overlap, returning NA."); return(corr) }

if ((length(x)>10000)||length(y)>10000) { warning("nexcf: more than 10000 obs in one time series - this will take a lot of memory and time to compute.")}

tx<-index(x);
ty<-index(y);



if (length(lag)==1){
dtlag<-max(mean(diff(index(x))),mean(diff(index(y))))
} else {
dtlag<-mean(diff(lag))
}

if ((T1-T0)/dtlag<=10) { warning("nexcf: time series don't overlap sufficiently, returning NA."); return(corr) }

## Check if autocorrelation should be computed

if ((length(tx)==length(ty))&&all(tx==ty)&&(x==y)) {is.acf=TRUE} else {is.acf=FALSE}
# find lag zero

### preprocessing: 
# time index dimensionless
tx<-tx/dtlag; x<-zoo(x,order.by=tx)
ty<-ty/dtlag; y<-zoo(y,order.by=ty)
normlag<-lag/dtlag

# normalize (assumption!)
x<-scale(x)
y<-scale(y)#(y-mean(y,na.omit=TRUE))/sd(y)

### Main part
### TBD: optimize for long time series (see nexcf_largedata.R in spielwiese)
# data product matrix xydist
tdist_xy<--outer(index(x),index(y),"-")# pairwise distances of a vector
pdist_xy<-(x)%*%t(y)

#tdist_xx<--outer(index(x),index(x),"-")# pairwise distances of a vector
#pdist_xx<-(x-mean(x))%*%t(x-mean(x))


#tdist_yy<--outer(index(y),index(y),"-")# pairwise distances of a vector
#pdist_yy<-(y-mean(y))%*%t(y-mean(y))

for (i in 1:length(normlag)){
# for each lag build new distance matrix corresponding to the considered lag
t_dist_d<-tdist_xy+normlag[i]

# matrix of weights (-> apply kernel)
WL_xy<-gaussfun(t_dist_d,h)
#WL_xx<-gaussfun(tdist_xx+normlag[i],h)
#WL_yy<-gaussfun(tdist_yy+normlag[i],h)

#Weight<-sum(WL_xy,na.rm=TRUE)
# check kernel content
kernel_cont<-(t_dist_d<=5*h)&(t_dist_d >=(-5*h))
if (sum(WL_xy[kernel_cont],na.rm=TRUE)<0.5){
warning("kernel content insufficient - maybe use larger kernel")
corr[i]<-NA
next
}


Cov_xy<-sum(pdist_xy*WL_xy)/sum(WL_xy,na.rm=TRUE)
#Cov_xx<-sum(pdist_xx*WL_xx)/sum(WL_xx,na.rm=TRUE)
#Cov_yy<-sum(pdist_yy*WL_yy)/sum(WL_yy,na.rm=TRUE)


#C_<-Cov_xy/(sqrt(Cov_xx)*sqrt(Cov_yy))
# weighted sum of products
#COV_xy<-(sum(pdist_xy*WL_xy)/Weight)

#VARxVary<-sd(x)*sd(y)

#if (VARxVary<COV_xy) {VARxVary=COV_xy}
C_<-Cov_xy#/VARxVary

corr[i]<-C_
}

#if (is.acf&&length(lag)>1){corr<-corr/corr[which(lag==0)]}
#if (enforce&&(max(abs(corr))>1)){corr<-corr/corr[which(lag==0)]}
mx<-max(abs(corr),na.rm=TRUE)

if (mx>1){ # need enforcement or warning
if ((enforce)&&is.acf&&(any(lag==0))){corr<-corr/corr[which(lag==0)]}
if ((enforce)&&!is.acf){corr<-corr/(max(abs(corr)))}

}


#if (any(abs(corr)>1)){warning("Numerical difficulties in nexcf: Abs(Corr)>1",call.=TRUE)}

return(corr)
}
#TBD: return list of values, and the timescale of the lag -> depends on the lag vector
