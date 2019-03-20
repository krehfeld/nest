nexcf_ci<-function(x,y,conflevel=0.05,dt=NULL,smoo=FALSE,detr=FALSE){
# NEFF_CI computes ccf effective sample size and confidence intervals for the Gaussian-Kernel-based cross correlation
# Formula adapted from N_eff derived from Dawdy and Matalas, 1964.
# Dawdy, D.R., and Matalas, N.C., 1964, Statistical and probability analysis of hydrologic data,part III: Analysis of variance, covariance and time series, in Ven Te Chow, ed., Handbook of applied hydrology, a compendium of water-resources technology: New York, McGraw-Hill Book Company, p. 8.68-8.90.
# Von Storch, Zwiers 1999
detrend_linear<-function(zoobj){
if (!is.zoo(zoobj)){return(zoobj)}
x<-index(zoobj);
y<-coredata(zoobj);
res<-lm(y~x); 
yprime<-y-res$coefficients[2]*x+res$coefficients[1];
zoobdetr<-zoo(yprime,order.by=x);
return(zoobdetr)
}
## NEW reference (17062015) via Legendre and Legendre, numerical ecology: 
#An important points is that in correlation or regression analysis, spatial correlation has a deleterious effect on tests of significance only when it is present in both variables. Simulation studies have shown that when spatial correlation was present in only one of the two variables, the test had a correct rate of type I error (Bivand/Legendre2002). These simulations have also shown that deterministic spatial structures present in both vatriables have the same effect as spatial autocorrelation. For example, with a deterministic structure in one of the variables and spatial autocorrelation in the other, tests of significance had inflated rates of type I error.
#Bivand, Roger. "A Monte Carlo study of correlation coefficient estimation with spatially autocorrelated observations." Quaestiones Geographicae 6 (1980): 5-10.
#Legendre, Pierre, et al. "The consequences of spatial structure for the design and analysis of ecological field surveys." Ecography 25.5 (2002): 601-615.
#

# conflevel -> symmetric confidence levels around zero
# p-value -> percentile of the correlation value on the t-distribution of an correlation estimate around zero (null hypothesis) p=0.08 -> neg. correl. p=0.98 -> pos. correl.
# neff -> conservative estimate of the effective degrees of freedom, corrected for (positive) autocorrelation.
# rxy -> correlation computed using nexcf

reslist<-list(rxy=NA,ci=NA,neff=NA,tau=c(NA,NA),pval=NA,ci.rxy=NA)

if (length(x)<10||length(y)<10) {warning("one or both time series are too short length");return(reslist)}
tx<-index(x);ty<-index(y)
x<-coredata(x);y<-coredata(y)

# timescale adapted to sampling of the "coarser" time series
dtx=mean(diff(tx))
dty=mean(diff(ty))
if (is.null(dt)) dt=max(dtx,dty)

if (dt<max(dtx,dty)){
warning("Chosen timescale is below the sampling rate of the time series")
} 

if (smoo==TRUE){
x<-gausssmooth(tx=tx,x=x,h=dt/6)
y<-gausssmooth(tx=y,x=y,h=dt/6)
}

if (detr){
x<-coredata(detrend_linear(zoo(x,order.by=tx)))
y<-coredata(detrend_linear(zoo(y,order.by=ty)))
}



# Input length equivalent to length of the "shorter" time series
nx=length(tx);ny=length(ty);
Rx<-range(tx)[2]-range(tx)[1]
Ry<-range(ty)[2]-range(ty)[1]

taux<-tauest(zoo(x,order.by=tx),deltas=dt)
tauy<-tauest(zoo(y,order.by=ty),deltas=dt)

## get effective number of observations
#neff<-max(nx/taux,ny/tauy); #this would be the standard formula given by von Storch/Zwiers, applicable if taux/tauy were in units of samples, not time.
neff<-min(max(Rx/taux,Ry/tauy,na.rm=TRUE),max(nx,ny))


## Cross-Correlation
rxy=nexcf(zoo(x,order.by=tx),zoo(y,order.by=ty),lag=c(-dt,0,dt),h=0.25,enforce=TRUE)[2]

if (is.na(rxy)){return(reslist)}

## Get the test value for the significance test
testval=qnorm(1-conflevel/2);
CI1=testval/sqrt(neff-2);

## Confidence interval around zero
ci=c(-CI1,CI1);

## Confidence interval around the correlation coefficient using Fisher's transform
fishTransf<-function(rxy){0.5*(log(1+rxy)-log(1-rxy))}
sigmaz<-1/sqrt(neff-3)
ldo<-fishTransf(rxy)-(testval*sigmaz)
lup<-fishTransf(rxy)+(testval*sigmaz)
ci.rxy<-c(NA,NA)
ci.rxy[1]<-(exp(2*ldo)-1)/(exp(2*ldo)+1)
ci.rxy[2]<-(exp(2*lup)-1)/(exp(2*lup)+1)

##################################################################################
# compare this to r.con(rxy,neff,p=pval) from the psych package
# get p-value from the correlation estimate by using the fact that the r2 value is t-distributed with neff degrees of freedom. From the t-distribution we can then get the p-value corresponding to the correlation estimate. This does not explicitly use the fact that the time series may be irregularly sampled. If the correlation estimate was obtained via interpolation, the value should be compared against estimates from surrogates.
lowtail=FALSE # take the p-value that a correlation GREATER than the one estimated is found from normally distributed uncorrelated data
if (rxy<0){lowtail=TRUE}
rtot<-rxy * sqrt((neff - 2)/(1 - rxy^2))
#pval<-pt(r2t(rxy,neff),neff-2,lower.tail=lowtail)
pval<-pt(rtot,neff-2,lower.tail=lowtail)

reslist<-list(rxy=rxy,ci=ci,neff=neff,tau=c(taux,tauy),pval=pval,ci.rxy=ci.rxy)

# return everything 
return(reslist)

}

