ixcf <-
function(x,y=x,lag=0) {
# interpolated correlation
# initialize output
#corr<-matrix(NA,length(lag))

### prechecks
tx<-index(x)
ty<-index(y)

if (class(x)!=class(y)) stop("ixcf: x and y don't have the same class")
# check length of time series x and y, also wrt time vector

#browser()

#check ts class? (is.ts(x)&is.
#coerce from ts class?
T0<-min(tx[1],ty[1]) #index(x)[1].... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T1<-max(tx[length(x)],ty[length(y)])

# check overlap in time axes, if no overlap -> error
if ((T1-T0)<=0) stop("ixcf: time series don't overlap")

if (length(lag)==1){
dtlag=max(mean(diff(tx),mean(diff(ty))))
} else {dtlag=mean(diff(lag))
}
normlag<-lag/dtlag

tx<-tx/dtlag
ty<-ty/dtlag


### preprocessing: 
### interpolation
if (isTRUE(all.equal(tx,ty))&&(dtlag=1)){#((length(x)!=length(y))&&(index(x)!=index(y))&&(dtlag!=1)) { # if anything has to be interpolated
ix<-x
iy<-y
} else { # perform interpolation #((length(x)!=length(y))&&(index(x)!=index(y))&&(dtlag!=1))
ti<-seq(from=T0,to=T1,by=dtlag)/dtlag# common time axis interpolated to ...|
#ti<-seq(from=T0,to=T1,by=1)# common time axis interpolated to ...|

ax<-approx(tx,x,method="linear",xout=ti)
#attributes(ax)
ay<-approx(ty,y,method="linear",xout=ti)
ix<-scale(zoo(ax$y,order.by=ti))
iy<-scale(zoo(ay$y,order.by=ti))
}

#browser()
### Compute correlation
C_<-ccf(drop(ix),drop(iy),lag.max=ceiling(max(normlag)),type="correlation",plot=FALSE,na.action=na.omit)
#print(C_$lag)
corr<-C_$acf

### return
return(corr)
}
