
gaussdetr<-function(X, tsc.in = mean(diff(index(X)))*10,prune=FALSE) 
{
h<-tsc.in/6
if (!is.zoo(X)){stop("time series must be zoo object")}

tx<-index(X)
x<-coredata(X)
DT<-mean(diff(index(X)))

    if (is.ts(x)) {
        Xsmooth = detr = x
    }
    else {
        Xsmooth = rnorm(length(tx))
        detr <- Xsmooth
    }
    Xsmooth[1] = mean(x[1:floor(h/2)])
    for (k in 2:length(tx)) {
        WL = 1/sqrt(2 * pi * h^2) * (exp(-((tx - tx[k])^2)/(2 * 
            h^2)))
        Xsmooth[k] = sum(WL * x)/sum(WL)
    }
    Xsmooth[1] = Xsmooth[2]
    
    if (prune){
    Xsmooth[c(1:ceiling(h/(DT)),length(Xsmooth):(length(Xsmooth)-ceiling(h/(DT))))]<-NA
    }

    detr <- x - Xsmooth

    

    result <- list(Xsmooth = zoo(Xsmooth,order.by=tx), detr = zoo(detr,order.by=tx))
    class(result) <- "smoothingobject"
    return(result)
}
