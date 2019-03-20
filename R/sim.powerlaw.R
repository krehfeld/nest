#' Simulate a power law time series
#' 
#' Simulate a power-law time series
#' 
#' 
#' @usage sim.powerlaw(beta, N)
#' @param beta power-law exponent
#' @param N length of the time series
#' @note The Filter used in the power law generation is currently unscaled. To
#' force Variance(X)=1, a scaling is performed on the final timeseries. This is
#' problematic for short time series and if the distribution of variances is of
#' interest. A workaround is to generate a very long powerlaw series (with
#' overall variance=1) and cut it into many shorter ones.
#' x<-matrix(sim.powerlaw(betas[b],Nsur*Nlength),ncol=Nsur,nrow=Nlength)
#' @keywords ~kwd1 ~kwd2
sim.powerlaw <-
function(beta,N)
{
#time series with power law PSD. variance 1, slope beta, length N
    Norg<-N
    N<-ceiling(N/2)*2
    df  = 1/(N);
    f=seq(from=df,to=1/(2),by=df)
    Filter=sqrt(1/(f^beta));
    Filter = c(max(Filter), Filter,rev(Filter))
    #  Filter = c(Filter,rev(Filter))
    x   = scale(rnorm(N+1,1))  ### Check if the scale is necessary here
    fx  =fft(x)
    ffx =fx*Filter;
    result<-scale(Re(fft(ffx,inverse=TRUE))) 
    ### To look at a realistic distribution, this scaling has to be removed - however, at present the Filter is unscaled - therefore the variance explodes.
    return(result[1:Norg])
}
