#' Periodogram from autocorrelation function
#' 
#' Periodogram from autocorrelation function
#' 
#' 
#' @usage pspek(acf, dt = 1, norm = FALSE, varin = 1)
#' @param acf %% ~~Describe \code{acf} here~~
#' @param dt %% ~~Describe \code{dt} here~~
#' @param norm %% ~~Describe \code{norm} here~~
#' @param varin %% ~~Describe \code{varin} here~~
pspek <-
function(acf,dt=1,norm=FALSE,varin=1){
# compute psd from acf

NFFT=length(acf)


Fs=1/dt;#sampling frequency
T=1/Fs;#sampling period

Y=fft(acf)*Fs #multiply w/ sampling freq

f=Fs*seq(from=0,to=0.5,length.out=floor(NFFT/2)+1)[-1]
#f=Fs*seq(from=0,to=1,length.out=NFFT/2+1)[-1]

power<-abs(Y)[1:floor(NFFT/2)+1]
# mean(power for freq=0:0.5)=0.5*variance(X)

if (norm==TRUE) power<-power/mean(power)*varin*dt


return(list(p=power,f=f,norm=norm))

}
