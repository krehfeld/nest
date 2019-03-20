#quality_check_data

# INPUT: prxlist, start, end dates, maxhiat, length.min



#' Quality checks and convenient splitting of time series
#' 
#' Quality checks and convenient splitting of time series
#' 
#' 
#' @usage quality_check(prxlist, meta, T0, T1, maxhiat = 800, length.min = 15,
#' min.res = 500, min.range=(T1-T0)/3)
#' @param prxlist List of zoo objects
#' @param meta metadata
#' @param T0 Vector of start points
#' @param T1 Vector of end points
#' @param maxhiat Maximum hiatus length in time series
#' @param length.min Minimum acceptable length of the time series
#' @param min.res Minimum acceptable resolution
#' @param min.range Minimum time covered (in time units)
#' @return Returns a list of \item{WBT.out}{Window-by-time-object}
#' \item{metafilt}{Filtered metadata object} \item{ind.ok}{Indices of
#' acceptable objects}
#' @note TBD: Floating window definition (e.g. Eemian by thermal max. in
#' window)
#' @seealso \code{\link{window_by_time}}
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' 
#' # Generate a list of long time series time series, one for each network node
#' Xlist<-lapply(replicate(3,zoo(x=rnorm(1000)*seq(1,1000)),simplify=FALSE),sample,250)
#' metafake<-data.frame(c(-75,30,75),c(120,120,120),c("a","b","c"))
#' colnames(metafake)<-c("Lat","Lon","Name")
#' 
#' # window each time series
#' #OUT<-lapply(Xlist,window_by_time,T0=0,T1=1000,shift=0.75)
#' 
#' OUT<-window_by_time(Xlist,T0=c(0,500),T1=c(500,1000))
#' 
#' OUT<-quality_check(Xlist,metafake,T0=0,T1=1000,maxhiat=10,length.min=100,min.res=10)
#' 
#' 
#' @export quality_check
quality_check<-function(prxlist,meta,T0,T1,maxhiat=800,length.min=15,min.res=500,min.range=(T1-T0)/3){

# min.range<-(T1-T0)/3

make_index_unique<-function(zoobj){
if ((length(zoobj)==0)||is.na(zoobj)) return(zoobj)
if (!is.zoo(zoobj)) zoobj<-as.zoo(zoobj)
tx<-index(zoobj)
dind<-which(diff(index(zoobj))==0)
if (length(dind)>0){
zoobj<-zoobj[-dind]
}
if (any(is.na(zoobj))){zoobj<-na.omit(zoobj)}
if (any(is.na(index(zoobj)))){ind<-which(is.na(index(zoobj)));zoobj<-zoobj[-ind]}

if (mean(diff(index(zoobj)))<0) {warning("decreasing index - check time definition")}
return(zoobj)
}

prxlist<-lapply(prxlist,make_index_unique)


##### _______________

# check if a time series is consistently well sampled or if it has a hiatus of a length > maxhiat. Return the time series for the longest segment if this is sufficiently long (minlength).
cut.hiatus<-function(ts.temp,maxhiat,minlength,forceoneseg=FALSE){
  res<-rle(diff(index(ts.temp))<maxhiat) # find the longest-running stretch with a maximum hiatus of maxhiat
  
  if (all(res$lengths<minlength)) {
  
    ts.out<-NULL
    length.ts<-NA
    hiatus<-NA
    res.out<-NA
  } else { # there is at least one segment which is longer than minlength
  
    id.seg<-which(max(res$lengths)==res$lengths&res$values==TRUE) # find which segment is the longest
    if (length(id.seg)==0) {return(list(hiat=NA,ts.out=NA,length=NA,res=NA))}

    if (length(res$lengths)==1) {# if we have one consecutive time series
      hiatus=FALSE # no hiatus
      if (res$lengths[1]>=length.min){ # if this one piece is long enough
	ts.out<-ts.temp
	length.ts<-length(ts.temp)
	res.out<-mean(diff(index(ts.out)))
      } else { ts.out<-NULL;length.ts<-res$lengths[1];res.out<-NA } # if it isn't long enough, don't spit out the ts at all
    }
    else { # if we have a spacing that is larger than we want (i.e. a hiatus)
      hiatus=TRUE;
      # check if we have a piece that is long enough
      if (res$lengths[id.seg]>=minlength){ #if yes, return the section that is longest
	# indices for starting point
	temp<- c(0,cumsum(res$lengths))+1
	ts.out<-ts.temp[temp[id.seg]:temp[id.seg+1]]
	length.ts<-res$lengths[id.seg]
	res.out<-mean(diff(index(ts.out)))
      } else { #ts not long enough -- don't return the longest segment
	ts.out<-NULL
	length.ts<-res$lengths[id.seg]
	res.out<-NA
      }
    }
  }
  outlist<-list(hiat=hiatus,ts.out=ts.out,length=length.ts,res=res.out)
  return(outlist)
}


####_________________



if (!length(T0)==length(T1)) stop("start/end dates must be the same")
# 
# for (i in 1:length(T0)){
# # cut the time series windows
# 
# }
WBT<-window_by_time(prxlist,T0=T0,T1=T1,shift=NULL,wwidth=NULL)
WBT.filt<-list()

WBT.strip<-lapply(WBT,function(x){x$tsplit}) #list of time series only

chk.res<-ind.ok<-chk.range<-matrix(NA,nrow=length(WBT),ncol=length(T0))
# now find the n for which all k windows satisfy the criteria
# WBT[[n]]$tsplit[[k]]
for (n in 1:length(WBT.strip)){
#ts.hol.in<-lapply(prxlist,window,start=start.hol,end=end.hol)
WBT.filt[[n]]<-lapply(WBT.strip[[n]],cut.hiatus,maxhiat=800,minlength=length.min)
chk.res[n,]<-sapply(WBT.filt[[n]],function(xl){xl$res})<=min.res #& sapply(temp.hol.filt,function(xl){xl$res})<=min.res

chk.range[n,]<-as.numeric(as.character(sapply(WBT.filt[[n]],function(xl){r<-(range(index(xl$ts.out))[2]-range(index(xl$ts.out))[1])})))>min.range
chk.range[is.infinite(chk.range)]<-NA

#ind.ok<-which(chk.res&chk.length)
}


# find all which satisfy the criteria in all subwindows
ind.ok<-which(apply(chk.res,1,all)&apply(chk.range,1,all))

WBT.out<-lapply(WBT.filt,lapply,function(xl){xl$ts.out})
tmid=T0+(T1-T0)/2
WBT.out<-lapply(WBT.filt,function(X,tmid){list(tmid=tmid,tsplit=lapply(X,function(xl){xl$ts.out}))},tmid=tmid)
#lapply(WBT.out,function(x){list(tsplit=x,tmid=})

names(WBT.out)<-meta[["Name"]]
class(WBT.out)<-"window-by-time-object"

#ts.lgm<-lapply(temp.lgm.filt[ind.ok],function(xl){xl$ts.out})
#ts.hol<-lapply(temp.hol.filt[ind.ok],function(xl){xl$ts.out})
#ts.all<-prxlist[ind.ok]
metafilt<-meta[ind.ok,]


return(list(WBT.out=WBT.out,metafilt=metafilt,ind.ok=ind.ok))

}

