


#####################################################################################################################
#####################################################################################################################
#####                                                                                                           #####
#####    Function needed to correct the time of the device to a preferred reference time (usually Etc/GMT-1)    #####
#####                                                                                                           #####
#####################################################################################################################
#####################################################################################################################


library(ibts)
library(data.table)

shift_dt <- function(x,ibts,tz="Etc/GMT-1",Plot=FALSE){
	offset <- as.numeric(x$RefTime - x$DevTime,units="secs")
	times <- with_tz(x$RefTime,tz=tz)
	a <- offset[-length(offset)]
	b <- (offset[-1] - a)/as.numeric(times[-1] - times[-length(times)],units="secs")
	st_out <- st_in <- st(ibts)
	et_out <- et_in <- et(ibts)
	ind <- findInterval(st_in,times,all.inside=TRUE)
	for(i in unique(ind)){
	  st_sub <- st_in[ind == i]
	  st_out[ind == i] <- st_sub + round(a[i] + b[i]*as.numeric(st_sub - times[i],units="secs"), 5)
	  et_sub <- et_in[ind == i]
	  et_out[ind == i] <- et_sub + round(a[i] + b[i]*as.numeric(et_sub - times[i],units="secs"), 5)
	}
	if(Plot){
	  par(mfrow=c(2,1))
	  d_st <- st_out-st_in
	  d_et <- et_out-et_in
	  plot(ibts[,1],blank=TRUE,ylim=range(d_st),ylab=attr(d_st,"units"),main="st")
	  lines(st_in,d_st)
	  plot(ibts[,1],blank=TRUE,ylim=range(d_et),ylab=attr(d_et,"units"),main="et")
	  lines(et_in,d_et)
	}
	attr(ibts,"st") <- st_out
	attr(ibts,"et") <- et_out
	ibts
	 }