
######################################################
######################################################
#####                                            #####
#####    Data treatment weather station WS700    #####
#####                                            #####
######################################################
######################################################



###########################
### libraries and Paths ###
###########################
	
	library(ibts)
	library(data.table)
	library(readxl)

	PathData <- "Path to /data"		
	PathRSaves <- "Path to /RSaves"
	source('https://raw.githubusercontent.com/hafl-gel/gel-scripts/main/weatherstation-functions.r')


#################
### Functions ###
#################

	shift_dt <- function(x,ibts,tz="Etc/GMT-1",Plot=FALSE){
		versatz <- as.numeric(x$RefTime - x$GFTime,units="secs")
		zeiten <- with_tz(x$RefTime,tz=tz)
		a <- versatz[-length(versatz)]
		b <- (versatz[-1] - a)/as.numeric(zeiten[-1] - zeiten[-length(zeiten)],units="secs")
		st_out <- st_in <- st(ibts)
		et_out <- et_in <- et(ibts)
		ind <- findInterval(st_in,zeiten,all.inside=TRUE)
		for(i in unique(ind)){
		  st_sub <- st_in[ind == i]
		  st_out[ind == i] <- st_sub + round(a[i] + b[i]*as.numeric(st_sub - zeiten[i],units="secs"), 5)
		  et_sub <- et_in[ind == i]
		  et_out[ind == i] <- et_sub + round(a[i] + b[i]*as.numeric(et_sub - zeiten[i],units="secs"), 5)
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

###################
### time offset ###
###################

	time_offset_raw <- read_excel(file.path(dirname(PathRSaves),"Other/time_offsets.xlsx"))
	time_offset <- as.data.table(time_offset_raw)
	setnames(time_offset,c("Inst","RefCET","GF"))
	time_offset[,RefTime := parse_date_time3(RefCET,tz="CET")]
	time_offset[,GFTime := parse_date_time3(GF,tz="Etc/GMT-1")]
	setkey(time_offset, Inst)

####################
### read WS700_2 ###
####################

	File_Path <- paste0(PathData,"/Wetterstation/WS2")
	WS700_2 <- readOML(File_Path, as_ibts = TRUE)
	## WDs with 0.0 (vectorial wind speed WS_0480 0.0) set to NA.
	# WS700_2[[1]]$data[which(WS700_2[[1]]$data[,"WS_0580"] == 0.0)]
	WS700_2[[1]]$data[which(WS700_2[[1]]$data[,"WS_0480"] == 0.0),"WS_0580"] <- NA_real_
	## Deklination korrigieren https://www.swisstopo.admin.ch/de/online/calculation-services/declination.html
	## 20.01.2021 = 2.81168, 26.03.2021 = 2.84343 ---> Der magnetische NP liegt 2.8° östlich vom True North. östlich ist positiv, westilich ist negativ.
	# Soll also eine auf magnetisch Nord bezogene Richtungsangabe in eine auf geographisch Nord bezogene Richtung umgewandelt werden (falsch → wahr),
	# so ist die Deklination unter Berücksichtigung ihres angegebenen Vorzeichens zu addieren. Soll ein auf geographisch Nord bezogener Kurs in einen 
	# auf magnetisch Nord bezogenen Kurs umgewandelt werden (wahr → falsch), so ist die Deklination unter Umkehrung ihres Vorzeichens zu addieren.
	## ---> 2.8° addieren
	WS700_2[[3]]$data$Mag_north_corr <-	(WS700_2[[3]]$data + 2.8) %% 360
	# Spalte für korigierte WD einfügen
	WS700_2[[1]]$data$WD_corr <- NA_real_
	# WD Korrektur Vorsondierung
	MKVor <- "07.12.2020 00:00 - 15.12.2020 00:00"
	WD_Offset_WS2_MKVor <- 360 - mean(WS700_2[[3]]$data[MKVor]$Mag_north_corr,na.rm=TRUE)
	WS700_2[[1]]$data[MKVor]$WD_corr <- (WS700_2[[1]]$data[MKVor]$WS_0580 - WD_Offset_WS2_MKVor) %% 360
	# WD Korrektur QV
	MKQV <- "18.01.2021 00:00 - 26.01.2021 12:00"
	WD_Offset_WS2_MKQV <- 360 - mean(WS700_2[[3]]$data[MKQV]$Mag_north_corr,na.rm=TRUE)
	WS700_2[[1]]$data[MKQV]$WD_corr <- (WS700_2[[1]]$data[MKQV]$WS_0580 - WD_Offset_WS2_MKQV) %% 360
	# WD Korrektur MKQV # da sich alles um 1° gedreht hat.
	MKMK <- "26.01.2021 12:00 - 26.03.2021 12:00"
	WD_Offset_WS2_MKMK <- 360 - mean(WS700_2[[3]]$data[MKMK]$Mag_north_corr,na.rm=TRUE)
	WS700_2[[1]]$data[MKMK]$WD_corr <- (WS700_2[[1]]$data[MKMK]$WS_0580 - WD_Offset_WS2_MKMK) %% 360
	#
	colClasses(WS700_2[[1]]$data)[c("WS_0580","WD_corr")] <- "circ"
	WS700_2_d1 <- shift_dt(time_offset["WS2"], WS700_2[[1]]$data)
	WS700_2_d2 <- shift_dt(time_offset["WS2"], WS700_2[[2]]$data)
	WS700_2_d3 <- shift_dt(time_offset["WS2"], WS700_2[[3]]$data)
	WS2 <- WS700_2
	WS2[[1]]$data <- WS700_2_d1["07.12.2020 - 26.03.2021 12:00"]
	WS2[[2]]$data <- WS700_2_d2["07.12.2020 - 26.03.2021 12:00"]
	WS2[[3]]$data <- WS700_2_d3["07.12.2020 - 26.03.2021 12:00"]

#################
### save data ###
#################

	saveRDS(WS2,file=paste0(PathRSaves,"/WS700.rds"))

