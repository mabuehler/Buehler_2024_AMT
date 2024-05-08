

#############################################
#############################################
#####                                   #####
#####    Datenaufbereitung GasFinder    #####
#####                                   #####
#############################################
#############################################


#################################
### Header (Pfade, Libraries) ###
#################################

	library(ibts)
	library(data.table)
	library(readxl)

PathData <- "Path to /data"		
PathRSaves <- "Path to /RSaves"
source("https://raw.githubusercontent.com/hafl-gel/gel-scripts/main/gasfinder-functions.r")


###################
### Zeitversatz ###
###################

	Zeitversatz_raw <- read_excel(file.path(dirname(PfadRSaves),"Andere Dateien","Zeitversatz_Stockerematt.xlsx"))
	Zeitversatz <- as.data.table(Zeitversatz_raw)
	setnames(Zeitversatz,c("Inst","RefCET","GF"))
	Zeitversatz[,RefTime := parse_date_time3(RefCET,tz="CET")]
	Zeitversatz[,GFTime := parse_date_time3(GF,tz="Etc/GMT-1")]
	setkey(Zeitversatz, Inst)
  	# Zeitversatz[,.(Inst, RefTime - GFTime)]

	### Funktionen
	shift_dt <- function(x,ibts,tz="Etc/GMT-1",Plot=FALSE){
		versatz <- as.numeric(x$RefTime - x$GFTime,units="secs")
		zeiten <- with_tz(x$RefTime,tz="Etc/GMT-1")
		a <- versatz[-length(versatz)]
		b <- (versatz[-1] - a)/as.numeric(zeiten[-1] - zeiten[-length(zeiten)],units="secs")
		st_out <- st_in <- st(ibts)
		et_out <- et_in <- et(ibts)
		ind <- findInterval(st_in,zeiten,all.inside=TRUE)
		for(i in unique(ind)){
			st_sub <- st_in[ind == i]
			st_out[ind == i] <- st_sub + a[i] + b[i]*as.numeric(st_sub - zeiten[i],units="secs")
			et_sub <- et_in[ind == i]
			et_out[ind == i] <- et_sub + a[i] + b[i]*as.numeric(et_sub - zeiten[i],units="secs")
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


##############################################
### Wetterstation laden für Temp und Press ###
##############################################

	### lade WS700_2
	WS700_2 <- readRDS(paste0(PfadRSaves,"/WS700_2_STO.rds"))
	# extract only Temp & Press
	WS2 <- merge(WS700_2[[1]]$data[,"WS_0160"],WS700_2[[2]]$data[,"WS_0360"])
	names(WS2) <- c("T2m_deg","P_hPa")
	### lade Pfadlängen
	load(file.path(PfadRSaves,"pathLengths_STO.RData"))


################################
### Messkampagnen definieren ###
################################

	StartMeas <- "19.01.2021"
	StopMeas <- "26.03.2021 12:00"
	QV1 <- " to 28.01.2021 12:00"
	QV2 <- "05.03.2021 to 10.03.2021 18:00"
	QV1u2 <- " to 10.03.2021 18:00"
	MK <- "18.03.2021 11:00 - 21.03.2021 14:00"
	QV3	<- "21.03.2021 14:00 to "


###########################
### Lese GF-30016 Daten ###
###########################

	File_Path <- paste0(PfadDaten,"/GasFinder/GF16")
	GF_16_Raw <- read.GF3(File_Path,From=StartMeas,To=StopMeas)
	# Korrektur Zeitversatz
	GF_16_X <- shift_dt(Zeitversatz["GF16"],GF_16_Raw)
	dummy <- round(st(GF_16_X)[1], "secs")
	st(GF_16_X) <- dummy + round(as.numeric(st(GF_16_X) - dummy), 3)
	et(GF_16_X) <- dummy + round(as.numeric(et(GF_16_X) - dummy), 3)
	# Selektion gemäss status_code (siehe GasFinder Manual Seiten 38 bis 40)
	# table(GF_16_X[,"status_code"])
	GF_16_X <- GF_16_X[grepl("1$|3$|5$|7$",GF_16_X$"status_code"),c(1,2,4)]
	##### Füge Temperatur und Druck Daten hinzu
	GF_16 <- WS2[GF_16_X,c("T2m_deg","P_hPa")][,c(names(GF_16_X),c("T2m_deg","P_hPa"))] # Mit Temperatur von Wetterstation
	# rm(GF_16_X)
	##### Füge Pfadlänge hinzu
	GF_16$PathLength <- pathLengths_QV1u2["GF16"]
	GF_16[MK,"PathLength"] <- pathLengths_MK["GF16"]
	GF_16[QV3,"PathLength"] <- pathLengths_QV3["GF16"]
	##### berechne mg/m3:
	GF_16$CH4_mgm3 <- GF_16[["CH4 (ppm-m)"]]/GF_16$PathLength*C_pT("GF-16",GF_16[["P_hPa"]],GF_16[["T2m_deg"]],units="mgm3")
	##### berechne ppm:
	GF_16$CH4_ppm <- GF_16[["CH4 (ppm-m)"]]/GF_16$PathLength*C_pT("GF-16",GF_16[["P_hPa"]],GF_16[["T2m_deg"]],units="ppm")


###########################
### Lese GF-30017 Daten ###
###########################

	File_Path <- paste0(PfadDaten,"/GasFinder/GF17")
	GF_17_Raw <- read.GF3(File_Path,From=StartMeas,To=StopMeas)
	# Korrektur Zeitversatz
	GF_17_X <- shift_dt(Zeitversatz["GF17"],GF_17_Raw)
	dummy <- round(st(GF_17_X)[1], "secs")
	st(GF_17_X) <- dummy + round(as.numeric(st(GF_17_X) - dummy), 3)
	et(GF_17_X) <- dummy + round(as.numeric(et(GF_17_X) - dummy), 3)
	# Selektion gemäss status_code (siehe GasFinder Manual Seiten 38 bis 40)
	# table(GF_17_X[,"status_code"])
	GF_17_X <- GF_17_X[grepl("1$|3$|5$|7$",GF_17_X$"status_code"),c(1,2,4)]
	##### Füge Temperatur und Druck Daten hinzu
	# GF_17 <- Meteo[GF_17_X,c("T2m_deg","P_hPa")][,c(names(GF_17_X),c("T2m_deg","P_hPa"))] # Mit Temperatur von Wetterstation
	GF_17 <- WS2[GF_17_X,c("T2m_deg","P_hPa")][,c(names(GF_17_X),c("T2m_deg","P_hPa"))] # Mit Temperatur von Wetterstation
	# rm(GF_17_X)
	##### Füge Pfadlänge hinzu
	GF_17$PathLength <- pathLengths_QV1u2["GF17"]
	GF_17[MK,"PathLength"] <- pathLengths_MK["GF17"]
	GF_17[QV3,"PathLength"] <- pathLengths_QV3["GF17"]
	##### berechne mg/m3:
	GF_17$CH4_mgm3 <- GF_17[["CH4 (ppm-m)"]]/GF_17$PathLength*C_pT("GF-17",GF_17[["P_hPa"]],GF_17[["T2m_deg"]],units="mgm3")
	##### berechne ppm:
	GF_17$CH4_ppm <- GF_17[["CH4 (ppm-m)"]]/GF_17$PathLength*C_pT("GF-17",GF_17[["P_hPa"]],GF_17[["T2m_deg"]],units="ppm")


###########################
### Lese GF-30018 Daten ###
###########################

	File_Path <- paste0(PfadDaten,"/GasFinder/GF18")
	GF_18_Raw <- read.GF3(File_Path,From=StartMeas,To=StopMeas)
	# Korrektur Zeitversatz 
	GF_18_X <- shift_dt(Zeitversatz["GF18"],GF_18_Raw)
	dummy <- round(st(GF_18_X)[1], "secs")
	st(GF_18_X) <- dummy + round(as.numeric(st(GF_18_X) - dummy), 3)
	et(GF_18_X) <- dummy + round(as.numeric(et(GF_18_X) - dummy), 3)
	# Selektion gemäss status_code (siehe GasFinder Manual Seiten 38 bis 40)
	# table(GF_18_X[,"status_code"])
	GF_18_X <- GF_18_X[grepl("1$|3$|5$|7$",GF_18_X$"status_code"),c(1,2,4)]
	##### Füge Temperatur und Druck Daten hinzu
	# GF_18 <- Meteo[GF_18_X,c("T2m_deg","P_hPa")][,c(names(GF_18_X),c("T2m_deg","P_hPa"))] # Mit Temperatur von Wetterstation
	GF_18 <- WS2[GF_18_X,c("T2m_deg","P_hPa")][,c(names(GF_18_X),c("T2m_deg","P_hPa"))] # Mit Temperatur von Wetterstation
	# rm(GF_18_X)
	##### Füge Pfadlänge hinzu
	GF_18$PathLength <- pathLengths_QV1u2["GF18"]
	GF_18[MK,"PathLength"] <- pathLengths_MK["GF18"]
	GF_18[QV3,"PathLength"] <- pathLengths_QV3["GF18"]
	##### berechne mg/m3:
	GF_18$CH4_mgm3 <- GF_18[["CH4 (ppm-m)"]]/GF_18$PathLength*C_pT("GF-18",GF_18[["P_hPa"]],GF_18[["T2m_deg"]],units="mgm3")
	##### berechne ppm:
	GF_18$CH4_ppm <- GF_18[["CH4 (ppm-m)"]]/GF_18$PathLength*C_pT("GF-18",GF_18[["P_hPa"]],GF_18[["T2m_deg"]],units="ppm")


###########################
### Lese GF-30025 Daten ###
###########################

	File_Path <- paste0(PfadDaten,"/GasFinder/GF25")
	GF_25_Raw <- read.GF3(File_Path,From=StartMeas,To=StopMeas)
	# Korrektur Zeitversatz
	GF_25_X <- shift_dt(Zeitversatz["GF25"],GF_25_Raw)
	dummy <- round(st(GF_25_X)[1], "secs")
	st(GF_25_X) <- dummy + round(as.numeric(st(GF_25_X) - dummy), 3)
	et(GF_25_X) <- dummy + round(as.numeric(et(GF_25_X) - dummy), 3)
	# Selektion gemäss status_code (siehe GasFinder Manual Seiten 38 bis 40)
	# table(GF_25_X[,"status_code"])
	GF_25_X <- GF_25_X[grepl("1$|3$|5$|7$",GF_25_X$"status_code"),c(1,2,4)]
	##### Füge Temperatur und Druck Daten hinzu
	GF_25 <- WS2[GF_25_X,c("T2m_deg","P_hPa")][,c(names(GF_25_X),c("T2m_deg","P_hPa"))] # Mit Temperatur von Wetterstation
	# rm(GF_25_X)
	##### Füge Pfadlänge hinzu
	GF_25$PathLength <- pathLengths_QV1u2["GF25"]
	GF_25[MK,"PathLength"] <- pathLengths_MK["GF25"]
	GF_25[QV3,"PathLength"] <- pathLengths_QV3["GF25"]
	##### berechne mg/m3:
	GF_25$CH4_mgm3 <- GF_25[["CH4 (ppm-m)"]]/GF_25$PathLength*C_pT("GF-25",GF_25[["P_hPa"]],GF_25[["T2m_deg"]],units="mgm3")
	##### berechne ppm:
	GF_25$CH4_ppm <- GF_25[["CH4 (ppm-m)"]]/GF_25$PathLength*C_pT("GF-25",GF_25[["P_hPa"]],GF_25[["T2m_deg"]],units="ppm")
		

###########################
### Lese GF-30026 Daten ###
###########################

	File_Path <- paste0(PfadDaten,"/GasFinder/GF26")
	GF_26_Raw <- read.GF3(File_Path,From=StartMeas,To=StopMeas)
	# Korrektur Zeitversatz
	GF_26_X <- shift_dt(Zeitversatz["GF26"],GF_26_Raw)
	dummy <- round(st(GF_26_X)[1], "secs")
	st(GF_26_X) <- dummy + round(as.numeric(st(GF_26_X) - dummy), 3)
	et(GF_26_X) <- dummy + round(as.numeric(et(GF_26_X) - dummy), 3)
	# Selektion gemäss status_code (siehe GasFinder Manual Seiten 38 bis 40)
	# table(GF_26_X[,"status_code"])
	GF_26_X <- GF_26_X[grepl("1$|3$|5$|7$",GF_26_X$"status_code"),c(1,2,4)]
	##### Füge Temperatur und Druck Daten hinzu
	GF_26 <- WS2[GF_26_X,c("T2m_deg","P_hPa")][,c(names(GF_26_X),c("T2m_deg","P_hPa"))] # Mit Temperatur von Wetterstation
	# rm(GF_26_X)
	##### Füge Pfadlänge hinzu
	GF_26$PathLength <- pathLengths_QV1u2["GF26"]
	GF_26[MK,"PathLength"] <- pathLengths_MK["GF26"]
	GF_26[QV3,"PathLength"] <- pathLengths_QV3["GF26"]
	##### berechne mg/m3:
	GF_26$CH4_mgm3 <- GF_26[["CH4 (ppm-m)"]]/GF_26$PathLength*C_pT("GF-26",GF_26[["P_hPa"]],GF_26[["T2m_deg"]],units="mgm3")
	##### berechne ppm:
	GF_26$CH4_ppm <- GF_26[["CH4 (ppm-m)"]]/GF_26$PathLength*C_pT("GF-26",GF_26[["P_hPa"]],GF_26[["T2m_deg"]],units="ppm")


####################
### Time Shifts ####
####################
		
	# GF_16 <- readRDS(file.path(PfadRSaves,"STO_GF_16.rds"))
	# GF_17 <- readRDS(file.path(PfadRSaves,"STO_GF_17.rds"))
	# GF_18 <- readRDS(file.path(PfadRSaves,"STO_GF_18.rds"))
	# GF_25 <- readRDS(file.path(PfadRSaves,"STO_GF_25.rds"))
	# GF_26 <- readRDS(file.path(PfadRSaves,"STO_GF_26.rds"))

	#### check time shifts:
	cov_fun <- function(x1,x2){
	  n <- length(x1)
	  x1 <- fft(x1)/n
	  x2 <- fft(x2)/n
	  if(n%%2){
	    # n ungerade:
	    Re(fft(Conj(x2) * x1, inverse=TRUE))[c(((n+1)/2+1):n,1:((n+1)/2))]*n/(n-1)
	  } else {
	    # n gerade:
	    Re(fft(Conj(x2) * x1, inverse=TRUE))[c((n/2+1):n,1:(n/2))]*n/(n-1)
	  }
	}

	find_dynlag <- function(x,dyn){
	  n <- length(x)
	  if(n%%2){
	    # ungerade
	    m <- (n+1)/2
	  } else {
	    # gerade
	    m <- n/2 + 1
	  }
	  ind <- seq(dyn[1],dyn[2]) + m
	  # find max:
	  maxis <- ind[which.max(abs(x[ind]))]
	  maxis - m
	}

	checkTimeDiff <- function(xy,dyn = c(-1, 1)*1E4, Plot = FALSE, xlim = dl + c(-1, 1)*50, ...){
	  dt_xy <- as.numeric(median(diff(st(xy))), units = "secs") 
	  cv <- cov_fun(xy[[1]], xy[[2]])
	  n <- length(cv)
	  if((n %% 2) == 1) { 
	    xx <- seq(-(n - 1)/2,(n - 1)/2,by=dt_xy) * dt_xy
	  } else {
	    xx <- seq(-n/2 + 1,n/2,by=dt_xy) * dt_xy
	  }
	  dl <- find_dynlag(cv, dyn)
	  if(Plot){
	    plot(xx, cv, type = "l", xlim = xlim, ...)
	  }
	  dl * dt_xy
	}

	dt_days <- 1*24*3600
	starts <- seq(parse_date_time3(StartMeas, tz = "Etc/GMT-1"),
		parse_date_time3(StopMeas, tz = "Etc/GMT-1"), by = dt_days)
	ends <- starts + dt_days

	#########################

	##### GF26 vs. GF16 (ohne MK)
	dyn <- c(-1,1)*200

	xy_2616 <- na.omit(cbind(
		pool(GF_26[-MK, "CH4_mgm3"], "1secs", StartMeas, StopMeas), 
		pool(GF_16[-MK, "CH4_mgm3"], "1secs", StartMeas, StopMeas)
		))

	ts_16 <- numeric(length(starts)) * NA_real_

	for(i in seq_along(starts)){
	  ind <- deparse_timerange(starts[i], ends[i])
	  if(nrow(xy_2616[ind]) > diff(dyn)) ts_16[i] <- checkTimeDiff(xy_2616[ind], dyn = dyn)
	}
	median(ts_16, na.rm = TRUE)
	ts_16

	xs <- starts+12*3600
	mod16 <- MASS::rlm(ts_16 ~ xs)
	x11()
	plot(xs, ts_16)
	abline(mod16, col = "red")

	# shift GF-16
	GF_16_backup <- GF_16
	attr(GF_16, "st") <- attr(GF_16, "st") + predict(mod16, list(xs=attr(GF_16, "st")))
	attr(GF_16, "et") <- attr(GF_16, "et") + predict(mod16, list(xs=attr(GF_16, "et")))
	
	#########################


	#### GF26 vs GF17 (ohne MK)
	dyn <- c(-1,1)*200

	xy_2617 <- na.omit(cbind(
		pool(GF_26[-MK, "CH4_mgm3"], "1secs", StartMeas, StopMeas), 
		pool(GF_17[-MK, "CH4_mgm3"], "1secs", StartMeas, StopMeas)
		))

	ts_17 <- numeric(length(starts)) * NA_real_

	for(i in seq_along(starts)){
	  ind <- deparse_timerange(starts[i], ends[i])
	  if(nrow(xy_2617[ind]) > diff(dyn)) ts_17[i] <- checkTimeDiff(xy_2617[ind], dyn = dyn)
	}
	median(ts_17, na.rm = TRUE)
	ts_17

	xs <- starts+12*3600
	mod17 <- MASS::rlm(ts_17 ~ xs)
	x11()
	plot(xs, ts_17)
	abline(mod17, col = "red")

	# shift GF-17
	GF_17_backup <- GF_17
	attr(GF_17, "st") <- attr(GF_17, "st") + predict(mod17, list(xs=attr(GF_17, "st")))
	attr(GF_17, "et") <- attr(GF_17, "et") + predict(mod17, list(xs=attr(GF_17, "et")))


	#########################


	#### GF26 vs GF18 (ohne MK)
	dyn <- c(-1,1)*200

	xy_2618 <- na.omit(cbind(
		pool(GF_26[-MK, "CH4_mgm3"], "1secs", StartMeas, StopMeas), 
		pool(GF_18[-MK, "CH4_mgm3"], "1secs", StartMeas, StopMeas)
		))

	ts_18 <- numeric(length(starts)) * NA_real_

	for(i in seq_along(starts)){
	  ind <- deparse_timerange(starts[i], ends[i])
	  if(nrow(xy_2618[ind]) > diff(dyn)) ts_18[i] <- checkTimeDiff(xy_2618[ind], dyn = dyn)
	}
	median(ts_18, na.rm = TRUE)
	ts_18

	xs <- starts+12*3600
	mod18 <- MASS::rlm(ts_18 ~ xs)
	x11()
	plot(xs, ts_18)
	abline(mod18, col = "red")

	# shift GF-18
	GF_18_backup <- GF_18
	attr(GF_18, "st") <- attr(GF_18, "st") + predict(mod18, list(xs=attr(GF_18, "st")))
	attr(GF_18, "et") <- attr(GF_18, "et") + predict(mod18, list(xs=attr(GF_18, "et")))
	
	#########################

	
	#### GF26 vs. GF25 (ohne MK)
	dyn <- c(-1,1)*200

	xy_2625 <- na.omit(cbind(
		pool(GF_26[-MK, "CH4_mgm3"], "1secs", StartMeas, StopMeas), 
		pool(GF_25[-MK, "CH4_mgm3"], "1secs", StartMeas, StopMeas)
		))

	ts_25 <- numeric(length(starts)) * NA_real_

	for(i in seq_along(starts)){
	  ind <- deparse_timerange(starts[i], ends[i])
	  if(nrow(xy_2625[ind]) > diff(dyn)) ts_25[i] <- checkTimeDiff(xy_2625[ind], dyn = dyn)
	}
	median(ts_25, na.rm = TRUE)
	ts_25

	xs <- starts+12*3600
	mod25 <- MASS::rlm(ts_25 ~ xs)
	x11()
	plot(xs, ts_25)
	abline(mod25, col = "red")

	# shift GF-25
	GF_25_backup <- GF_25
	attr(GF_25, "st") <- attr(GF_25, "st") + predict(mod25, list(xs=attr(GF_25, "st")))
	attr(GF_25, "et") <- attr(GF_25, "et") + predict(mod25, list(xs=attr(GF_25, "et")))
	
	
#########################
### Daten abspeichern ###
#########################

	# Daten sichern:
	saveRDS(GF_16,file=paste0(PfadRSaves,"/STO_GF_16.rds"))
	# Daten sichern:
	saveRDS(GF_17,file=paste0(PfadRSaves,"/STO_GF_17.rds"))
	# Daten sichern:
	saveRDS(GF_18,file=paste0(PfadRSaves,"/STO_GF_18.rds"))
	# Daten sichern:
	saveRDS(GF_25,file=paste0(PfadRSaves,"/STO_GF_25.rds"))
	# Daten sichern:
	saveRDS(GF_26,file=paste0(PfadRSaves,"/STO_GF_26.rds"))

