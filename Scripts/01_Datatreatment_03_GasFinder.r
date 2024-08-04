
##########################################
##########################################
#####                                #####
#####    Data treatment GasFinder    #####
#####                                #####
##########################################
##########################################

# Author: Marcel Bühler
# Date: August 4, 2024
# Contact: mb@bce.au.dk or Christoph Häni christoph.haeni@bfh.ch
# Description: This script reads in the GasFinder data and makes it ready for further use.
#
# Note: This code was written by Marcel Bühler and is intended to follow the publication 'Applicability of the inverse dispersion method to measure emissions from animal housing' in AMT. 
# Please feel free to use and modify it, but attribution is appreciated.


#################
### Libraries ###
#################

library(ibts)
library(data.table)
library(readxl)


#############
### Paths ###
#############

PathData <- "Path to /data"   
PathRSaves <- "Path to /RSaves"


#################
### Functions ###
#################

source("https://raw.githubusercontent.com/hafl-gel/gel-scripts/main/gasfinder-functions.r")
source(file.path(file.path(dirname(PathRSaves),"Other/shift_dt.r")))

cov_fun <- function(x1,x2){
  n <- length(x1)
  x1 <- fft(x1)/n
  x2 <- fft(x2)/n
  if(n%%2){
    # n odd:
    Re(fft(Conj(x2) * x1, inverse=TRUE))[c(((n+1)/2+1):n,1:((n+1)/2))]*n/(n-1)
  } else {
    # n even:
    Re(fft(Conj(x2) * x1, inverse=TRUE))[c((n/2+1):n,1:(n/2))]*n/(n-1)
  }
}

find_dynlag <- function(x,dyn){
  n <- length(x)
  if(n%%2){
    # odd
    m <- (n+1)/2
  } else {
    # even
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


###################
### time offset ###
###################

time_offset_raw <- read_excel(file.path(dirname(PathRSaves),"Other/time_offset.xlsx"))
time_offset <- as.data.table(time_offset_raw)
setnames(time_offset,c("Inst","RefCET","DeviceTime"))
time_offset[,RefTime := parse_date_time3(RefCET,tz="CET")]
time_offset[,DevTime := parse_date_time3(DeviceTime,tz="Etc/GMT-1")]
setkey(time_offset, Inst)


##############################################################
### load weather station data for Temperature and Pressure ###
##############################################################

### load WS700_2
WS700_2 <- readRDS(file.path(PathRSaves,"WS700.rds"))
# extract only Temp & Press
WS2 <- merge(WS700_2[[1]]$data[,"WS_0160"],WS700_2[[2]]$data[,"WS_0360"])
names(WS2) <- c("T2m_deg","P_hPa")

### load path lengths of GasFinder
load(file.path(PathRSaves,'Geometry.RData'))


########################
### Define Campaigns ###
########################

StartMeas <- "04.03.2021" # to ensure that all data is read in
StopMeas <- "26.03.2021 12:00"
IC1 <- "05.03.2021 to 10.03.2021 18:00"
MC <- "18.03.2021 11:00 - 21.03.2021 14:00"
IC2	<- "21.03.2021 14:00 to "


#####################################
### read GF-30016 or OP-8.6h data ###
#####################################

	File_Path <- paste0(PathData,"/GasFinder/GF16")
	GF_16_Raw <- read.GF3(File_Path,From=StartMeas,To=StopMeas)
	# apply time offset
	GF_16_X <- shift_dt(time_offset["GF16"],GF_16_Raw)
	dummy <- round(st(GF_16_X)[1], "secs")
	st(GF_16_X) <- dummy + round(as.numeric(st(GF_16_X) - dummy), 3)
	et(GF_16_X) <- dummy + round(as.numeric(et(GF_16_X) - dummy), 3)
	# Selection according to status_code (see GasFinder Manual pages 38 to 40)
	# table(GF_16_X[,"status_code"])
	GF_16_X <- GF_16_X[grepl("1$|3$|5$|7$",GF_16_X$"status_code"),c(1,2,4)]
	##### Add temperature and pressure
	GF_16 <- WS2[GF_16_X,c("T2m_deg","P_hPa")][,c(names(GF_16_X),c("T2m_deg","P_hPa"))]
	# rm(GF_16_X)
	##### Add path lengths
	GF_16$PathLength <- pathlengths['IC1',"GF16"]
	GF_16[MC,"PathLength"] <- pathlengths['MC',"GF16"]
	GF_16[IC2,"PathLength"] <- pathlengths['IC2',"GF16"]
	##### calculate mg/m3:
	GF_16$CH4_mgm3 <- GF_16[["CH4 (ppm-m)"]]/GF_16$PathLength*C_pT("GF-16",GF_16[["P_hPa"]],GF_16[["T2m_deg"]],units="mgm3")
	##### calculate ppm:
	GF_16$CH4_ppm <- GF_16[["CH4 (ppm-m)"]]/GF_16$PathLength*C_pT("GF-16",GF_16[["P_hPa"]],GF_16[["T2m_deg"]],units="ppm")


#####################################
### read GF-30017 or OP-2.0h data ###
#####################################

	File_Path <- paste0(PathData,"/GasFinder/GF17")
	GF_17_Raw <- read.GF3(File_Path,From=StartMeas,To=StopMeas)
	# correct time offset
	GF_17_X <- shift_dt(time_offset["GF17"],GF_17_Raw)
	dummy <- round(st(GF_17_X)[1], "secs")
	st(GF_17_X) <- dummy + round(as.numeric(st(GF_17_X) - dummy), 3)
	et(GF_17_X) <- dummy + round(as.numeric(et(GF_17_X) - dummy), 3)
	# Selection according to status_code (see GasFinder Manual pages 38 to 40)
	# table(GF_17_X[,"status_code"])
	GF_17_X <- GF_17_X[grepl("1$|3$|5$|7$",GF_17_X$"status_code"),c(1,2,4)]
	##### Add temperature and pressure
	# GF_17 <- Meteo[GF_17_X,c("T2m_deg","P_hPa")][,c(names(GF_17_X),c("T2m_deg","P_hPa"))]
	GF_17 <- WS2[GF_17_X,c("T2m_deg","P_hPa")][,c(names(GF_17_X),c("T2m_deg","P_hPa"))]
	# rm(GF_17_X)
	##### Add path lengths
	GF_17$PathLength <- pathlengths['IC1',"GF17"]
	GF_17[MC,"PathLength"] <- pathlengths['MC',"GF17"]
	GF_17[IC2,"PathLength"] <- pathlengths['IC2',"GF17"]
	##### calculate mg/m3:
	GF_17$CH4_mgm3 <- GF_17[["CH4 (ppm-m)"]]/GF_17$PathLength*C_pT("GF-17",GF_17[["P_hPa"]],GF_17[["T2m_deg"]],units="mgm3")
	##### calculate ppm:
	GF_17$CH4_ppm <- GF_17[["CH4 (ppm-m)"]]/GF_17$PathLength*C_pT("GF-17",GF_17[["P_hPa"]],GF_17[["T2m_deg"]],units="ppm")


#####################################
### read GF-30018 or OP-5.3h data ###
#####################################

	File_Path <- paste0(PathData,"/GasFinder/GF18")
	GF_18_Raw <- read.GF3(File_Path,From=StartMeas,To=StopMeas)
	# correct time offset 
	GF_18_X <- shift_dt(time_offset["GF18"],GF_18_Raw)
	dummy <- round(st(GF_18_X)[1], "secs")
	st(GF_18_X) <- dummy + round(as.numeric(st(GF_18_X) - dummy), 3)
	et(GF_18_X) <- dummy + round(as.numeric(et(GF_18_X) - dummy), 3)
	# Selection according to status_code (see GasFinder Manual pages 38 to 40)
	# table(GF_18_X[,"status_code"])
	GF_18_X <- GF_18_X[grepl("1$|3$|5$|7$",GF_18_X$"status_code"),c(1,2,4)]
	##### Add temperature and pressure
	# GF_18 <- Meteo[GF_18_X,c("T2m_deg","P_hPa")][,c(names(GF_18_X),c("T2m_deg","P_hPa"))]
	GF_18 <- WS2[GF_18_X,c("T2m_deg","P_hPa")][,c(names(GF_18_X),c("T2m_deg","P_hPa"))]
	# rm(GF_18_X)
	##### Add path lengths
	GF_18$PathLength <- pathlengths['IC1',"GF18"]
	GF_18[MC,"PathLength"] <- pathlengths['MC',"GF18"]
	GF_18[IC2,"PathLength"] <- pathlengths['IC2',"GF18"]
	##### calculate mg/m3:
	GF_18$CH4_mgm3 <- GF_18[["CH4 (ppm-m)"]]/GF_18$PathLength*C_pT("GF-18",GF_18[["P_hPa"]],GF_18[["T2m_deg"]],units="mgm3")
	##### calculate ppm:
	GF_18$CH4_ppm <- GF_18[["CH4 (ppm-m)"]]/GF_18$PathLength*C_pT("GF-18",GF_18[["P_hPa"]],GF_18[["T2m_deg"]],units="ppm")


####################################
### read GF-30025 or OP-12h data ###
####################################

	File_Path <- paste0(PathData,"/GasFinder/GF25")
	GF_25_Raw <- read.GF3(File_Path,From=StartMeas,To=StopMeas)
	# correct time offset
	GF_25_X <- shift_dt(time_offset["GF25"],GF_25_Raw)
	dummy <- round(st(GF_25_X)[1], "secs")
	st(GF_25_X) <- dummy + round(as.numeric(st(GF_25_X) - dummy), 3)
	et(GF_25_X) <- dummy + round(as.numeric(et(GF_25_X) - dummy), 3)
	# Selection according to status_code (see GasFinder Manual pages 38 to 40)
	# table(GF_25_X[,"status_code"])
	GF_25_X <- GF_25_X[grepl("1$|3$|5$|7$",GF_25_X$"status_code"),c(1,2,4)]
	##### Add temperature and pressure
	GF_25 <- WS2[GF_25_X,c("T2m_deg","P_hPa")][,c(names(GF_25_X),c("T2m_deg","P_hPa"))]
	# rm(GF_25_X)
	##### Add path lengths
	GF_25$PathLength <- pathlengths['IC1',"GF25"]
	GF_25[MC,"PathLength"] <- pathlengths['MC',"GF25"]
	GF_25[IC2,"PathLength"] <- pathlengths['IC2',"GF25"]
	##### calculate mg/m3:
	GF_25$CH4_mgm3 <- GF_25[["CH4 (ppm-m)"]]/GF_25$PathLength*C_pT("GF-25",GF_25[["P_hPa"]],GF_25[["T2m_deg"]],units="mgm3")
	##### calculate ppm:
	GF_25$CH4_ppm <- GF_25[["CH4 (ppm-m)"]]/GF_25$PathLength*C_pT("GF-25",GF_25[["P_hPa"]],GF_25[["T2m_deg"]],units="ppm")
		

###################################
### read GF-30026 or OP-UW data ###
###################################

	File_Path <- paste0(PathData,"/GasFinder/GF26")
	GF_26_Raw <- read.GF3(File_Path,From=StartMeas,To=StopMeas)
	# correct time offset
	GF_26_X <- shift_dt(time_offset["GF26"],GF_26_Raw)
	dummy <- round(st(GF_26_X)[1], "secs")
	st(GF_26_X) <- dummy + round(as.numeric(st(GF_26_X) - dummy), 3)
	et(GF_26_X) <- dummy + round(as.numeric(et(GF_26_X) - dummy), 3)
	# Selection according to status_code (see GasFinder Manual pages 38 to 40)
	# table(GF_26_X[,"status_code"])
	GF_26_X <- GF_26_X[grepl("1$|3$|5$|7$",GF_26_X$"status_code"),c(1,2,4)]
	##### Add temperature and pressure
	GF_26 <- WS2[GF_26_X,c("T2m_deg","P_hPa")][,c(names(GF_26_X),c("T2m_deg","P_hPa"))]
	# rm(GF_26_X)
	##### Add path lengths
	GF_26$PathLength <- pathlengths['IC1',"GF26"]
	GF_26[MC,"PathLength"] <- pathlengths['MC',"GF26"]
	GF_26[IC2,"PathLength"] <- pathlengths['IC2',"GF26"]
	##### calculate mg/m3:
	GF_26$CH4_mgm3 <- GF_26[["CH4 (ppm-m)"]]/GF_26$PathLength*C_pT("GF-26",GF_26[["P_hPa"]],GF_26[["T2m_deg"]],units="mgm3")
	##### calculate ppm:
	GF_26$CH4_ppm <- GF_26[["CH4 (ppm-m)"]]/GF_26$PathLength*C_pT("GF-26",GF_26[["P_hPa"]],GF_26[["T2m_deg"]],units="ppm")


###########################################
### correct for really small timeshifts ###
###########################################

dt_days <- 1*24*3600
starts <- seq(parse_date_time3(StartMeas, tz = "Etc/GMT-1"),
	parse_date_time3(StopMeas, tz = "Etc/GMT-1"), by = dt_days)
ends <- starts + dt_days


	##### GF26 vs. GF16 (without MC)
	dyn <- c(-1,1)*200

	xy_2616 <- na.omit(cbind(
		pool(GF_26[-MC, "CH4_mgm3"], "1secs", StartMeas, StopMeas), 
		pool(GF_16[-MC, "CH4_mgm3"], "1secs", StartMeas, StopMeas)
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


	#### GF26 vs GF17 (without MC)
	dyn <- c(-1,1)*200

	xy_2617 <- na.omit(cbind(
		pool(GF_26[-MC, "CH4_mgm3"], "1secs", StartMeas, StopMeas), 
		pool(GF_17[-MC, "CH4_mgm3"], "1secs", StartMeas, StopMeas)
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


	#### GF26 vs GF18 (without MC)
	dyn <- c(-1,1)*200

	xy_2618 <- na.omit(cbind(
		pool(GF_26[-MC, "CH4_mgm3"], "1secs", StartMeas, StopMeas), 
		pool(GF_18[-MC, "CH4_mgm3"], "1secs", StartMeas, StopMeas)
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

	
	#### GF26 vs. GF25 (without MC)
	dyn <- c(-1,1)*200

	xy_2625 <- na.omit(cbind(
		pool(GF_26[-MC, "CH4_mgm3"], "1secs", StartMeas, StopMeas), 
		pool(GF_25[-MC, "CH4_mgm3"], "1secs", StartMeas, StopMeas)
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
	
	
#################
### save data ###
#################

save(GF_16,GF_17,GF_18,GF_25,GF_26,file=file.path(PathRSaves,'GFraw.RData'))

