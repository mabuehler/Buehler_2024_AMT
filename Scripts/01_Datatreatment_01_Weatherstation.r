
######################################################
######################################################
#####                                            #####
#####    Data treatment weather station WS700    #####
#####                                            #####
######################################################
######################################################

# Author: Marcel Bühler
# Date: August 4, 2024
# Contact: mb@bce.au.dk or Christoph Häni christoph.haeni@bfh.ch
# Description: This script reads in the weather station data and makes it ready for further use.
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

source('https://raw.githubusercontent.com/hafl-gel/gel-scripts/main/weatherstation-functions.r')
source(file.path(file.path(dirname(PathRSaves),"Other/shift_dt.r")))


###################
### time offset ###
###################

time_offset_raw <- read_excel(file.path(dirname(PathRSaves),"Other/time_offset.xlsx")) # read in file with the different time offsets
time_offset <- as.data.table(time_offset_raw)
setnames(time_offset,c("Inst","RefCET","DeviceTime"))
time_offset[,RefTime := parse_date_time3(RefCET,tz="CET")]
time_offset[,DevTime := parse_date_time3(DeviceTime,tz="Etc/GMT-1")]
setkey(time_offset, Inst)


####################################
### read data and data treatment ###
####################################

File_Path <- paste0(PathData,"/Weatherstation/WS2/")
WS700_2 <- readOML(File_Path, as_ibts = TRUE) # note, this file contains more data than needed. The correct data will be selected in line 70-72.
  ## WDs with 0.0 (vectorial wind speed WS_0480 0.0) set to NA.
  ## WS700_2[[1]]$data[which(WS700_2[[1]]$data[,"WS_0580"] == 0.0)]
WS700_2[[1]]$data[which(WS700_2[[1]]$data[,"WS_0480"] == 0.0),"WS_0580"] <- NA_real_
  ## Correct declination with the help of the following website: https://www.swisstopo.admin.ch/de/online/calculation-services/declination.html
  ## 26.03.2021 = 2.84343 ---> The magnetic north pole lays 2.8° east of true north. Estward is positive, westward is negative. 
  ## "If a direction indicating magnetically North is to be converted into a direction indicating geographically North (false → true), 
  ## the declination has to be added taking into account its given sign. If a course related to geographic north is to be converted to
  ## a course related to magnetic north (true → false), the declination is to be added by reversing its sign."
  ## ---> adding 2.8°
WS700_2[[3]]$data$Mag_north_corr <-	(WS700_2[[3]]$data + 2.8) %% 360
  ## Add column for corrected WD
WS700_2[[1]]$data$WD_corr <- NA_real_
i_period <- "05.03.2021 - 26.03.2021 12:00"
WD_Offset_WS2 <- 360 - mean(WS700_2[[3]]$data[i_period]$Mag_north_corr,na.rm=TRUE)
WS700_2[[1]]$data[i_period]$WD_corr <- (WS700_2[[1]]$data[i_period]$WS_0580 - WD_Offset_WS2) %% 360
#
colClasses(WS700_2[[1]]$data)[c("WS_0580","WD_corr")] <- "circ" # needed for plotting and averaging
WS700_2_d1 <- shift_dt(time_offset["WS2"], WS700_2[[1]]$data)
WS700_2_d2 <- shift_dt(time_offset["WS2"], WS700_2[[2]]$data)
WS700_2_d3 <- shift_dt(time_offset["WS2"], WS700_2[[3]]$data)
WS2 <- WS700_2
WS2[[1]]$data <- WS700_2_d1["05.03.2021 11:10 - 26.03.2021 09:10"]
WS2[[2]]$data <- WS700_2_d2["05.03.2021 11:10 - 26.03.2021 09:10"]
WS2[[3]]$data <- WS700_2_d3["05.03.2021 - 27.03.2021"]


#################
### save data ###
#################

saveRDS(WS2,file=paste0(PathRSaves,"/WS700.rds"))

