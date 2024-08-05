
###############################################################
###############################################################
#####                                                     #####
#####    Data treatment MFC and KELLER pressure sensor    #####
#####                                                     #####
###############################################################
###############################################################

# Author: Marcel Bühler
# Date: August 5, 2024
# Contact: mb@bce.au.dk or Christoph Häni christoph.haeni@bfh.ch
# Description: This script reads in the mass flow controller (MFC) and pressure sensor data and makes it ready for further use.
#
# Note: This code was written by Marcel Bühler and is intended to follow the publication 'Applicability of the inverse dispersion method to measure emissions from animal housings' in AMT. 
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


########################
### Define Campaigns ###
########################

IC1 <- "05.03.2021 to 10.03.2021 18:00"
MC <- "18.03.2021 11:00 - 21.03.2021 14:00"
IC2	<- "21.03.2021 14:00 to "


######################################
######################################
#####                            #####
#####    Mass flow controller    #####
#####                            #####
######################################
######################################

MFC_dir <- dir(file.path(PathData,"MFC"), full.names = TRUE)
MFC_list <- lapply(MFC_dir[1:5], fread)
MFC_data <- rbindlist(MFC_list)

names(MFC_data) <- c("Date","Counter_value","Temp_MFC","Flow_MFC_lnmin","Setpoint_MFC_lnmin")
MFC_data[,Time_new := parse_date_time3(Date,tz="Etc/GMT-1")]
## calculate Massflow in kg/h
MFC_data[,Q_MFC := Flow_MFC_lnmin * 60 * 0.717 / 1000] # (Qmfc [nl/min] * 60 [min/h] * 0.717 [g/nl]) / 1000 [g/kg] --> [nl min g kg] / [min h nl g]
d_t <- MFC_data[, as.numeric(diff(Time_new)) / 2]
md_t <- median(d_t)
d_t[d_t > md_t*2] <- md_t
MFC1 <- as.ibts(MFC_data[,c(3:5,7)], st = MFC_data[, Time_new] - c(md_t, d_t), 
  et = MFC_data[, Time_new] + c(d_t, md_t))
MFC_1min <- pool(MFC1,granularity="1mins", st.to="09.03.2021 15:42")
## fill values with 0 at times MFC was not running
MFC_dt <- cbind(as.data.table(MFC_1min),st=st(MFC_1min))
MK_vect <- data.table(st=parse_date_time3("18.03.2021 11:00",tz="Etc/GMT-1") + c(0:6000)*60,Q_MFCex=0)
MFC_dt2 <- merge(MFC_dt,MK_vect,all=TRUE)
MFC_dt2[,et := st + 60]
MFC_dt2[!is.na(Q_MFC),Q_MFCex := Q_MFC]
# MFC_dt2[st>"2021-03-20 03:53" & st< "2021-03-20 04:05",.(st,Q_MFC)]
MFC_dt2[st> "2021-03-20 03:53" & st< "2021-03-20 04:05",Q_MFCex := 6.0225] # Flow was not logged
MFC <- as.ibts(MFC_dt2,st='st',et='et')


#################
### save data ###
#################

saveRDS(MFC,file.path(PathRSaves,"MFC.rds"))


####################################################################################################
### Check, if the released CH4 fits with the amout of gas originally inside the cyclidner bundle ###
####################################################################################################

# pool to 1 min
plot(MFC["19.03.2021 to ","Flow_MFC_lnmin"]) # check if timeperiod is correct
plot(MFC["22.03.2021 to ","Flow_MFC_lnmin"]) # check if timeperiod is correct
abline(v=parse_date_time3("22.03.2021	12:30",tz="Etc/GMT-1"),col="red")
colSums(MFC["19.03.2021 to 22.03.2021 12:30","Flow_MFC_lnmin"],na.rm=TRUE) # select the correct period and sum up the Flow
# --> 141456 nl ---> 141 m^3 of methane
# ---> In the bundle was supposted to be 151 m^3 methane inside by 15°C.
151*273.15/288.15
# --> 143.14 m3 # with norm conditions. 1 bar remaining inside are about 1 m3
# --> 142.5 m3 should be inside and 141 m3 was released. That is a difference of 1% that is fine!
colSums(MFC["22.03.2021 12:30 to ","Flow_MFC_lnmin"],na.rm=TRUE) # select the correct period and sum up the Flow
# --> 12610 nl ---> 12.6 m^3 of methane
# ---> In the bundle was supposted to be 13 m^3 methane inside by 15°C (I actually don't know. just divided the 151 by 12).
12.6*273.15/288.15
# --> 12 m3 # with norm conditions. 1 bar remaining inside are about 1 m3
# --> 11 m3 should be inside and 13 m3 was released. That is a bit too much difference. But I don't know the exact number of what was inside!
plot(MFC[" to 19.03.2021","Flow_MFC_lnmin"]) # check if timeperiod is correct
colSums(MFC[" to 19.03.2021","Flow_MFC_lnmin"],na.rm=TRUE) # select the correct period and sum up the Flow
# 11456.5 nl --> 11.4 m3. That is fine as we expect around 11 m3


#################################
#################################
#####                       #####
#####    Pressure sensor    #####
#####                       #####
#################################
#################################


CC_dir <- dir(file.path(PathData,"Pressuresensor"), full.names = TRUE)
CC_list <- lapply(CC_dir[1:5], fread,skip=5,sep=";",fill=TRUE,header=TRUE,na.strings=c(",,","#NV"))
CC_data <- rbindlist(CC_list)
CC_data[,V5 := NULL]
names(CC_data) <- c("Time_P1","Press_CC","Time_TOB1","Temp_CC")
CC_data[,Time := parse_date_time3(Time_P1,tz="Etc/GMT-1")]
CC_data[,dP_CC := c(NA, diff(Press_CC)*-1)] # calculate dP
# CC_data[,dT_JT := dP_CC * 0.41] # Joule-Thomson effect der Temperatur
# CC_data[dT_JT < 0 | dT_JT > 10,dT_JT := 0] # Joule-Thomson effect der Temperatur
CC_data[,V_Cylinder := 0.05] # add Volumne of 1 cylinder
d_t <- CC_data[, as.numeric(diff(Time)) / 2]
md_t <- median(d_t)
d_t[d_t > md_t*2] <- md_t
CC_Sensor <- as.ibts(CC_data[,c(-1,-3,-5)], st = CC_data[, Time] - c(md_t, d_t), 
  et = CC_data[, Time] + c(d_t, md_t))

## changing colClasses from avg to sum
colClasses(CC_Sensor)["dP_CC"] <- "sum"
# colClasses(CC_Sensor)["dT_JT"] <- "sum"
CC_Sensor["19.03.2021	- 22.03.2021 12:30","V_Cylinder"] <- 0.6


#################
### save data ###
#################

saveRDS(CC_Sensor,file.path(PfadRSaves,"CC_Sensor_raw.rds"))


#######################################
### check if things are alright etc ###
#######################################

## pool to 10min
CC_Sensor_10min <- pool(CC_Sensor,granularity="10mins",st="06.03.2021 08:40")

CC_Sensor_10min[,"Q_CC"] <- NA_real_
a <- 225*10^-3
b <- 42.8*10^-6
R <- 8.31446261815324
M <- 16.043*10^-3

for(i in 1:length(CC_Sensor_10min$V_Cylinder)){
	f <- function(m) { (((m/M) * R * (CC_Sensor_10min[i]$Temp_CC+273.15)) / (CC_Sensor_10min[i]$V_Cylinder - m/M * b) - 
		((m/M)^2 * a )/ (CC_Sensor_10min[i]$V_Cylinder)^2 - CC_Sensor_10min[i]$dP_CC*10^5*6)^2 }
	out <- optimize(f,interval=c(-100,10000))
	CC_Sensor_10min[i, "Q_CC"] <- as.numeric(out[1])
}


#####################
### some plotting ###
#####################

WS700_2 <- readRDS(file.path(PathRSaves,"WS700.rds"))

plot(MFC[MC,"Q_MFC"],ylim=c(-10,10))
lines(CC_Sensor_10min[,"Q_CC"],col="red")

plot(MFC["22.03.2021 to ","Q_MFC"],ylim=c(-10,10))
lines(CC_Sensor_10min[,"Q_CC"],col="red")

plot(MFC["19.03.2021 to ","Temp_MFC"])
lines(CC_Sensor_10min[,"Temp_CC"],col="blue")
lines(WS700_2[[1]]$data[,"WS_0160"],col="magenta")
# lines(CC_Sensor_10min[,"dT_JT"],col="red")
lines(CC_Sensor_10min[,"T_sum"],col="red")

plot(CC_Sensor_10min["19.03.2021 to ","dT_JT"],col="red",ylim=c(-5,5))


CC_Sensor["19.03.2021 12:00"]
# 168.68
CC_Sensor["19.03.2021 12:30"]
# 163.12
## 5.56 bar
5.56*10^5 * 0.6 * M / R / 278 * 2

## pool to 1min
CC_Sensor_1min <- pool(CC_Sensor,granularity="1mins",st="06.03.2021 08:40")

# CC_Sensor_1min[,"T_sum"] <- NA_real_
# CC_Sensor_1min[which(CC_Sensor_10min["19.03.2021 08:00 - 22.03.2021 12:30","dT_JT"] != "NA"),"T_sum"] <- cumsum(CC_Sensor_10min[which(CC_Sensor_10min["19.03.2021 08:00 - 22.03.2021 12:30","dT_JT"] != "NA"),"dT_JT"])

CC_Sensor_1min[,"Q_CC"] <- NA_real_

a <- 225*10^-3
b <- 42.8*10^-6
R <- 8.31446261815324
M <- 16.043*10^-3

for(i in 1:length(CC_Sensor_1min$V_Cylinder)){
	f <- function(m) { (((m/M) * R * (CC_Sensor_1min[i]$Temp_CC+273.15)) / (CC_Sensor_1min[i]$V_Cylinder - m/M * b) - 
		((m/M)^2 * a )/ (CC_Sensor_1min[i]$V_Cylinder)^2 - CC_Sensor_1min[i]$dP_CC*10^5*60)^2 } # 1 min
	out <- optimize(f,interval=c(-100,10000))
	CC_Sensor_1min[i, "Q_CC"] <- as.numeric(out[1])
}

plot(MFC["19.03.2021 to ","Q_MFC"],ylim=c(-10,100))
plot(MFC["22.03.2021 to ","Q_MFC"],ylim=c(-10,10))
lines(CC_Sensor_1min[,"Q_CC"],col="red")


## The temperature measurement of the pressure sensor is not representativ for the cylinder temperature or the point where the gas was expanded
## Here I calculate the expected temperature to match with the flow of the MFC

## ad MFC to CC_Sensor
CC_MFC <- merge(CC_Sensor_1min,MFC[,"Q_MFC"])

## use optimise as above but with the temp as unknown
CC_MFC[,"Temp_CC_exp"] <- NA_real_

a <- 225*10^-3
b <- 42.8*10^-6
R <- 8.31446261815324
M <- 16.043*10^-3

for(i in 1:length(CC_MFC$V_Cylinder)){
	f <- function(T) { (((CC_MFC[i]$Q_MFC/M) * R * T) / (CC_MFC[i]$V_Cylinder - CC_MFC[i]$Q_MFC/M * b) - 
		((CC_MFC[i]$Q_MFC/M)^2 * a )/ (CC_MFC[i]$V_Cylinder)^2 - CC_MFC[i]$dP_CC*10^5*60)^2 } # 1 min
	out <- optimize(f,interval=c(-100,10000))
	CC_MFC[i, "Temp_CC_exp"] <- as.numeric(out[1])
}


## Temp plotten
plot(CC_MFC[MC,"Temp_CC"],col="blue",ylim=c(-150,10))
lines(CC_MFC[MC,"Temp_CC_exp"]-273.15,col="red")

plot(CC_MFC[IC2,"Temp_CC"],col="blue",ylim=c(-150,10))
lines(CC_MFC[IC2,"Temp_CC_exp"]-273.15,col="red")

## such low expected temperatures are not realistic. The pressure must also be affected

## anyhow, recalculate flow with the expected temperature just to see if it worked
CC_MFC[,"Q_CC_exp"] <- NA_real_

a <- 225*10^-3
b <- 42.8*10^-6
R <- 8.31446261815324
M <- 16.043*10^-3

for(i in 1:length(CC_MFC$V_Cylinder)){
	f <- function(m) { (((m/M) * R * (CC_MFC[i]$Temp_CC_exp)) / (CC_MFC[i]$V_Cylinder - m/M * b) - 
		((m/M)^2 * a )/ (CC_MFC[i]$V_Cylinder)^2 - CC_MFC[i]$dP_CC*10^5*60)^2 } # 1 min
	out <- optimize(f,interval=c(-100,10000))
	CC_MFC[i, "Q_CC_exp"] <- as.numeric(out[1])
}

plot(MFC["19.03.2021 to ","Q_MFC"],ylim=c(-10,10))
plot(MFC_1min["22.03.2021 to ","Q_MFC"],ylim=c(-10,10))
lines(CC_MFC[,"Q_CC_exp"],col="red")
lines(CC_MFC[,"Q_CC"],col="blue")

## there are some really high values that are not correct
iQ <- which(CC_MFC[,"Q_CC_exp"] > 20)

plot(CC_MFC[MC,"Q_CC"],col="blue",ylim=c(-10,10),type="p",pch=19)
lines(CC_MFC[iQ,"Q_CC"],col="red",type="p",pch=19)

## check what temperature there was
plot(CC_MFC[MC,"Temp_CC"],col="blue",ylim=c(-10,10),type="p",pch=19)
lines(CC_MFC[iQ,"Temp_CC"],col="red",type="p",pch=19)
## all increasing temperature. Why does this cause an error?

## check pressure
plot(CC_MFC[MC,"Press_CC"],col="blue",type="p",pch=19)
lines(CC_MFC[iQ,"Press_CC"],col="red",type="p",pch=19)
# really low pressures (probably wrong dP but others are not special)
# plot(CC_MFC[MC,"dP_CC"],col="blue",type="p",pch=19)
plot(CC_MFC[MC,"dP_CC"],col="blue",type="p",pch=19,ylim=c(-0.5,0.5))
lines(CC_MFC[iQ,"dP_CC"],col="red",type="p",pch=19)
abline(h=0.15)
# large dPs, negative dPs, and dPs that are below 0.15 bar



plot(MFC_1min["19.03.2021 to ","Temp_MFC"])
lines(CC_Sensor_10min[,"Temp_CC"],col="blue")
lines(WS700_2[[1]]$data[,"WS_0160"],col="magenta")
# lines(CC_Sensor_10min[,"dT_JT"],col="red")
lines(CC_Sensor_10min[,"T_sum"],col="red")

plot(CC_Sensor_10min["19.03.2021 to ","dT_JT"],col="red",ylim=c(-5,5))


CC_Sensor_10min["19.03.2021 12:00 to ","dP_CC"] * 10E5

rowSums(CC_Sensor_10min["19.03.2021 to 22.03.2021 12:30","dT_JT"] )
colSums(CC_Sensor_10min["19.03.2021 to 22.03.2021 12:30"]$dT_JT)

CC_data[,T_sum := cumsum(dT_JT,na.rm=TRUE)] 

