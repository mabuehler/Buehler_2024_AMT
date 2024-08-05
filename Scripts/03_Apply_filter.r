
#####################################################
#####################################################
#####                                           #####
#####    Quality filter of the emission data    #####
#####                                           #####
#####################################################
#####################################################

# Author: Marcel Bühler
# Date: August 5, 2024
# Contact: mb@bce.au.dk or Christoph Häni christoph.haeni@bfh.ch
# Description: This script applies the quality filtering and makes the data ready for further use.
#
# Note: This code was written by Marcel Bühler and is intended to follow the publication 'Applicability of the inverse dispersion method to measure emissions from animal housings' in AMT. 
# Please feel free to use and modify it, but attribution is appreciated.


#################
### Libraries ###
#################

library(ibts)
library(RgoogleMaps)
library(bLSmodelR)
library(RColorBrewer)
library(ggplot2)


#############
### Paths ###
#############

PathData <- "Path to /data"		
PathRSaves <- "Path to /RSaves"
PathFigures <- 'Path to /Figures'


#################
### Functions ###
#################

source("https://raw.githubusercontent.com/hafl-gel/gel-scripts/main/wgs84-ch1903.r")

lines_sec2xy <- function(xyMK,sensor,node=1,wd,col="lightblue",lwd=2,...){
	GF <- xyMK[xyMK[,1] %in% sensor,]
	sens <- as.numeric(GF[GF[,3] == node,4:5])
	b <- tan((90 - wd)/180*pi)
	x <- if(wd <= 180) 600 else -600
	y <- sens[2] - (sens[1] - x)*b
	lines(c(sens[1],x),c(sens[2],y),col=col,lwd=lwd,...)
}


###################
### Definitions ###
###################

	Col_CH4 <- "#1e88e5"
	# display.brewer.all()
	CH4Cols <- brewer.pal(5,"Dark2")
	names(CH4Cols) <- paste0("GF",c(16:18,25:26))

	SonicCols <- brewer.pal(4,"BrBG")
	names(SonicCols) <- paste0("Sonic",c("A",2,"B","C"))

	GFnames <- c("GF17","GF18","GF16","GF25","GF26")
	Sonicnames <- paste0("Sonic",c("A",2,"B","C"))

	mgs_to_kgd <- (24*3600)/1E6
	mgs_to_kgh <- (3600)/1E6
	mgm2s_to_gm2h <- 3600/1E3


#################
### Campaigns ###
#################

	IC1 <- "05.03.2021 to 10.03.2021 18:00"
	MC <- "18.03.2021 11:00 - 21.03.2021 14:00"
	IC2	<- "21.03.2021 14:00 to "


#################
### load data ###
#################

	# Geometry
	load(file.path(PathRSaves,"Geometry.RData"))
	# Google map
	STO_Map <- ReadMapTile(file.path(PathFigures,"STO_GoogleMaps.png"))
	# Emission results
	Emiss_Result_raw <- readRDS(file = file.path(PathRSaves, "Emiss_orig_P23_10min.rds")) # <-------------------- the one to use

	## not needed. Only necessary if you wanna reproduce a plot in the initial submission
	# Emiss_Result_orig_wo <- readRDS(file = file.path(PfadRSaves, "Emiss_orig_wo_10min.rds"))
	# Emiss_Result_orig_P23 <- readRDS(file = file.path(PfadRSaves, "Emiss_orig_P23_10min.rds"))
	# Emiss_Result_orig_P6 <- readRDS(file = file.path(PfadRSaves, "Emiss_orig_P6_10min.rds"))
	# Emiss_Result_orig_P4 <- readRDS(file = file.path(PfadRSaves, "Emiss_orig_P4_10min.rds"))
	# Emiss_Result_WDvar_P23 <- readRDS(file = file.path(PfadRSaves, "STO_Emiss_WDvar_P23_10min.rds"))
	# Emiss_Result_WDvar_P23 <- readRDS(file = file.path(PfadRSaves, "STO_Emiss_WDvar_P6_10min.rds"))


################################
### Create complete data set ###
################################

# create complete time data set because of plotting (18.03.2021 12:50 - 22.03.2021 14:20)
time_vect <- data.table(st=seq(0,585)*600 + parse_date_time3("18.03.2021 12:50",tz="Etc/GMT-1"))
Emiss_Result <- merge(time_vect,Emiss_Result_raw,by="st",all=TRUE)
# MFC für Darstellung auf 0 setzen vor und nach Release
# Emiss_Result[Campaign=="MC" & is.na(Q_MFC),]

setkey(Emiss_Result,st)


##########################
### create XY Geometry ###
##########################

Sensors_MC_xy <- ch_to_map(STO_Map,Sensors_MC)
Sensors_IC2_xy <- ch_to_map(STO_Map,Sensors_IC2)
Source_xy <- ch_to_map(STO_Map,Source)


#########################
### plot and overview ###
#########################

plot(Q_GF16 ~ st, Emiss_Result[Campaign %in% c("MC")],type="l",col=CH4Cols["GF16"],ylim=c(-10,10))
lines(Q_GF17 ~ st, Emiss_Result[Campaign %in% c("MC")],col=CH4Cols["GF17"])
lines(Q_GF18 ~ st, Emiss_Result[Campaign %in% c("MC")],col=CH4Cols["GF18"])
lines(Q_GF25 ~ st, Emiss_Result[Campaign %in% c("MC")],col=CH4Cols["GF25"])
lines(Q_GF26 ~ st, Emiss_Result[Campaign %in% c("MC")],col=CH4Cols["GF26"])
lines(Q_MFC ~ st, Emiss_Result[Campaign %in% c("MC")],col="black",lwd=2)


##################################
### only use data with NE wind ###
##################################

Result_raw <- Emiss_Result[Wind_dir == "NE" | is.na(Wind_dir)]


################################################
### add missing times caused by power outage ###
################################################

# the difference to above is, that the name of the Source and the Sonic should be filled in as well
t_step <- as.numeric(median(Result_raw[,et-st],na.rm=TRUE))
st_fill <- parse_date_time3("2021-03-19 09:00",tz='Etc/GMT-1')
st_dt <- data.table(st=st_fill + seq(0,length.out=250/t_step,by=t_step) * 60, et=st_fill + seq(0,length.out=250/t_step,by=t_step) * 60 + t_step*60)

Result_ls <- lapply(Result_raw[,unique(Sonic)], function(x) {
	# browser()
	out <- merge(Result_raw[Sonic==x],st_dt,by=c('st','et'),all=TRUE)
	out[,Sonic := x]
	out[st %in% st_dt[1],Campaign := 'MC']
	out[,Source := 'Schopf']
})

Result <- rbindlist(Result_ls)

setkey(Result,st)


#######################################
### Wind direction filtering for MC ###
#######################################

## This is testing and might not be used

plot(Source,Sensors_MC)
# Define southwest and SE points of the source (Schopf)
# Schopf NW corner
ENW <- as.numeric(Source[Source[,3] == max(Source[,3]),2:3])
# Schopf NE corner
ENE <- as.numeric(Source[Source[,2] == max(Source[,2]),2:3])
# Schopf SW corner
ESW <- as.numeric(Source[Source[,2] == min(Source[,2]),2:3])
# Schopf SE corner
ESE <- as.numeric(Source[Source[,3] == min(Source[,3]),2:3])
## Determine middle of the GasFinder path
P16M <- as.numeric(colMeans(Sensors_MC[Sensors_MC[,1] == "GF16" ,c(4,5)]))
P17M <- as.numeric(colMeans(Sensors_MC[Sensors_MC[,1] == "GF17" ,c(4,5)]))
P18M <- as.numeric(colMeans(Sensors_MC[Sensors_MC[,1] == "GF18" ,c(4,5)]))
P25M <- as.numeric(colMeans(Sensors_MC[Sensors_MC[,1] == "GF25" ,c(4,5)]))

## Generate dummy Sensor in the middle of the path for checking angles
PMitte <- genSensors(data.frame(Name=c("GF16","GF17","GF18","GF25")
	,x = c(P16M[1],P17M[1],P18M[1],P25M[1])
	,y = c(P16M[2],P17M[2],P18M[2],P25M[2])
	,z=rep(1,4),rep(NA,4),rep(1,4)
	))

Mitte_xy <- ch_to_map(STO_Map,PMitte)

### Determine angles for variante 'edges' (I usually did that manually with the map) for the MC
iWD16 <- c(90 - round(atan(abs(Sensors_MC[Sensors_MC[,1] == "GF16" & Sensors_MC[,3] == 1,5]-ENW[2])/
		abs(Sensors_MC[Sensors_MC[,1] == "GF16" & Sensors_MC[,3] == 1,4]-ENW[1]))*180/pi,0),
	0 + round(atan(abs(Sensors_MC[Sensors_MC[,1] == "GF16" & Sensors_MC[,3] == 2,4]-ESE[1])/
		abs(Sensors_MC[Sensors_MC[,1] == "GF16" & Sensors_MC[,3] == 2,5]-ESE[2]))*180/pi,0))
iWD17 <- c(360 - round(atan(abs(Sensors_MC[Sensors_MC[,1] == "GF17" & Sensors_MC[,3] == 1,4]-ESW[1])/
		abs(Sensors_MC[Sensors_MC[,1] == "GF17" & Sensors_MC[,3] == 1,5]-ESW[2]))*180/pi,0),
	90 + round(atan(abs(Sensors_MC[Sensors_MC[,1] == "GF17" & Sensors_MC[,3] == 2,5]-ESE[2])/
		abs(Sensors_MC[Sensors_MC[,1] == "GF17" & Sensors_MC[,3] == 2,4]-ESE[1]))*180/pi,0))
iWD18 <- c(90 - round(atan(abs(Sensors_MC[Sensors_MC[,1] == "GF18" & Sensors_MC[,3] == 1,5]-ESW[2])/
		abs(Sensors_MC[Sensors_MC[,1] == "GF18" & Sensors_MC[,3] == 1,4]-ESW[1]))*180/pi,0),
	0 + round(atan(abs(Sensors_MC[Sensors_MC[,1] == "GF18" & Sensors_MC[,3] == 2,4]-ESE[1])/
		abs(Sensors_MC[Sensors_MC[,1] == "GF18" & Sensors_MC[,3] == 2,5]-ESE[2]))*180/pi,0))
iWD25 <- c(90 - round(atan(abs(Sensors_MC[Sensors_MC[,1] == "GF25" & Sensors_MC[,3] == 1,5]-ENW[2])/
		abs(Sensors_MC[Sensors_MC[,1] == "GF25" & Sensors_MC[,3] == 1,4]-ENW[1]))*180/pi,0),
	0 + round(atan(abs(Sensors_MC[Sensors_MC[,1] == "GF25" & Sensors_MC[,3] == 2,4]-ESE[1])/
		abs(Sensors_MC[Sensors_MC[,1] == "GF25" & Sensors_MC[,3] == 2,5]-ESE[2]))*180/pi,0))

### Determine angles for the variante 'middle' for the MC
iWD16M <- c(0 + round(atan(abs(P16M[1]-ENW[1])/abs(P16M[2]-ENW[2]))*180/pi,0),
	90 - round(atan(abs(P16M[2]-ESE[2])/abs(P16M[1]-ESE[1]))*180/pi,0))
iWD17M <- c(0 + round(atan(abs(P17M[1]-ENW[1])/abs(P17M[2]-ENW[2]))*180/pi,0),
	90 - round(atan(abs(P17M[2]-ESE[2])/abs(P17M[1]-ESE[1]))*180/pi,0))
iWD18M <- c(0 + round(atan(abs(P18M[1]-ENW[1])/abs(P18M[2]-ENW[2]))*180/pi,0),
	90 - round(atan(abs(P18M[2]-ESE[2])/abs(P18M[1]-ESE[1]))*180/pi,0))
iWD25M <- c(0 + round(atan(abs(P25M[1]-ENW[1])/abs(P25M[2]-ENW[2]))*180/pi,0),
	90 - round(atan(abs(P25M[2]-ESE[2])/abs(P25M[1]-ESE[1]))*180/pi,0))


########################################
### wind direction filtering for IC2 ###
########################################

## This is testing and might not be used

plot(Source,Sensors_IC2)
## Only determine for GF18 and then apply to all others, as GF18 is in the middle and differences are probably very small

# Determine SW and SE corners of the source
# Schopf corner NW
ENW <- as.numeric(Source[Source[,3] == max(Source[,3]),2:3])
# Schopf corner NE
ENE <- as.numeric(Source[Source[,2] == max(Source[,2]),2:3])
# Schopf corner SW
ESW <- as.numeric(Source[Source[,2] == min(Source[,2]),2:3])
# Schopf corner SE
ESE <- as.numeric(Source[Source[,3] == min(Source[,3]),2:3])

## Determine middle of GF18 path
PMIC2 <- as.numeric(colMeans(Sensors_IC2[Sensors_IC2[,1] == "GF18" ,c(4,5)]))
## create dummy Sensor in the middel of the path
PMitteIC2 <- genSensors(data.frame(Name=c("GFM")
	,x = PMIC2[1]
	,y = PMIC2[2]
	,z=1,NA,1
	))

MitteIC2_xy <- ch_to_map(STO_Map,PMitteIC2)

### Determine angles for variante 'edges' (I usually did that manually with the map) for the IC2
iWDIC2 <- c(360 - round(atan(abs(Sensors_IC2[Sensors_IC2[,1] == "GF18" & Sensors_IC2[,3] == 2,4]-ESW[1])/
		abs(Sensors_IC2[Sensors_IC2[,1] == "GF18" & Sensors_IC2[,3] == 2,5]-ESW[2]))*180/pi,0),
	90 + round(atan(abs(Sensors_IC2[Sensors_IC2[,1] == "GF18" & Sensors_IC2[,3] == 1,5]-ESE[2])/
		abs(Sensors_IC2[Sensors_IC2[,1] == "GF18" & Sensors_IC2[,3] == 1,4]-ESE[1]))*180/pi,0))

### Determine angles for the variante 'middle' for the IC2
iWDMIC2 <- c(0 + round(atan(abs(PMIC2[1]-ENW[1])/abs(PMIC2[2]-ENW[2]))*180/pi,0),
	90 - round(atan(abs(PMIC2[2]-ESE[2])/abs(PMIC2[1]-ESE[1]))*180/pi,0))


####################################
### check/validate angles on map ###
####################################

graphics.off()

## GF17
# iWD17neu <- c(22,68)
PlotOnStaticMap(STO_Map)
plot(Sensors_MC_xy[Sensors_MC_xy[,1] %in% c("GF17"),],sensor.text.args=list(labels="",cex=2),lines.args=list(lwd=2,col="orange",lty=3),points.args = list(pch = 20, cex = 2, col ="orange"),add=TRUE)
lines_sec2xy(Sensors_MC_xy,"GF17",1,iWD17[1],lty=2)
lines_sec2xy(Sensors_MC_xy,"GF17",2,iWD17[2],lty=2)
lines_sec2xy(Sensors_MC_xy,"GF17",2,iWD17[1])
lines_sec2xy(Sensors_MC_xy,"GF17",1,iWD17[2])
lines_sec2xy(Mitte_xy,"GF17",1,iWD17M[1],col="pink",lty=3)
lines_sec2xy(Mitte_xy,"GF17",1,iWD17M[2],col="pink",lty=3)
lines_sec2xy(Sensors_MC_xy,"GF17",1,iWD17M[2],col="pink")
lines_sec2xy(Sensors_MC_xy,"GF17",1,iWD17M[1],lty=2,col="pink")
lines_sec2xy(Sensors_MC_xy,"GF17",2,iWD17M[1],col="pink")
lines_sec2xy(Sensors_MC_xy,"GF17",2,iWD17M[2],lty=2,col="pink")

## GF18
# iWD18neu <- c(22,68)
PlotOnStaticMap(STO_Map)
plot(Sensors_MC_xy[Sensors_MC_xy[,1] %in% c("GF18"),],sensor.text.args=list(labels="",cex=2),lines.args=list(lwd=2,col="orange",lty=3),points.args = list(pch = 20, cex = 2, col ="orange"),add=TRUE)
lines_sec2xy(Sensors_MC_xy,"GF18",1,iWD18[1],lty=2)
lines_sec2xy(Sensors_MC_xy,"GF18",2,iWD18[2],lty=2)
lines_sec2xy(Sensors_MC_xy,"GF18",2,iWD18[1])
lines_sec2xy(Sensors_MC_xy,"GF18",1,iWD18[2])
lines_sec2xy(Mitte_xy,"GF18",1,iWD18M[1],col="pink",lty=3)
lines_sec2xy(Mitte_xy,"GF18",1,iWD18M[2],col="pink",lty=3)
lines_sec2xy(Sensors_MC_xy,"GF18",1,iWD18M[2],col="pink")
lines_sec2xy(Sensors_MC_xy,"GF18",1,iWD18M[1],lty=2,col="pink")
lines_sec2xy(Sensors_MC_xy,"GF18",2,iWD18M[1],col="pink")
lines_sec2xy(Sensors_MC_xy,"GF18",2,iWD18M[2],lty=2,col="pink")

## GF16
# iWD16neu <- c(22,68)
PlotOnStaticMap(STO_Map)
plot(Sensors_MC_xy[Sensors_MC_xy[,1] %in% c("GF16"),],sensor.text.args=list(labels="",cex=2),lines.args=list(lwd=2,col="orange",lty=3),points.args = list(pch = 20, cex = 2, col ="orange"),add=TRUE)
lines_sec2xy(Sensors_MC_xy,"GF16",1,iWD16[1],lty=2)
lines_sec2xy(Sensors_MC_xy,"GF16",2,iWD16[2],lty=2)
lines_sec2xy(Sensors_MC_xy,"GF16",2,iWD16[1])
lines_sec2xy(Sensors_MC_xy,"GF16",1,iWD16[2])
lines_sec2xy(Mitte_xy,"GF16",1,iWD16M[1],col="pink",lty=3)
lines_sec2xy(Mitte_xy,"GF16",1,iWD16M[2],col="pink",lty=3)
lines_sec2xy(Sensors_MC_xy,"GF16",1,iWD16M[2],col="pink")
lines_sec2xy(Sensors_MC_xy,"GF16",1,iWD16M[1],lty=2,col="pink")
lines_sec2xy(Sensors_MC_xy,"GF16",2,iWD16M[1],col="pink")
lines_sec2xy(Sensors_MC_xy,"GF16",2,iWD16M[2],lty=2,col="pink")

## GF25
# iWD25neu <- c(22,68)
PlotOnStaticMap(STO_Map)
plot(Sensors_MC_xy[Sensors_MC_xy[,1] %in% c("GF25"),],sensor.text.args=list(labels="",cex=2),lines.args=list(lwd=2,col="orange",lty=3),points.args = list(pch = 20, cex = 2, col ="orange"),add=TRUE)
lines_sec2xy(Sensors_MC_xy,"GF25",1,iWD25[1],lty=2)
lines_sec2xy(Sensors_MC_xy,"GF25",2,iWD25[2],lty=2)
lines_sec2xy(Sensors_MC_xy,"GF25",2,iWD25[1])
lines_sec2xy(Sensors_MC_xy,"GF25",1,iWD25[2])
lines_sec2xy(Mitte_xy,"GF25",1,iWD25M[1],col="pink",lty=3)
lines_sec2xy(Mitte_xy,"GF25",1,iWD25M[2],col="pink",lty=3)
lines_sec2xy(Sensors_MC_xy,"GF25",1,iWD25M[2],col="pink")
lines_sec2xy(Sensors_MC_xy,"GF25",1,iWD25M[1],lty=2,col="pink")
lines_sec2xy(Sensors_MC_xy,"GF25",2,iWD25M[1],col="pink")
lines_sec2xy(Sensors_MC_xy,"GF25",2,iWD25M[2],lty=2,col="pink")


########### IC2
PlotOnStaticMap(STO_Map)
plot(Sensors_IC2_xy[Sensors_IC2_xy[,1] %in% c("GF18","GF26"),],sensor.text.args=list(labels="",cex=2),lines.args=list(lwd=2,col="orange",lty=3),points.args = list(pch = 20, cex = 2, col ="orange"),add=TRUE)
lines_sec2xy(Sensors_IC2_xy,"GF18",1,iWDIC2[2],lty=2)
lines_sec2xy(Sensors_IC2_xy,"GF18",2,iWDIC2[1],lty=2)
lines_sec2xy(Sensors_IC2_xy,"GF18",1,iWDIC2[1])
lines_sec2xy(Sensors_IC2_xy,"GF18",2,iWDIC2[2])
lines_sec2xy(MitteIC2_xy,"GFM",1,iWDMIC2[1],col="pink",lty=3)
lines_sec2xy(MitteIC2_xy,"GFM",1,iWDMIC2[2],col="pink",lty=3)
lines_sec2xy(Sensors_IC2_xy,"GF18",1,iWDMIC2[2],lty=2,col="pink")
lines_sec2xy(Sensors_IC2_xy,"GF18",1,iWDMIC2[1],col="pink")
lines_sec2xy(Sensors_IC2_xy,"GF18",2,iWDMIC2[1],lty=2,col="pink")
lines_sec2xy(Sensors_IC2_xy,"GF18",2,iWDMIC2[2],col="pink")

graphics.off()


## filter according to map
iWD <- list(
	iWD_GF16 = iWD16
	,iWD_GF17 = iWD17
	,iWD_GF18 = iWD18
	,iWD_GF25 = iWD25
)

## according to middle point of the path (ARA Method)
iWDM <- list(
	iWD_GF16 = iWD16M
	,iWD_GF17 = iWD17M
	,iWD_GF18 = iWD18M
	,iWD_GF25 = iWD25M
)

##################################
### Emiss (unfiltered) vs. WD  ###
##################################

graphics.off()

par(mfrow=c(2,5),mar=c(4,4,0,0))
plot(Q_GF17 ~ I((WD_WS2 + 180) %% 360 - 180),Result[Campaign == "MC"],xlab="WD WS", col= CH4Cols["GF17"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD17 + 180) %% 360 - 180)
for(i in names(SonicCols)){
plot(Q_GF17 ~ I((WD_WS2 + 180) %% 360 - 180),Result[Campaign == "MC" & Sonic == i],xlab="WD WS", col= CH4Cols["GF17"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD17 + 180) %% 360 - 180)
}
plot(Q_GF17 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MC"],xlab="WD WS", col= CH4Cols["GF17"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD17 + 180) %% 360 - 180)
for(i in names(SonicCols)){
plot(Q_GF17 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MC" & Sonic == i],xlab=paste0("WD ", i), col= CH4Cols["GF17"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD17 + 180) %% 360 - 180)
}

par(mfrow=c(2,5),mar=c(4,4,0,0))
plot(Q_GF18 ~ I((WD_WS2 + 180) %% 360 - 180),Result[Campaign == "MC"],xlab="WD WS", col= CH4Cols["GF18"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD18 + 180) %% 360 - 180)
for(i in names(SonicCols)){
plot(Q_GF18 ~ I((WD_WS2 + 180) %% 360 - 180),Result[Campaign == "MC" & Sonic == i],xlab="WD WS", col= CH4Cols["GF18"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD18 + 180) %% 360 - 180)
}
plot(Q_GF18 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MC"],xlab="WD WS", col= CH4Cols["GF18"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD18 + 180) %% 360 - 180)
for(i in names(SonicCols)){
plot(Q_GF18 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MC" & Sonic == i],xlab=paste0("WD ", i), col= CH4Cols["GF18"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD18 + 180) %% 360 - 180)
}

par(mfrow=c(2,5),mar=c(4,4,0,0))
plot(Q_GF16 ~ I((WD_WS2 + 180) %% 360 - 180),Result[Campaign == "MC"],xlab="WD WS", col= CH4Cols["GF16"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD16 + 180) %% 360 - 180)
for(i in names(SonicCols)){
plot(Q_GF16 ~ I((WD_WS2 + 180) %% 360 - 180),Result[Campaign == "MC" & Sonic == i],xlab="WD WS", col= CH4Cols["GF16"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD16 + 180) %% 360 - 180)
}
plot(Q_GF16 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MC"],xlab="WD WS", col= CH4Cols["GF16"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD16 + 180) %% 360 - 180)
for(i in names(SonicCols)){
plot(Q_GF16 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MC" & Sonic == i],xlab=paste0("WD ", i), col= CH4Cols["GF16"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD16 + 180) %% 360 - 180)
}

par(mfrow=c(2,5),mar=c(4,4,0,0))
plot(Q_GF25 ~ I((WD_WS2 + 180) %% 360 - 180),Result[Campaign == "MC"],xlab="WD WS", col= CH4Cols["GF25"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD25 + 180) %% 360 - 180)
for(i in names(SonicCols)){
plot(Q_GF25 ~ I((WD_WS2 + 180) %% 360 - 180),Result[Campaign == "MC" & Sonic == i],xlab="WD WS", col= CH4Cols["GF25"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD25 + 180) %% 360 - 180)
}
plot(Q_GF25 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MC"],xlab="WD WS", col= CH4Cols["GF25"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD25 + 180) %% 360 - 180)
for(i in names(SonicCols)){
plot(Q_GF25 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MC" & Sonic == i],xlab=paste0("WD ", i), col= CH4Cols["GF25"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD25 + 180) %% 360 - 180)
}


---> WD filtering seems to work better if the wind direction of the individual sonics are used instead of the weather staion.


### WD filter according to map and corresponding Soinc as reference for the wind direction
# MC
for(x in c("GF16","GF17","GF18","GF25")){
	if(x != "GF17"){
	Result[Campaign == "MC" & (WD < iWD[[paste0("iWD_",x)]][1] | WD > iWD[[paste0("iWD_",x)]][2]),
		grep(paste0("Q_",x),names(Result),value=TRUE)	 := NA_real_ ]
		} else {
	Result[Campaign == "MC" & (WD < iWD[[paste0("iWD_",x)]][1] & WD > iWD[[paste0("iWD_",x)]][2]),
		grep(paste0("Q_",x),names(Result),value=TRUE)	 := NA_real_ ]
		}
}
# IC2
for(x in c("GF16","GF17","GF18","GF25","GF26")){
	Result[Campaign == "IC2" & (WD < iWDIC2[1] - 5 & WD > iWDIC2[2] + 5),
		grep(paste0("Q_",x),names(Result),value=TRUE)	 := NA_real_ ]
}


################################
### bLS model quality filter ###
################################

y_lim <- c(-10,30)
y_lim2 <- c(-10,200)
y_lim3 <- range(Result[Q_MFC > 0,.(Q_GF16,Q_GF17,Q_GF18,Q_GF25)],na.rm=TRUE)

### Ustar
# MC
par(mfrow=c(2,2))
for(i in names(SonicCols)){
plot(Q_GF17 ~ Ustar,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("u* ",i),ylim=y_lim,xlim=c(0,0.5))
points(Q_GF18 ~ Ustar,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ Ustar,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ Ustar,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF17","GF18","GF16","GF25")])
abline(h=0,v=c(0.05,0.1,0.15),col="blue")
}
# IC2
for(i in names(SonicCols)){
plot(Q_GF26 ~ Ustar,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF26"],ylab="Q [kg/h]",xlab=paste0("u* ",i),ylim=y_lim,xlim=c(0,0.5))
points(Q_GF17 ~ Ustar,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF17"])
points(Q_GF18 ~ Ustar,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ Ustar,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ Ustar,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(0.05,0.1,0.15),col="blue")
}
# only when MFC has a flow
for(i in names(SonicCols)){
plot(Q_GF17 ~ Ustar,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("u* ",i),ylim=y_lim,xlim=c(0,0.5))
points(Q_GF18 ~ Ustar,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ Ustar,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ Ustar,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(0.05,0.1,0.15),col="blue")
}
# ---> actually no u* filter necessary. If, then without any issue at 0.1 or 0.15
---> Just take 0.15 and nothing else

### sUu
par(mfrow=c(2,2))
for(i in names(SonicCols)){
plot(Q_GF17 ~ sUu,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("sUu ",i),ylim=y_lim,xlim=c(0,16))
points(Q_GF18 ~ sUu,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sUu,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sUu,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF17","GF18","GF16","GF25")])
abline(h=0,v=c(4,5,6),col="blue")
}
# IC2
for(i in names(SonicCols)){
plot(Q_GF26 ~ sUu,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF26"],ylab="Q [kg/h]",xlab=paste0("sUu",i),ylim=y_lim,xlim=c(0,16))
points(Q_GF17 ~ sUu,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF17"])
points(Q_GF18 ~ sUu,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sUu,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sUu,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(4,5,6),col="blue")
}
# only when MFC has a flow
for(i in names(SonicCols)){
plot(Q_GF17 ~ sUu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("sUu",i),ylim=y_lim,xlim=c(0,16))
points(Q_GF18 ~ sUu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sUu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sUu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(4,5,6),col="blue")
}

# ---> 4 or 6. 4 seems to be OK, as there are only values higher with SonicA and SonicC
---> do not filter

### sVu
par(mfrow=c(2,2))
for(i in names(SonicCols)){
plot(Q_GF17 ~ sVu,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("sVu ",i),ylim=y_lim,xlim=c(0,16))
points(Q_GF18 ~ sVu,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sVu,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sVu,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF17","GF18","GF16","GF25")])
abline(h=0,v=c(4,5,6),col="blue")
}
# IC2
for(i in names(SonicCols)){
plot(Q_GF26 ~ sVu,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF26"],ylab="Q [kg/h]",xlab=paste0("sVu",i),ylim=y_lim,xlim=c(0,16))
points(Q_GF17 ~ sVu,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF17"])
points(Q_GF18 ~ sVu,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sVu,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sVu,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(4,5,6),col="blue")
}
# only when MFC has a flow
for(i in names(SonicCols)){
plot(Q_GF17 ~ sVu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("sVu",i),ylim=y_lim,xlim=c(0,16))
points(Q_GF18 ~ sVu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sVu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sVu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(4,5,6),col="blue")
}

# ---> 4 or 5. 4 seems to be OK, as there are only values higher with SonicA and SonicC
---> do not filter

### sWu
par(mfrow=c(2,2))
for(i in names(SonicCols)){
plot(Q_GF17 ~ sWu,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("sWu ",i),ylim=y_lim,xlim=c(0,10))
points(Q_GF18 ~ sWu,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sWu,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sWu,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF17","GF18","GF16","GF25")])
abline(h=0,v=c(2,3,4),col="blue")
}
# IC2
for(i in names(SonicCols)){
plot(Q_GF26 ~ sWu,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF26"],ylab="Q [kg/h]",xlab=paste0("sWu",i),ylim=y_lim,xlim=c(0,10))
points(Q_GF17 ~ sWu,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF17"])
points(Q_GF18 ~ sWu,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sWu,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sWu,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(2,3,4),col="blue")
}
# only when MFC has a flow
for(i in names(SonicCols)){
plot(Q_GF17 ~ sWu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("sWu",i),ylim=y_lim,xlim=c(0,10))
points(Q_GF18 ~ sWu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sWu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sWu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(2,3,4),col="blue")
}

# ---> 2 oder 4. Better 4 or not filter at all, as with SonicA there seems to many good values over the threshold
---> do not filter

### C0
par(mfrow=c(2,2))
for(i in names(SonicCols)){
plot(Q_GF17 ~ C0,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("C0 ",i),ylim=y_lim,xlim=c(0,50))
points(Q_GF18 ~ C0,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ C0,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ C0,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF17","GF18","GF16","GF25")])
abline(h=0,v=c(7,10,30),col="blue")
}
# IC2
for(i in names(SonicCols)){
plot(Q_GF26 ~ C0,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF26"],ylab="Q [kg/h]",xlab=paste0("C0",i),ylim=y_lim,xlim=c(0,50))
points(Q_GF17 ~ C0,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF17"])
points(Q_GF18 ~ C0,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ C0,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ C0,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(7,10,30),col="blue")
}
# only when MFC has a flow
for(i in names(SonicCols)){
plot(Q_GF17 ~ C0,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("C0",i),ylim=y_lim,xlim=c(0,50))
points(Q_GF18 ~ C0,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ C0,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ C0,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(7,10,30),col="blue")
}

# ---> some values are above 10, thus maybe use 30 even though not plausible. SonicA and SonicC affected
----> do not filter

### L
# MC
for(i in names(SonicCols)){
plot(Q_GF17 ~ I(1/L),Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF17"],xlab=paste0("1/L ",i),ylab="Q [kg/h]",ylim=y_lim,xlim=c(-2,2))
points(Q_GF18 ~ I(1/L),Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ I(1/L),Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ I(1/L),Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF17","GF18","GF16","GF25")])
abline(h=0,v=c(-0.5,0.5),col="blue")
}
# IC2
for(i in names(SonicCols)){
plot(Q_GF26 ~ I(1/L),Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF26"],xlab=paste0("1/L ",i),ylab="Q [kg/h]",ylim=y_lim,xlim=c(-2,2))
points(Q_GF17 ~ I(1/L),Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF17"])
points(Q_GF18 ~ I(1/L),Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ I(1/L),Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ I(1/L),Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(-0.5,0.5),col="blue")
}
# only when MFC has a flow
for(i in names(SonicCols)){
plot(Q_GF17 ~ I(1/L),Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF17"],xlab=paste0("1/L ",i),ylab="Q [kg/h]",ylim=y_lim,xlim=c(-2,2))
points(Q_GF18 ~ I(1/L),Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ I(1/L),Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ I(1/L),Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(-0.5,0.5),col="blue")
}

--> abs(0.5) or could even take abs(0.2)
--> do not filter

### Zo
# MC
par(mfrow=c(2,2))
for(i in names(SonicCols)){
plot(Q_GF17 ~ Zo,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("Zo ",i),ylim=y_lim,xlim=c(0,1))
points(Q_GF18 ~ Zo,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ Zo,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ Zo,Result[Campaign == "MC" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF17","GF18","GF16","GF25")])
abline(h=0,v=c(0.05,0.07,0.1),col="blue")
}
# IC2
for(i in names(SonicCols)){
plot(Q_GF26 ~ Zo,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF26"],ylab="Q [kg/h]",xlab=paste0("Zo",i),ylim=y_lim,xlim=c(0,1))
points(Q_GF17 ~ Zo,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF17"])
points(Q_GF18 ~ Zo,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ Zo,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ Zo,Result[Campaign == "IC2" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(0.05,0.07,0.1),col="blue")
}
# only when MFC has a flow
for(i in names(SonicCols)){
plot(Q_GF17 ~ Zo,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("Zo",i),ylim=y_lim,xlim=c(0,0.2))
points(Q_GF18 ~ Zo,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ Zo,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ Zo,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(0.05,0.07,0.1),col="blue")
}

# ---> 0.07 or 0.1 resp. not necessary to filter
---> do not filter


##### For plotting, copy object so I can plot unfiltered data:
Res_unfil <- copy(Result)

# first some rough filter and only some parameters. This was kind of a recursive process
# Result[sUu >= 4 | sVu >= 4 | sWu >= 2 | abs(1/L) >= 0.5 | Zo >= 0.07 | C0 >= 10, grep("Q_GF",names(Result),value = TRUE) := NA_real_]
# Result[ abs(1/L) >= 0.5 |  C0 >= 10, grep("Q_GF",names(Result),value = TRUE) := NA_real_] # all the other parameters are actually not necessary and a bit arbitrary.
# above was some old filter. I rather exclude C0 and take the others.
Result[Ustar < 0.15, grep("Q_GF",names(Result),value = TRUE) := NA_real_]


## all Sonics and GFs
par(mfrow=c(2,2))
plot(Q_GF16 ~ Ustar,Res_unfil[],pch=19,col="black",ylab="Q [kg/h]",type="n",ylim=c(-5,10))
abline(h=0,v=0.15)
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Res_unfil) ~ Ustar,Res_unfil[],pch=19,col="black")}
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Result) ~ Ustar,Result[],pch=19,col="#F56868")}

plot(Q_GF16 ~ sUu,Res_unfil[],pch=19,col="black",ylab="Q [kg/h]",type="n",ylim=c(-5,10))
abline(h=0)
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Res_unfil) ~ sUu,Res_unfil[],pch=19,col="black")}
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Result) ~ sUu,Result[],pch=19,col="#F56868")}

plot(Q_GF16 ~ C0,Res_unfil[],pch=19,col="black",ylab="Q [kg/h]",type="n",ylim=c(-5,10),xlim=c(0,100))
abline(h=0)
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Res_unfil) ~ C0,Res_unfil[],pch=19,col="black")}
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Result) ~ C0,Result[],pch=19,col="#F56868")}

plot(Q_GF16 ~ Zo,Res_unfil[],pch=19,col="black",ylab="Q [kg/h]",type="n",ylim=c(-5,10),xlim=c(0,0.6))
abline(h=0)
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Res_unfil) ~ Zo,Res_unfil[],pch=19,col="black")}
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Result) ~ Zo,Result[],pch=19,col="#F56868")}


# png(file.path(PfadFigures,"Sonic_filtering.png"),width=24,height=12,unit="in",res=300)
par(mfrow=c(2,3), mar=c(4,5,2,1))
## SonicA
plot(Q_GF16 ~ Ustar,Res_unfil[Sonic == "SonicA"],pch=19,col="black",ylab="Q [kg/h]",type="n",ylim=c(-5,10),xlim=c(0,0.5))
abline(h=0,v=0.15)
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Res_unfil[Sonic == "SonicA"]) ~ Ustar,Res_unfil[Sonic == "SonicA"],pch=19,col="black")}
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Result) ~ Ustar,Result[],pch=19,col=SonicCols["SonicA"])}

plot(Q_GF16 ~ sUu,Res_unfil[Sonic == "SonicA"],pch=19,col="black",ylab="Q [kg/h]",type="n",ylim=c(-5,10),xlim=c(0,13))
abline(h=0)
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Res_unfil[Sonic == "SonicA"]) ~ sUu,Res_unfil[Sonic == "SonicA"],pch=19,col="black")}
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Result) ~ sUu,Result[],pch=19,col=SonicCols["SonicA"])}

plot(Q_GF16 ~ C0,Res_unfil[Sonic == "SonicA"],pch=19,col="black",ylab="Q [kg/h]",type="n",ylim=c(-5,10),xlim=c(0,100))
abline(h=0)
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Res_unfil[Sonic == "SonicA"]) ~ C0,Res_unfil[Sonic == "SonicA"],pch=19,col="black")}
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Result) ~ C0,Result[],pch=19,col=SonicCols["SonicA"])}

## SonicB
plot(Q_GF16 ~ Ustar,Res_unfil[Sonic == "SonicB"],pch=19,col="black",ylab="Q [kg/h]",type="n",ylim=c(-5,10),xlim=c(0,0.5))
abline(h=0,v=0.15)
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Res_unfil[Sonic == "SonicB"]) ~ Ustar,Res_unfil[Sonic == "SonicB"],pch=19,col="black")}
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Result) ~ Ustar,Result[],pch=19,col=SonicCols["SonicB"])}

plot(Q_GF16 ~ sUu,Res_unfil[Sonic == "SonicB"],pch=19,col="black",ylab="Q [kg/h]",type="n",ylim=c(-5,10),xlim=c(0,13))
abline(h=0)
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Res_unfil[Sonic == "SonicB"]) ~ sUu,Res_unfil[Sonic == "SonicB"],pch=19,col="black")}
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Result) ~ sUu,Result[],pch=19,col=SonicCols["SonicB"])}

plot(Q_GF16 ~ C0,Res_unfil[Sonic == "SonicB"],pch=19,col="black",ylab="Q [kg/h]",type="n",ylim=c(-5,10),xlim=c(0,100))
abline(h=0)
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Res_unfil[Sonic == "SonicB"]) ~ C0,Res_unfil[Sonic == "SonicB"],pch=19,col="black")}
for(i in c("Q_GF16","Q_GF17","Q_GF18","Q_GF25")){
points(get(i,Result) ~ C0,Result[],pch=19,col=SonicCols["SonicB"])}

# dev.off()


#########################
### look at emissions ###
#########################

### all emissions with MFC. MC and IC2
par(mfrow=c(3,1))
# Emiss
plot(Q_GF16 ~ st, Result[Campaign %in% c("MC","IC2")],type="l",col=CH4Cols["GF16"],ylim=c(-10,10),ylab="Q [kg/h]")
lines(Q_GF17 ~ st, Result[Campaign %in% c("MC","IC2")],col=CH4Cols["GF17"])
lines(Q_GF18 ~ st, Result[Campaign %in% c("MC","IC2")],col=CH4Cols["GF18"])
lines(Q_GF25 ~ st, Result[Campaign %in% c("MC","IC2")],col=CH4Cols["GF25"])
lines(Q_GF26 ~ st, Result[Campaign %in% c("MC","IC2")],col=CH4Cols["GF26"])
lines(Q_MFC ~ st, Result[Campaign %in% c("MC","IC2")],col="black",lwd=2)
# wind speed
plot(U_sonic ~ st, Result[Campaign %in% c("MC","IC2") & Sonic == "SonicB"],type="l",col="blue",ylab="U_sonic [m/s]")
# wind direction
plot(WD_WS2 ~ st, Result[Campaign %in% c("MC","IC2") & Sonic == "SonicB"],type="l",col="red",ylab="wind direction")

### all emissions with MFC. only MC
par(mfrow=c(3,1))
# Emiss
plot(Q_GF16 ~ st, Result[Campaign %in% c("MC")],type="l",col=CH4Cols["GF16"],ylim=c(-10,10),ylab="Q [kg/h]")
lines(Q_GF17 ~ st, Result[Campaign %in% c("MC")],col=CH4Cols["GF17"])
lines(Q_GF18 ~ st, Result[Campaign %in% c("MC")],col=CH4Cols["GF18"])
lines(Q_GF25 ~ st, Result[Campaign %in% c("MC")],col=CH4Cols["GF25"])
lines(Q_GF26 ~ st, Result[Campaign %in% c("MC")],col=CH4Cols["GF26"])
lines(Q_MFC ~ st, Result[Campaign %in% c("MC")],col="black",lwd=2)
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])

# wind speed
plot(U_sonic ~ st, Result[Campaign %in% c("MC") & Sonic == "SonicB"],type="l",col="blue",ylab="U_sonic [m/s]")
# wind direction
plot(WD_WS2 ~ st, Result[Campaign %in% c("MC") & Sonic == "SonicB"],type="l",col="red",ylab="wind direction")

graphics.off()


####################################################
### Filter again the WD as it is not good enough ###
####################################################

par(mfrow=c(2,5),mar=c(4,4,0,0))
plot(Q_GF17 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MC"],xlab="WD all Sonics", col= CH4Cols["GF17"], ylim = c(-20,10),xlim=c(-10,110))
abline(h=0, v=(iWD17 + 180) %% 360 - 180)
abline(v=c(21),lty=3,col="red")
for(i in names(SonicCols)){
plot(Q_GF17 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MC" & Sonic == i],xlab=paste0("WD ", i), col= CH4Cols["GF17"], ylim = c(-20,10),xlim=c(-10,110))
abline(h=0, v=(iWD17 + 180) %% 360 - 180)
abline(v=c(21),lty=3,col="red")
}
plot(Q_GF18 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MC"],xlab="WD all Sonics", col= CH4Cols["GF18"], ylim = c(-20,10),xlim=c(-10,110))
abline(h=0, v=(iWD18 + 180) %% 360 - 180)
abline(v=c(21,71),lty=3,col="red")
for(i in names(SonicCols)){
plot(Q_GF18 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MC" & Sonic == i],xlab=paste0("WD ", i), col= CH4Cols["GF18"], ylim = c(-20,10),xlim=c(-10,110))
abline(h=0, v=(iWD18 + 180) %% 360 - 180)
abline(v=c(21,71),lty=3,col="red")
}

par(mfrow=c(2,5),mar=c(4,4,0,0))
plot(Q_GF16 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MC"],xlab="WD all Sonics", col= CH4Cols["GF16"], ylim = c(-20,10),xlim=c(-10,110))
abline(h=0, v=(iWD16 + 180) %% 360 - 180)
abline(v=c(30,64),lty=3,col="red")
for(i in names(SonicCols)){
plot(Q_GF16 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MC" & Sonic == i],xlab=paste0("WD ", i), col= CH4Cols["GF16"], ylim = c(-20,10),xlim=c(-10,110))
abline(h=0, v=(iWD16 + 180) %% 360 - 180)
abline(v=c(30,64),lty=3,col="red")
}
plot(Q_GF25 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MC"],xlab="WD all Sonics", col= CH4Cols["GF25"], ylim = c(-20,10),xlim=c(-10,110))
abline(h=0, v=(iWD25 + 180) %% 360 - 180)
abline(v=c(30,61),lty=3,col="red")
for(i in names(SonicCols)){
plot(Q_GF25 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MC" & Sonic == i],xlab=paste0("WD ", i), col= CH4Cols["GF25"], ylim = c(-20,10),xlim=c(-10,110))
abline(h=0, v=(iWD25 + 180) %% 360 - 180)
abline(v=c(30,61),lty=3,col="red")
}

## according to map
iWD_refined <- list(
	iWD_GF16 = c(30,64)
	,iWD_GF17 = c(21,iWD17[2])
	,iWD_GF18 = c(21,71)
	,iWD_GF25 = c(30,61)
)


### Filter according to map and corresponding sonics as a reference for the wind direction

# MC
for(x in c("GF16","GF17","GF18","GF25")){
	Result[Campaign == "MC" & (WD < iWD_refined[[paste0("iWD_",x)]][1] | WD > iWD_refined[[paste0("iWD_",x)]][2]),
		grep(paste0("Q_",x),names(Result),value=TRUE)	 := NA_real_ ]
}

# # IC2 # this is not done and probably not necessary
# for(x in c("GF16","GF17","GF18","GF25","GF26")){
# 	Result[Campaign == "IC2" & (WD < iWDIC2[1] - 5 & WD > iWDIC2[2] + 5),
# 		grep(paste0("Q_",x),names(Result),value=TRUE)	 := NA_real_ ]
# }

---> maybe filter also CQ?


#########################
### Look at CQ values ###
#########################

CQ_thresh <- 8E-5

par(mfrow=c(2,5),mar=c(4,4,0,0))
plot(Q_GF17 ~ CQ_GF17,Result[Campaign == "MC"],xlab=paste0("CQ_GF17 ",i), col= CH4Cols["GF17"], ylim = c(-20,10),xlim=c(0,0.0015))
abline(h=0, v=CQ_thresh)
for(i in names(SonicCols)){
plot(Q_GF17 ~ CQ_GF17,Result[Campaign == "MC" & Sonic == i],xlab=paste0("CQ_GF17 ", i), col= CH4Cols["GF17"], ylim = c(-20,10),xlim=c(0,0.0015))
abline(h=0, v=CQ_thresh)
}
plot(Q_GF18 ~ CQ_GF18,Result[Campaign == "MC"],xlab=paste0("CQ_GF18 ",i), col= CH4Cols["GF18"], ylim = c(-20,10),xlim=c(0,0.0015))
abline(h=0, v=CQ_thresh)
for(i in names(SonicCols)){
plot(Q_GF18 ~ CQ_GF18,Result[Campaign == "MC" & Sonic == i],xlab=paste0("CQ_GF18 ", i), col= CH4Cols["GF18"], ylim = c(-20,10),xlim=c(0,0.0015))
abline(h=0, v=CQ_thresh)
}

par(mfrow=c(2,5),mar=c(4,4,0,0))
plot(Q_GF16 ~ CQ_GF16,Result[Campaign == "MC"],xlab=paste0("CQ_GF16 ",i), col= CH4Cols["GF16"], ylim = c(-20,10),xlim=c(0,0.0015))
abline(h=0, v=CQ_thresh)
for(i in names(SonicCols)){
plot(Q_GF16 ~ CQ_GF16,Result[Campaign == "MC" & Sonic == i],xlab=paste0("CQ_GF16 ", i), col= CH4Cols["GF16"], ylim = c(-20,10),xlim=c(0,0.0015))
abline(h=0, v=CQ_thresh)
}
plot(Q_GF25 ~ CQ_GF25,Result[Campaign == "MC"],xlab=paste0("CQ_GF25 ",i), col= CH4Cols["GF25"], ylim = c(-20,10),xlim=c(0,0.0015))
abline(h=0, v=CQ_thresh)
for(i in names(SonicCols)){
plot(Q_GF25 ~ CQ_GF25,Result[Campaign == "MC" & Sonic == i],xlab=paste0("CQ_GF25 ", i), col= CH4Cols["GF25"], ylim = c(-20,10),xlim=c(0,0.0015))
abline(h=0, v=CQ_thresh)
}

---> CQ is not reasonable as for each GF and Sonic another threshold needs to be used, which is too arbitary.

graphics.off()


#################
### save data ###
#################

saveRDS(Result,file=file.path(PathRSaves,"Results_orig_P23.rds"))



#################################################################
#################################################################
#####                                                       #####
#####    Apply filter to wind direction variation object    #####
#####                                                       #####
#################################################################
#################################################################

## this is only needed to make the plot in the initial submission and has otherwise no furhter use.
## It can also be used to filter the different concentration offset corrections.
## Belwo some examples are given. Just adapt the code to whatever you calculated in the script 02_Calculation_03_Emissions.r


########################################
### bind the different runs together ###
########################################

## for Concentration offset correction additional column is needed
Emiss_Result_orig_wo[,Conc_corr := "wo"]
Emiss_Result_orig_P23[,Conc_corr := "P23"]
Emiss_Result_orig_P6[,Conc_corr := "P6"]
Emiss_Result_orig_P4[,Conc_corr := "P4"]
Emiss_Result_WDvar_P23[,Conc_corr := "P23"]
Emiss_Result_WDvar_P6[,Conc_corr := "P6"]

Emiss_Result_var <- rbind(Emiss_Result_orig_wo,Emiss_Result_orig_P23,Emiss_Result_orig_P6,Emiss_Result_orig_P4,Emiss_Result_WDvar_P23,Emiss_Result_WDvar_P6,fill=TRUE)
setkey(Emiss_Result_var,st)

## make a column to indicate the start and end of the CH4 release
Emiss_Result_var[,REL := FALSE]
Emiss_Result_var[st > parse_date_time3("19.03.2021 09:00",tz="Etc/GMT-1") & st < parse_date_time3("20.03.2021 09:00",tz="Etc/GMT-1"),REL := TRUE]

## only use data with NE wind
Result_var <- Emiss_Result_var[Wind_dir == "NE"]

#########################
### Plot and overview ###
#########################

Sonic_conc <- c("SonicA","SonicB","SonicC","Sonic2")
Sonic_WDvar <- paste0(Sonic_conc,rep(c(paste0("_m",1:10),"",paste0("_p",1:10)),4))

#### Concentration offset
plot(Q_GF16 ~ st, Result_var[REL == TRUE & Sonic %in% Sonic_conc],type="l",col=CH4Cols["GF16"],ylim=c(-10,10))
lines(Q_GF17 ~ st, Result_var[REL == TRUE & Sonic %in% Sonic_conc],col=CH4Cols["GF17"])
lines(Q_GF18 ~ st, Result_var[REL == TRUE & Sonic %in% Sonic_conc],col=CH4Cols["GF18"])
lines(Q_GF25 ~ st, Result_var[REL == TRUE & Sonic %in% Sonic_conc],col=CH4Cols["GF25"])
lines(Q_GF26 ~ st, Result_var[REL == TRUE & Sonic %in% Sonic_conc],col=CH4Cols["GF26"])
lines(Q_MFC ~ st, Result_var[REL == TRUE & Sonic %in% Sonic_conc],col="black",lwd=2)

#### Wind direction offset
plot(Q_GF16 ~ st, Result_var[REL == TRUE & Sonic %in% Sonic_WDvar & Conc_corr == "P23"],type="l",col=CH4Cols["GF16"],ylim=c(-10,10))
lines(Q_GF17 ~ st, Result_var[REL == TRUE & Sonic %in% Sonic_WDvar & Conc_corr == "P23"],col=CH4Cols["GF17"])
lines(Q_GF18 ~ st, Result_var[REL == TRUE & Sonic %in% Sonic_WDvar & Conc_corr == "P23"],col=CH4Cols["GF18"])
lines(Q_GF25 ~ st, Result_var[REL == TRUE & Sonic %in% Sonic_WDvar & Conc_corr == "P23"],col=CH4Cols["GF25"])
lines(Q_GF26 ~ st, Result_var[REL == TRUE & Sonic %in% Sonic_WDvar & Conc_corr == "P23"],col=CH4Cols["GF26"])
lines(Q_MFC ~ st, Result_var[REL == TRUE & Sonic %in% Sonic_WDvar & Conc_corr == "P23"],col="black",lwd=2)


#####################
### Apply filters ###
#####################

------> use the filtering from the original version. There should not be any changes necessary.

iWD16 <- c(30,64)
iWD17 <- c(21,94)
iWD18 <- c(21,71)
iWD25 <- c(30,61)

iWD <- list(
	iWD_GF16 = c(iWD16[1], iWD16[2])
	,iWD_GF17 = c(iWD17[1], iWD17[2])
	,iWD_GF18 = c(iWD18[1], iWD18[2])
	,iWD_GF25 = c(iWD25[1], iWD25[2])
)

iWDIC2 <- c(357,93)

### Wind direction
# MC
for(x in c("GF16","GF17","GF18","GF25")){
	Result_var[Campaign == "MC" & (WD < iWD[[paste0("iWD_",x)]][1] | WD > iWD[[paste0("iWD_",x)]][2]),
		grep(paste0("Q_",x),names(Result_var),value=TRUE)	 := NA_real_ ]
}

# IC2
for(x in c("GF16","GF17","GF18","GF25","GF26")){
	Result_var[Campaign == "IC2" & (WD_WS2 < iWDIC2[1] - 5 & WD_WS2 > iWDIC2[2] + 5),
		grep(paste0("Q_",x),names(Result_var),value=TRUE)	 := NA_real_ ]
}


## bLS Filtern
Result_var[Ustar < 0.15, grep("Q_GF",names(Result_var),value = TRUE) := NA_real_]


#########################
### Plotting the data ###
#########################

#########################
### Concentration offset:

par(mfrow=c(4,1),mar=c(3,4,2,0.2))
# GF16
plot(Q_GF16 ~ st, Result_var[REL == TRUE & Sonic %in% "SonicC" & Conc_corr == "wo"],type="l",col="red",ylim=c(-2,8),main="GF16",ylab="Q [kg/h]")
lines(Q_GF16 ~ st, Result_var[REL == TRUE & Sonic %in% "SonicC" & Conc_corr == "P23"],col="green")
lines(Q_GF16 ~ st, Result_var[REL == TRUE & Sonic %in% "SonicC" & Conc_corr == "P6"],col="blue")
lines(Q_GF16 ~ st, Result_var[REL == TRUE & Sonic %in% "SonicC" & Conc_corr == "P4"],col="pink")
lines(Q_MFC ~ st, Result_var[],col="black",lwd=2)
abline(h=0)
legend("topright",legend=c("wo","P23","P6",'P4'),fill=c("red","green","blue",'pink'))
# GF17
plot(Q_GF17 ~ st, Result_var[REL == TRUE & Sonic %in% "SonicC" & Conc_corr == "wo"],type="l",col="red",ylim=c(-2,8),main="GF17",ylab="Q [kg/h]")
lines(Q_GF17 ~ st, Result_var[REL == TRUE & Sonic %in% "SonicC" & Conc_corr == "P23"],col="green")
lines(Q_GF17 ~ st, Result_var[REL == TRUE & Sonic %in% "SonicC" & Conc_corr == "P6"],col="blue")
lines(Q_GF17 ~ st, Result_var[REL == TRUE & Sonic %in% "SonicC" & Conc_corr == "P4"],col="blue")
lines(Q_MFC ~ st, Result_var[],col="black",lwd=2)
abline(h=0)
legend("topright",legend=c("wo","P23","P6",'P4'),fill=c("red","green","blue",'pink'))
# GF18
plot(Q_GF18 ~ st, Result_var[REL == TRUE & Sonic %in% "SonicC" & Conc_corr == "wo"],type="l",col="red",ylim=c(-2,8),main="GF18",ylab="Q [kg/h]")
lines(Q_GF18 ~ st, Result_var[REL == TRUE & Sonic %in% "SonicC" & Conc_corr == "P23"],col="green")
lines(Q_GF18 ~ st, Result_var[REL == TRUE & Sonic %in% "SonicC" & Conc_corr == "P6"],col="blue")
lines(Q_GF18 ~ st, Result_var[REL == TRUE & Sonic %in% "SonicC" & Conc_corr == "P4"],col="blue")
lines(Q_MFC ~ st, Result_var[],col="black",lwd=2)
abline(h=0)
legend("topright",legend=c("wo","P23","P6",'P4'),fill=c("red","green","blue",'pink'))
# GF25
plot(Q_GF25 ~ st, Result_var[REL == TRUE & Sonic %in% "SonicC" & Conc_corr == "wo"],type="l",col="red",ylim=c(-2,8),main="GF25",ylab="Q [kg/h]")
lines(Q_GF25 ~ st, Result_var[REL == TRUE & Sonic %in% "SonicC" & Conc_corr == "P23"],col="green")
lines(Q_GF25 ~ st, Result_var[REL == TRUE & Sonic %in% "SonicC" & Conc_corr == "P6"],col="blue")
lines(Q_GF25 ~ st, Result_var[REL == TRUE & Sonic %in% "SonicC" & Conc_corr == "P4"],col="blue")
lines(Q_MFC ~ st, Result_var[],col="black",lwd=2)
abline(h=0)
legend("topright",legend=c("wo","P23","P6",'P4'),fill=c("red","green","blue",'pink'))

# ----> does not matter if P23 or P6 is used. Therefore, us P23. Above only the Emission with bLS SonicC (UA-UW)are given
# --> looking at IC2 is not necessary. There were no changes.


#############################
### Wind direction variation:

WDmCols <- colorRampPalette(c("#00ff00", "#003300"))(10)
WDpCols <- colorRampPalette(c("#00ffff", "#0033ff"))(10)
WDvarCols <- c(WDmCols,"red",WDpCols)
names(WDvarCols) <- paste0("SonicB",c(paste0("_m",1:10),"_0",paste0("_p",1:10)))

## GF17
par(mfrow=c(4,1),mar=c(3,4,2,0.2))
for(j in c("SonicA","Sonic2","SonicB","SonicC")){
plot(Q_GF17 ~ st, Result_var[REL == TRUE & Sonic %in% j & Conc_corr == "P23"],type="l",col="red",ylim=c(-2,10),main=j
	,ylab="Q GF17 [kg/h]")
for(i in 1:10){
	lines(Q_GF17 ~ st, Result_var[REL == TRUE & Sonic %in% paste0(j,"_m",i) & Conc_corr == "P23"],type="o",col=WDmCols[i])
	lines(Q_GF17 ~ st, Result_var[REL == TRUE & Sonic %in% paste0(j,"_p",i) & Conc_corr == "P23"],type="o",col=WDpCols[i])
}
lines(Q_MFC ~ st, Result_var,col="black",lwd=2)
abline(h=0)
legend("bottomright",legend=c("negative","neutral","positive"),fill=c(WDmCols[8],"red",WDpCols[8]))
}

## GF18
par(mfrow=c(4,1),mar=c(3,4,2,0.2))
for(j in c("SonicA","Sonic2","SonicB","SonicC")){
plot(Q_GF18 ~ st, Result_var[REL == TRUE & Sonic %in% j & Conc_corr == "P23"],type="l",col="red",ylim=c(-2,10),main=j
	,ylab="Q GF18 [kg/h]")
for(i in 1:10){
	lines(Q_GF18 ~ st, Result_var[REL == TRUE & Sonic %in% paste0(j,"_m",i) & Conc_corr == "P23"],col=WDmCols[i])
	lines(Q_GF18 ~ st, Result_var[REL == TRUE & Sonic %in% paste0(j,"_p",i) & Conc_corr == "P23"],col=WDpCols[i])
}
lines(Q_MFC ~ st, Result_var,col="black",lwd=2)
abline(h=0)
legend("bottomright",legend=c("negative","neutral","positive"),fill=c(WDmCols[8],"red",WDpCols[8]))
}

## GF16
par(mfrow=c(4,1),mar=c(3,4,2,0.2))
for(j in c("SonicA","Sonic2","SonicB","SonicC")){
plot(Q_GF16 ~ st, Result_var[REL == TRUE & Sonic %in% j & Conc_corr == "P23"],type="l",col="red",ylim=c(-2,10),main=j
	,ylab="Q GF16 [kg/h]")
for(i in 1:10){
	lines(Q_GF16 ~ st, Result_var[REL == TRUE & Sonic %in% paste0(j,"_m",i) & Conc_corr == "P23"],col=WDmCols[i])
	lines(Q_GF16 ~ st, Result_var[REL == TRUE & Sonic %in% paste0(j,"_p",i) & Conc_corr == "P23"],col=WDpCols[i])
}
lines(Q_MFC ~ st, Result_var,col="black",lwd=2)
abline(h=0)
legend("bottomright",legend=c("negative","neutral","positive"),fill=c(WDmCols[8],"red",WDpCols[8]))
}

## GF25
par(mfrow=c(4,1),mar=c(3,4,2,0.2))
for(j in c("SonicA","Sonic2","SonicB","SonicC")){
plot(Q_GF25 ~ st, Result_var[REL == TRUE & Sonic %in% j & Conc_corr == "P23"],type="l",col="red",ylim=c(-2,10),main=j
	,ylab="Q GF25 [kg/h]")
for(i in 1:10){
	lines(Q_GF25 ~ st, Result_var[REL == TRUE & Sonic %in% paste0(j,"_m",i) & Conc_corr == "P23"],col=WDmCols[i])
	lines(Q_GF25 ~ st, Result_var[REL == TRUE & Sonic %in% paste0(j,"_p",i) & Conc_corr == "P23"],col=WDpCols[i])
}
lines(Q_MFC ~ st, Result_var,col="black",lwd=2)
abline(h=0)
legend("bottomright",legend=c("negative","neutral","positive"),fill=c(WDmCols[8],"red",WDpCols[8]))
}


## all GasFinders and SonicB
par(mfrow=c(4,1),mar=c(3,4,2,0.2))
for(j in c("GF17","GF18","GF16","GF25")){
plot(get(paste0("Q_",j),Result_var[REL == TRUE & Sonic=="SonicB" & Conc_corr == "P23"]) ~ st
	, Result_var[REL == TRUE & Sonic == "SonicB" & Conc_corr == "P23"],type="o",col="red",ylim=c(-2,10), ylab=paste0("Q ",j," [kg/h]"))
for(i in 1:10){
	lines(get(paste0("Q_",j),Result_var[REL == TRUE & Sonic == paste0("SonicB_m",i) & Conc_corr == "P23"]) ~ st
		, Result_var[REL == TRUE & Sonic == paste0("SonicB_m",i) & Conc_corr == "P23"],type="o",col=WDmCols[i])
	lines(get(paste0("Q_",j),Result_var[REL == TRUE & Sonic == paste0("SonicB_p",i) & Conc_corr == "P23"]) ~ st
		, Result_var[REL == TRUE & Sonic == paste0("SonicB_p",i) & Conc_corr == "P23"],type="o",col=WDpCols[i])
}
lines(Q_MFC ~ st, Result_var,col="black",lwd=2)
abline(h=0)
# legend("bottomright",legend=c("negative","neutral","positive"),fill=c(WDmCols[8],"red",WDpCols[8]))
}


## all GasFinders and SonicC
par(mfrow=c(4,1),mar=c(3,4,2,0.2))
for(j in c("GF17","GF18","GF16","GF25")){
plot(get(paste0("Q_",j),Result_var[REL == TRUE & Sonic=="SonicC" & Conc_corr == "P23"]) ~ st
	, Result_var[REL == TRUE & Sonic == "SonicC" & Conc_corr == "P23"],type="o",col="red",ylim=c(-2,10), ylab=paste0("Q ",j," [kg/h]"))
for(i in 1:10){
	lines(get(paste0("Q_",j),Result_var[REL == TRUE & Sonic == paste0("SonicC_m",i) & Conc_corr == "P23"]) ~ st
		, Result_var[REL == TRUE & Sonic == paste0("SonicC_m",i) & Conc_corr == "P23"],type="o",col=WDmCols[i])
	lines(get(paste0("Q_",j),Result_var[REL == TRUE & Sonic == paste0("SonicC_p",i) & Conc_corr == "P23"]) ~ st
		, Result_var[REL == TRUE & Sonic == paste0("SonicC_p",i) & Conc_corr == "P23"],type="o",col=WDpCols[i])
}
lines(Q_MFC ~ st, Result_var,col="black",lwd=2)
abline(h=0)
# legend("bottomright",legend=c("negative","neutral","positive"),fill=c(WDmCols[8],"red",WDpCols[8]))
}


#################
### save data ###
#################

saveRDS(Result_var,file.path(PfadRSaves,"Result_variation.rds"))

