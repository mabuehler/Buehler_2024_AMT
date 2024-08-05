
##############################################################
##############################################################
#####                                                    #####
#####    Further data treatment of CH4 concentrations    #####
#####                                                    #####
##############################################################
##############################################################

# Author: Marcel Bühler
# Date: August 5, 2024
# Contact: mb@bce.au.dk or Christoph Häni christoph.haeni@bfh.ch
# Description: This script treats the unprocessed concentration data. It removes false concentrations, applies an intercalibration, and makes the data ready for further use.
#
# Note: This code was written by Marcel Bühler and is intended to follow the publication 'Applicability of the inverse dispersion method to measure emissions from animal housings' in AMT. 
# Please feel free to use and modify it, but attribution is appreciated.


#################
### Libraries ###
#################

library(ibts)
library(RgoogleMaps)
library(bLSmodelR)
library(ggplot2)
library(ggthemes)
library(ggpointdensity)
library(ggpubr)
library(RColorBrewer)


#############
### Paths ###
#############

PathData <- "Path to /data"		
PathRSaves <- "Path to /RSaves"
PathFigures <- "~/Path to /Figures"


#################
### Functions ###
#################

source("https://raw.githubusercontent.com/hafl-gel/gel-scripts/main/wgs84-ch1903.r")

lines_sec2xy <- function(sensor,wd,col="lightblue",lwd=2,...){
	sens <- as.numeric(Sensors_xy[Sensors_xy[,1] %in% sensor,2:3])
	b <- tan((90 - wd)/180*pi)
	x <- if(wd <= 180) 600 else -600
	y <- sens[2] - (sens[1] - x)*b
	lines(c(sens[1],x),c(sens[2],y),col=col,lwd=lwd,...)
}


#################
### Campaigns ###
#################

StartMeas <- "05.03.2021"
StopMeas <- "26.03.2021 12:00"
IC1 <- "05.03.2021 to 10.03.2021 18:00"
MC <- "18.03.2021 11:00 - 21.03.2021 14:00"
IC2	<- "21.03.2021 14:00 to "
indMCs <- deparse_timerange(StartMeas, StopMeas)

#################
### load data ###
#################

	# Geometry
	load(file.path(PathRSaves,"Geometry.RData"))
	# Google map
	center <- c(47.043053, 7.226895)
	STO_Map <- ReadMapTile(file.path(PathFigures,"STO_GoogleMaps.png"))
	# GF data
	load(file.path(PathRSaves,"GFraw.RData"))
	# Sonic data
	load(file.path(PathRSaves,"Sonics_10min.RData"))
	# WS700
	WS700 <- readRDS(file=file.path(PathRSaves,"WS700.rds"))


############################
### Colours for plotting ###
############################

# display.brewer.all()
CH4Cols <- brewer.pal(5,"Dark2")
names(CH4Cols) <- paste0("GF",c(16:18,25,26))

VerCols <- brewer.pal(6,"Set1")
names(VerCols) <- c("IC1","IC2","IC12","IC12E","IC2E","MC")


#########################
### Quality filtering ###
#########################

	rp_gf16 <- 70
	rp_gf17_IC1 <- 600
	rp_gf17_MC <- 400
	rp_gf17_IC2 <- 100
	rp_gf18 <- 120
	rp_gf25 <- 50
	rp_gf26 <- 200 # 50 would also be possible
	r2_thresh <- 98


#############
### GF-16 ###
#############

	GF16_revised <- GF_16[which(
		GF_16$received_power > rp_gf16 & 
		GF_16$R2 > r2_thresh
		),][indMCs]

	## IC1
	# x11(width=12)
	# par(oma=c(0,0,0,2))
	# plot(GF_16[IC1,c("CH4_mgm3","received_power")],col2="lightgrey",col=NA,type2="l",ylim=c(1,2))
	# lines(GF_16[IC1,"CH4_mgm3"], col = "indianred")
	# lines(GF16_revised[IC1,"CH4_mgm3"],col="blue")
	# abline(v=parse_date_time3("05.03.2021 15:05",tz="Etc/GMT-1"))
	# abline(v=parse_date_time3(c("07.03.2021 12:00","07.03.2021 18:00")),col="red") # in comparison with the others there are rather low conc. but nothing suspiceous here
	# abline(v=parse_date_time3(c("08.03.2021 12:00","08.03.2021 18:00")),col="red") # in comparison with the others there are rather low conc. but nothing suspiceous here
	# abline(v=parse_date_time3(c("09.03.2021 12:00","09.03.2021 18:00")),col="red") # in comparison with the others there are rather low conc. but nothing suspiceous here
	# --> remove start
	GF16_revised["05.03.2021 to 05.03.2021 15:05",c("CH4_mgm3","CH4_ppm")] <- NA_real_
	## MC
	# x11(width=12)
	# par(oma=c(0,0,0,2))
	# plot(GF_16[MC,c("CH4_mgm3","received_power")],col2="lightgrey",col=NA,type2="l",ylim=c(1,2))
	# lines(GF_16[MC,"CH4_mgm3"], col = "indianred")
	# lines(GF16_revised[MC,"CH4_mgm3"],col="blue")
	# --> remove start ???  ---> wahrscheinlich nicht. Alle anderen zeigen das gleiche Verhalten. Haben wir dort Gas freigesetzt??
	# x11(width=12)
	# par(oma=c(0,0,0,2))
	# plot(GF_16[IC1,c("CH4_mgm3","received_power")],col2="lightgrey",col=NA,type2="l",ylim=c(1,2))
	# lines(GF_16[IC1,"CH4_mgm3"], col = "indianred")
	# lines(GF16_revised[IC1,"CH4_mgm3"],col="blue")
	# # --> looks good


#############
### GF-17 ###
#############

	## IC1
	GF17_revised_IC1 <- GF_17[which(
		GF_17$received_power > rp_gf17_IC1 & 
		GF_17$R2 > r2_thresh
		),][IC1]
	# x11(width=12)
	# par(oma=c(0,0,0,2))
	# plot(GF_17[IC1,c("CH4_mgm3","received_power")],col2="lightgrey",col=NA,type2="l")
	# lines(GF_17[IC1,"CH4_mgm3"], col = "indianred")
	# lines(GF17_revised_IC1[,"CH4_mgm3"], col="blue")
	## noch mehr ausschliessen. Siehe GF16
	# abline(v=parse_date_time3(c("05.03.2021","05.03.2021 15:05",tz="Etc/GMT-1")))
	# GF17_revised_QV1["05.03.2021 to 05.03.2021 15:05",c("CH4_mgm3","CH4_ppm")] <- NA_real_
	## MC
	GF17_revised_MC <- GF_17[which(
		GF_17$received_power > rp_gf17_MC & 
		GF_17$R2 > r2_thresh
		),][MC]
	#
	# x11(width=12)
	# par(oma=c(0,0,0,2))
	# plot(GF_17[MC,c("CH4_mgm3","received_power")],col2="lightgrey",col=NA,type2="l",ylim=c(1,2))
	# lines(GF_17[MC,"CH4_mgm3"], col = "indianred")
	# lines(GF17_revised[MC,"CH4_mgm3"],col="blue")
	# --> remove start
	# abline(v=parse_date_time3("18.03.2021 13:50",tz="Etc/GMT-1"))
	GF17_revised_MC["17.03.2021 to 18.03.2021 13:50",c("CH4_mgm3","CH4_ppm")] <- NA_real_
	## IC2
	GF17_revised_IC2 <- GF_17[which(
		GF_17$received_power > rp_gf17_IC2 & 
		GF_17$R2 > r2_thresh
		),][IC2]
	# #
	# x11(width=12)
	# par(oma=c(0,0,0,2))
	# plot(GF_17[IC2,c("CH4_mgm3","received_power")],col2="lightgrey",col=NA,type2="l",ylim=c(1,2))
	# lines(GF_17[IC2,"CH4_mgm3"], col = "indianred")
	# lines(GF17_revised_IC2[,"CH4_mgm3"],col="blue")
	# # einiges ausschliessen. höherer Threshold ist aber nicht sinnvoll, da sonst zuviel wegfällt.
	# abline(v=parse_date_time3("21.03.2021 15:00",tz="Etc/GMT-1"))
	# abline(v=parse_date_time3(c("24.03.2021 07:00","24.03.2021 09:00"),tz="Etc/GMT-1"))
	# abline(v=parse_date_time3(c("25.03.2021 08:00","25.03.2021 09:00"),tz="Etc/GMT-1"))
	# abline(v=parse_date_time3(c("26.03.2021 07:00","26.03.2021 08:15"),tz="Etc/GMT-1"))
	GF17_revised_IC2[" to 21.03.2021 15:00",c("CH4_mgm3","CH4_ppm")] <- NA_real_
	GF17_revised_IC2["24.03.2021 07:00 to 24.03.2021 09:00",c("CH4_mgm3","CH4_ppm")] <- NA_real_
	GF17_revised_IC2["25.03.2021 07:00 to 25.03.2021 09:00",c("CH4_mgm3","CH4_ppm")] <- NA_real_
	GF17_revised_IC2["26.03.2021 07:00 to 26.03.2021 08:25",c("CH4_mgm3","CH4_ppm")] <- NA_real_
	## add together
	GF17_revised <- rbind(GF17_revised_IC1,GF17_revised_MC,GF17_revised_IC2)


#############
### GF-18 ###
#############

	GF18_revised <- GF_18[which(
		GF_18$received_power > rp_gf18 & 
		GF_18$R2 > r2_thresh
		),][indMCs]
	
	## IC1
	# x11(width=12)
	# par(oma=c(0,0,0,2))
	# plot(GF_18[IC1,c("CH4_mgm3","received_power")],col2="lightgrey",col=NA,type2="l")
	# lines(GF_18[IC1,"CH4_mgm3"], col = "indianred")
	# lines(GF18_revised[IC1,"CH4_mgm3"], col="blue")
	# ---> alles ok
	## MC
	# x11(width=12)
	# par(oma=c(0,0,0,2))
	# plot(GF_18[MC,c("CH4_mgm3","received_power")],col2="lightgrey",col=NA,type2="l",ylim=c(1,2))
	# lines(GF_18[MC,"CH4_mgm3"], col = "indianred")
	# lines(GF18_revised[MC,"CH4_mgm3"],col="blue")
	# abline(v=parse_date_time3("18.03.2021 13:35",tz="Etc/GMT-1"))
	# --> erste paar Intervalle ausschliessen, der Anfang ist aber bei jedem Gerät so
	GF18_revised["18.03.2021 to 18.03.2021 13:35",c("CH4_mgm3","CH4_ppm")] <- NA_real_
	# x11(width=12)
	# par(oma=c(0,0,0,2))
	# plot(GF_18[IC2,c("CH4_mgm3","received_power")],col2="lightgrey",col=NA,type2="l",ylim=c(1,2))
	# lines(GF_18[IC2,"CH4_mgm3"], col = "indianred")
	# lines(GF18_revised[IC2,"CH4_mgm3"],col="blue")
	# # einiges ausschliessen. höherer Threshold ist aber nicht sinnvoll, da sonst zuviel wegfällt.
	# abline(v=parse_date_time3("21.03.2021 15:00",tz="Etc/GMT-1"))
	# abline(v=parse_date_time3(c("24.03.2021 07:00","24.03.2021 08:50"),tz="Etc/GMT-1"))
	# abline(v=parse_date_time3(c("25.03.2021 08:00","25.03.2021 09:00"),tz="Etc/GMT-1"))
	# abline(v=parse_date_time3(c("26.03.2021 08:55","26.03.2021 10:00"),tz="Etc/GMT-1"))
	GF18_revised["21.03.2021 to 21.03.2021 15:00",c("CH4_mgm3","CH4_ppm")] <- NA_real_
	GF18_revised["24.03.2021 07:00 to 24.03.2021 08:50",c("CH4_mgm3","CH4_ppm")] <- NA_real_
	GF18_revised["25.03.2021 08:00 to 25.03.2021 09:00",c("CH4_mgm3","CH4_ppm")] <- NA_real_
	GF18_revised["26.03.2021 08:40 to 26.03.2021 10:00",c("CH4_mgm3","CH4_ppm")] <- NA_real_


#############
### GF-25 ###
#############

	GF25_revised <- GF_25[which(
		GF_25$received_power > rp_gf25 & 
		GF_25$R2 > r2_thresh
		),][indMCs]
	
	## IC1
	# x11(width=12)
	# par(oma=c(0,0,0,2))
	# plot(GF_25[IC1,c("CH4_mgm3","received_power")],col2="lightgrey",col="indianred",type2="l")
	# lines(GF25_revised[IC1,"CH4_mgm3"], col="blue")
	# abline(v=parse_date_time3("05.03.2021 15:25",tz="Etc/GMT-1"))
	GF25_revised["05.03.2021 to 05.03.2021 15:25",c("CH4_mgm3","CH4_ppm")] <- NA_real_
	## MC
	# x11(width=12)
	# par(oma=c(0,0,0,2))
	# plot(GF_25[MC,c("CH4_mgm3","received_power")],col2="lightgrey",col=NA,type2="l",ylim=c(1,2))
	# lines(GF_25[MC,"CH4_mgm3"], col = "indianred")
	# lines(GF25_revised[MC,"CH4_mgm3"],col="blue")
	# --> sieht gut aus
	## IC2
	# x11(width=12)
	# par(oma=c(0,0,0,2))
	# plot(GF_25[IC2,c("CH4_mgm3","received_power")],col2="lightgrey",col=NA,type2="l",ylim=c(1,2))
	# lines(GF_25[IC2,"CH4_mgm3"], col = "indianred")
	# lines(GF25_revised[IC2,"CH4_mgm3"],col="blue")
	# # --> remove start und es scheint ein paar Sprünge zu haben...
	# abline(v=parse_date_time3("21.03.2021 15:00",tz="Etc/GMT-1"))
	# abline(v=parse_date_time3("23.03.2021 08:00",tz="Etc/GMT-1")) # Da hat es einen Sprung gegeben. Sieht man, wenn man es mit den anderen GFs vergleicht.
	GF25_revised["21.03.2021 to 21.03.2021 15:00",c("CH4_mgm3","CH4_ppm")] <- NA_real_


#############
### GF-26 ###
#############

	GF26_revised <- GF_26[which(
		GF_26$received_power > rp_gf26 & 
		GF_26$R2 > r2_thresh
		),][indMCs]
	
	## IC1
	# x11(width=12)
	# par(oma=c(0,0,0,2))
	# plot(GF_26[IC1,c("CH4_mgm3","received_power")],col2="lightgrey",col=NA,type2="l")
	# lines(GF_26[IC1,"CH4_mgm3"], col = "indianred")
	# lines(GF26_revised[IC1,"CH4_mgm3"], col="blue")
	# plot(GF_26[indMCs,"R2"],type="l",col="lightgrey")
	# plot(GF26_revised[,"CH4_mgm3"])
	# --> hat 2-3 Intervalle, die man rausnehmen könnte. Aber sind zu wenige und spielen so keine Rolle.
	# --> hat es Sprünge
	# lines(GF16_revised[,"CH4_mgm3"],col="orange")
	# --> hat es nicht
	## MC
	# x11(width=12)
	# par(oma=c(0,0,0,2))
	# plot(GF_26[MC,c("CH4_mgm3","received_power")],col2="lightgrey",col=NA,type2="l",ylim=c(1,2))
	# lines(GF_26[MC,"CH4_mgm3"], col = "indianred")
	# lines(GF26_revised[MC,"CH4_mgm3"],col="blue")
	# --> sieht gut aus. Aber am 21.03 ca. 2:50 Uhr könnte es einen Sprung haben
	# lines(GF16_revised[,"CH4_mgm3"],col="orange")
	# --> hat es nicht
	## IC2
	# x11(width=12)
	# par(oma=c(0,0,0,2))
	# plot(GF_26[IC2,c("CH4_mgm3","received_power")],col2="lightgrey",col=NA,type2="l",ylim=c(1,2))
	# lines(GF_26[IC2,"CH4_mgm3"], col = "indianred")
	# lines(GF26_revised[IC2,"CH4_mgm3"],col="blue")
	#  # # einiges ausschliessen. höherer Threshold ist aber nicht sinnvoll, da sonst zuviel wegfällt.
	# abline(v=parse_date_time3("21.03.2021 14:30",tz="Etc/GMT-1"))
	# abline(v=parse_date_time3(c("24.03.2021 07:00","24.03.2021 08:50"),tz="Etc/GMT-1"))
	# abline(v=parse_date_time3(c("25.03.2021 08:00","25.03.2021 08:40"),tz="Etc/GMT-1"))
	GF26_revised["21.03.2021 to 21.03.2021 14:30",c("CH4_mgm3","CH4_ppm")] <- NA_real_
	GF26_revised["24.03.2021 07:00 to 24.03.2021 08:50",c("CH4_mgm3","CH4_ppm")] <- NA_real_
	GF26_revised["25.03.2021 08:00 to 25.03.2021 08:40",c("CH4_mgm3","CH4_ppm")] <- NA_real_


#################
### pool data ###
#################

	Cov_thresh <- 0.75

	GF16_10min <- pool(GF16_revised,"10mins",st.to=StartMeas,et.to=StopMeas) %>cNA% Cov_thresh
	GF17_10min <- pool(GF17_revised,GF16_10min) %>cNA% Cov_thresh
	GF18_10min <- pool(GF18_revised,GF16_10min) %>cNA% Cov_thresh
	GF25_10min <- pool(GF25_revised,GF16_10min) %>cNA% Cov_thresh
	GF26_10min <- pool(GF26_revised,GF16_10min) %>cNA% Cov_thresh
	WS_10min <- pool(WS700[[1]]$data[,c("WS_0480","WD_corr")],GF16_10min)

	## add all concentrations toghether inclusive wind speed and wind direction
	GFs10min <- cbind(GF16=GF16_10min[,"CH4_mgm3"],GF17=GF17_10min[,"CH4_mgm3"],GF18=GF18_10min[,"CH4_mgm3"],GF25=GF25_10min[,"CH4_mgm3"],GF26=GF26_10min[,"CH4_mgm3"],
						WS=WS_10min[,"WS_0480"],WD=WS_10min[,"WD_corr"])
	
	GFs10min[,"MC"] <- NA_real_
	GFs10min[IC1,"MC"] <-"IC1"
	GFs10min[MC,"MC"] <-"MC"
	# GFs10min["19.03.2021 09:00 - 20.03.2021 09:00","MC"] <-NA_real_
	GFs10min[IC2,"MC"] <-"IC2"
	GFs10min$MC <- factor(GFs10min$MC, levels=c("QV1","IC1","MC","IC2"))

	## save, so it can be used for plots.
	# saveRDS(GFs10min,file=file.path(PathRSaves,"Conc_uncorr_10min.rds"))

	
################
### plotting ###
################

	x11(width=12)
	# png(file.path(PathFigures,"Conc_uncorr_IC1.png"),width=13,height=7,unit="in",res=300)
	plot(GFs10min[IC1,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min[IC1,1:5],na.rm=TRUE),ylab="CH4 [mg/m3]")
	lines(GFs10min[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min[,"GF26"],col=CH4Cols["GF26"])
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")
	title("GasFinder-Intercomparison (10-min averages uncorr) - IC1",outer=TRUE,line=-1)
	---> GF16 is sometimes a bit low??
	abline(v=parse_date_time3(c("07.03.2021 12:00","07.03.2021 18:00")))
	abline(v=parse_date_time3(c("08.03.2021 12:00","08.03.2021 18:00")))
	abline(v=parse_date_time3(c("09.03.2021 12:00","09.03.2021 18:00")))
	# dev.off()

	x11(width=12)
	# png(file.path(PathFigures,"Conc_uncorr_MC.png"),width=13,height=7,unit="in",res=300)
	plot(GFs10min[MC,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min[MC,1:5],na.rm=TRUE),ylab="CH4 [mg/m3]")
	lines(GFs10min[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min[,"GF26"],col=CH4Cols["GF26"])
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")
	title("GasFinder-Intercomparison (10-min averages uncorr) - MC",outer=TRUE,line=-1)
	--> between the releases GF25 a bit high and after the second release GF18 a bit high
	----> this could be outgasing of the methane that was still in the stable. after the first realease the wind direction was off so we do not see anything
	# dev.off()

	x11(width=12)
	# png(file.path(PathFigures,"Conc_uncorr_IC2.png"),width=13,height=7,unit="in",res=300)
	plot(GFs10min[IC2,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min[IC2,1:5],na.rm=TRUE),ylab="CH4 [mg/m3]")
	lines(GFs10min[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min[,"GF26"],col=CH4Cols["GF26"])
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")
	title("GasFinder-Intercomparison (10-min averages uncorr) - IC2",outer=TRUE,line=-1)
	abline(v=parse_date_time3("23.03.2021 11:00",tz="Etc/GMT-1"))
	--> not really consistent concentration differences. GF25 had probably a drop
	---> outgasing not possible to see as there was no upwind sensor
	# dev.off()

graphics.off()


################################################################################
### Check, if GF26 or GF25 has a concentration step/jump on 23/03/2021 08:00 ###
################################################################################

# only use data with a certain amount of wind speed
iWS <- which(GFs10min[,"WS"] > 0.5)


plot(GF16 ~ GF26, GFs10min[IC2][iWS], stats = "deming")
points(GF16 ~ GF26, GFs10min["23.03.2021 08:00 to "][iWS], col="red")
plot(GF16 ~ GF25, GFs10min[IC2][iWS], stats = "deming")
points(GF16 ~ GF25, GFs10min["23.03.2021 08:00 to "][iWS], col="red")

plot(GF17 ~ GF26, GFs10min[IC2][iWS], stats = "deming")
points(GF17 ~ GF26, GFs10min["23.03.2021 08:00 to "][iWS], col="red")
plot(GF17 ~ GF25, GFs10min[-MC], stats = "deming")
points(GF17 ~ GF25, GFs10min["23.03.2021 08:00 to "][iWS], col="red")

plot(GF18 ~ GF26, GFs10min[IC2][iWS], stats = "deming")
points(GF18 ~ GF26, GFs10min["23.03.2021 08:00 to "][iWS], col="red")
plot(GF18 ~ GF25, GFs10min[IC2][iWS], stats = "deming")
points(GF18 ~ GF25, GFs10min["23.03.2021 08:00 to "][iWS], col="red")

plot(GF25 ~ GF26, GFs10min[IC2][iWS], stats = "deming")
points(GF25 ~ GF26, GFs10min["23.03.2021 08:00 to "][iWS], col="red")
--> GF25 ~ GF26 seem really off. Otherwise it is rather GF26..?


## deltaC plotten
dt_GF <- cbind(as.data.table(GFs10min[]),st=st(GFs10min[]))

dC16 <- dt_GF[MC %in% c("IC1","MC","IC2"),{
	ggplot(.SD, aes(x=st)) +
	geom_smooth(aes(y=GF17-GF16),col=CH4Cols["GF17"]) +
	geom_line(aes(y=GF17-GF16),alpha=0.5,col=CH4Cols["GF17"]) +
	geom_smooth(aes(y=GF18-GF16),col=CH4Cols["GF18"]) +
	geom_line(aes(y=GF18-GF16),alpha=0.5,col=CH4Cols["GF18"]) +
	geom_smooth(aes(y=GF25-GF16),col=CH4Cols["GF25"]) +
	geom_line(aes(y=GF25-GF16),alpha=0.5,col=CH4Cols["GF25"]) +
	geom_smooth(aes(y=GF26-GF16),col=CH4Cols["GF26"]) +
	geom_line(aes(y=GF26-GF16),alpha=0.5,col=CH4Cols["GF26"]) +
	ylab("GFX - GF16 [mg/m3]") +
	xlab(NULL) +
	ylim(-0.1,0.2) +
	facet_grid(~ MC, scales="free_x", space="free_x") +
	theme_bw(base_size = 20) +
    theme(legend.position = "none", text = element_text(size=30))
}]

dC17 <- dt_GF[MC %in% c("IC1","MC","IC2"),{
	ggplot(.SD, aes(x=st)) +
	geom_smooth(aes(y=GF16-GF17),col=CH4Cols["GF16"]) +
	geom_point(aes(y=GF16-GF17),alpha=0.5,col=CH4Cols["GF16"]) +
	geom_smooth(aes(y=GF18-GF17),col=CH4Cols["GF18"]) +
	geom_point(aes(y=GF18-GF17),alpha=0.5,col=CH4Cols["GF18"]) +
	geom_smooth(aes(y=GF25-GF17),col=CH4Cols["GF25"]) +
	geom_point(aes(y=GF25-GF17),alpha=0.5,col=CH4Cols["GF25"]) +
	geom_smooth(aes(y=GF26-GF17),col=CH4Cols["GF26"]) +
	geom_point(aes(y=GF26-GF17),alpha=0.5,col=CH4Cols["GF26"]) +
	ylab("GFX - GF17") +
	xlab(NULL) + 
	ylim(-0.1,0.2) +
	facet_grid(~ MC, scales="free_x", space="free_x") +
	theme_bw(base_size = 20) +
    theme(legend.position = "none", text = element_text(size=30))
}]

dC18 <- dt_GF[MC %in% c("IC1","MC","IC2"),{
	ggplot(.SD, aes(x=st)) +
	geom_smooth(aes(y=GF17-GF18),col=CH4Cols["GF17"]) +
	geom_line(aes(y=GF17-GF18),alpha=0.5, col=CH4Cols["GF17"]) +
	geom_smooth(aes(y=GF16-GF18),col=CH4Cols["GF16"]) +
	geom_line(aes(y=GF16-GF18),alpha=0.5, col=CH4Cols["GF16"]) +
	geom_smooth(aes(y=GF25-GF18),col=CH4Cols["GF25"]) +
	geom_line(aes(y=GF25-GF18),alpha=0.5, col=CH4Cols["GF25"]) +
	geom_smooth(aes(y=GF26-GF18),col=CH4Cols["GF26"]) +
	geom_line(aes(y=GF26-GF18),alpha=0.5, col=CH4Cols["GF26"]) +
	ylab("GFX - GF18") +
	xlab(NULL) +
	ylim(-0.1,0.2) +
	facet_grid(~ MC, scales="free_x", space="free_x") +
	theme_bw(base_size = 20) +
    theme(legend.position = "none", text = element_text(size=30))
}]

dC25 <- dt_GF[MC %in% c("IC1","MC","IC2"),{
	ggplot(.SD, aes(x=st)) +
	geom_point(aes(y=GF16-GF25),alpha=0.5, col=CH4Cols["GF16"]) +
	geom_smooth(aes(y=GF16-GF25),col=CH4Cols["GF16"]) +
	geom_point(aes(y=GF17-GF25),alpha=0.5, col=CH4Cols["GF17"]) +
	geom_smooth(aes(y=GF17-GF25),col=CH4Cols["GF17"]) +
	geom_point(aes(y=GF18-GF25),alpha=0.5, col=CH4Cols["GF18"]) +
	geom_smooth(aes(y=GF18-GF25),col=CH4Cols["GF18"]) +
	geom_point(aes(y=GF26-GF25),alpha=0.5, col=CH4Cols["GF26"]) +
	geom_smooth(aes(y=GF26-GF25),col=CH4Cols["GF26"]) +
	ylab("GFX - GF25") +
	ylim(-0.2,0.1) +
	facet_grid(~ MC, scales="free_x", space="free_x") +
	theme_bw(base_size = 20) +
    theme(legend.position = "none", text = element_text(size=30))
}]

dC26 <- dt_GF[MC %in% c("IC1","MC","IC2"),{
	ggplot(.SD, aes(x=st)) +
	geom_point(aes(y=GF16-GF26),alpha=0.5, col=CH4Cols["GF16"]) +
	geom_smooth(aes(y=GF16-GF26),col=CH4Cols["GF16"]) +
	geom_point(aes(y=GF17-GF26),alpha=0.5, col=CH4Cols["GF17"]) +
	geom_smooth(aes(y=GF17-GF26),col=CH4Cols["GF17"]) +
	geom_point(aes(y=GF18-GF26),alpha=0.5, col=CH4Cols["GF18"]) +
	geom_smooth(aes(y=GF18-GF26),col=CH4Cols["GF18"]) +
	geom_point(aes(y=GF25-GF26),alpha=0.5, col=CH4Cols["GF25"]) +
	geom_smooth(aes(y=GF25-GF26),col=CH4Cols["GF25"]) +
	ylab("GFX - GF26 [mg/m3]") +
	xlab(NULL) +
	ylim(-0.2,0.1) +
	facet_grid(~ MC, scales="free_x", space="free_x") +
	theme_bw(base_size = 20) +
    theme(legend.position = "none", text = element_text(size=30))
}]

dCGF_uncorr <- ggarrange(dC18,dC26,ncol=1,align="hv")
dCGF_uncorr
# ggsave(file.path(PathFigures,"dCGFs_uncorr.png"), width=12*1.5,height=12)


#### only concentrations
# dt_GF[MC %in% c("IC1","MC","IC2"),{
Conc_uncorr <-	dt_GF[MC != "QV1",{
	ggplot(.SD, aes(x=st)) +
	geom_line(aes(y=GF16),col=CH4Cols["GF16"], lwd=1) +
	geom_line(aes(y=GF17),col=CH4Cols["GF17"], lwd=1) +
	geom_line(aes(y=GF18),col=CH4Cols["GF18"], lwd=1) +
	geom_line(aes(y=GF25),col=CH4Cols["GF25"], lwd=1) +
	geom_line(aes(y=GF26),col=CH4Cols["GF26"], lwd=1) +
	ylab("CH4 Conc [mg / m3]") +
	xlab(NULL) + 
	facet_grid(~ MC, scales="free_x", space="free_x") +
	theme_bw(base_size = 20) +
    theme(legend.position = "none", text = element_text(size=30))
}]

Conc_uncorr
# ggsave(file.path(PathFigures,"Conc_uncorr.png"),width=12*2.5, height=12)


GF_uncorr <- ggarrange(Conc_uncorr,dC26,ncol=1,align="hv")
GF_uncorr
# ggsave(file.path(PathFigures,"GF_uncorr.png"),GF_uncorr,width=12*1.5, height=12)


## GF18 is the best. All the others have some kind of patterns/structures in the records, which one can not filter out.
## -> if I got it right, GF26 had an offset and not GF25
## but look at plots, once the offset and span are corrected.

x11()
par(mfrow=c(2,2))
	mod2616 <- plot(GF26 ~ GF16, GFs10min[-MC][iWS], stats = "deming",main="IC1 & IC2")
	mod2616_end <- plot(GF26 ~ GF16, GFs10min[-MC][-"23.03.2021 to 27.03.2021"][iWS], stats = "deming",main="IC1 & IC2 & -End")
	mod2516 <- plot(GF25 ~ GF16, GFs10min[-MC][iWS], stats = "deming",main="IC1 & IC2")
	mod2516_end <- plot(GF25 ~ GF16, GFs10min[-MC][-"23.03.2021 to 27.03.2021"][iWS], stats = "deming",main="IC1 & IC2 & -End")
--> GF25 better

x11()
par(mfrow=c(2,2))
	mod2617 <- plot(GF26 ~ GF17, GFs10min[-MC][iWS], stats = "deming",main="IC1 & IC2")
	mod2617_end <- plot(GF26 ~ GF17, GFs10min[-MC][-"23.03.2021 to 27.03.2021"][iWS], stats = "deming",main="IC1 & IC2 & -End")
	mod2517 <- plot(GF25 ~ GF17, GFs10min[-MC][iWS], stats = "deming",main="IC1 & IC2")
	mod2517_end <- plot(GF25 ~ GF17, GFs10min[-MC][-"23.03.2021 to 27.03.2021"][iWS], stats = "deming",main="IC1 & IC2 & -End")
--> similar

x11()
par(mfrow=c(2,2))
	mod2618 <- plot(GF26 ~ GF18, GFs10min[-MC][iWS], stats = "deming",main="IC1 & IC2")
	mod2618_end <- plot(GF26 ~ GF18, GFs10min[-MC][-"23.03.2021 to 27.03.2021"][iWS], stats = "deming",main="IC1 & IC2 & -End")
	mod2518 <- plot(GF25 ~ GF18, GFs10min[-MC][iWS], stats = "deming",main="IC1 & IC2")
	mod2518_end <- plot(GF25 ~ GF18, GFs10min[-MC][-"23.03.2021 to 27.03.2021"][iWS], stats = "deming",main="IC1 & IC2 & -End")
--> GF25 better

x11()
par(mfrow=c(2,2))
	mod2625 <- plot(GF26 ~ GF25, GFs10min[-MC][iWS], stats = "deming",main="IC1 & IC2")
	mod2625_end <- plot(GF26 ~ GF25, GFs10min[-MC][-"23.03.2021 to 27.03.2021"][iWS], stats = "deming",main="IC1 & IC2 & -End")

# ----> as I don't know if GF25 or GF26 had a step change, I match GF18 on GF25 (without the period at the end) and then everyting to GF18 


##################################
##################################
#####                        #####
#####    Intercalibration    #####
#####                        #####
##################################
##################################

# only use data with more than 0.5 m/s wind speed
iWS <- which(GFs10min[,"WS"] > 0.5)

# match GF18 on GF25 without the period at the end
par(mfrow=c(3,2))
	mod2518_IC1 <- plot(GF25 ~ GF18, GFs10min[IC1][iWS], stats = "deming",main="IC1")
	mod2518_IC2 <- plot(GF25 ~ GF18, GFs10min[IC2][iWS], stats = "deming",main="IC2")
	# mod2518 <- plot(GF25 ~ GF18, GFs10min[-MC][iWS], stats = "deming",main="-MC")
	mod2518_IC12 <- plot(GF25 ~ GF18, GFs10min[-MC][iWS], stats = "deming",main="IC1 & IC2")
	mod2518_IC12E <- plot(GF25 ~ GF18, GFs10min[-MC][-"23.03.2021 to 27.03.2021"][-"21.03.2021 18:30 to 22.03.2021 01:50"][iWS], stats = "deming",main="IC1 & IC2 -Start and -End")
	mod2518_IC2E <- plot(GF25 ~ GF18, GFs10min[IC2][-"21.03.2021 18:30 to 22.03.2021 01:50"][-"23.03.2021 to 27.03.2021"][iWS], stats = "deming",main="IC2 & -End")


ggplot(GFs10min[iWS],aes(x=GF18,y=GF25)) +
	geom_point(data=GFs10min[IC1],col=paste0(VerCols["IC1"],50)) +
	geom_abline(intercept=coef(mod2518_IC1)[1],slope=coef(mod2518_IC1)[2],lwd=1.6,col=VerCols["IC1"]) +
	annotate(geom="text", x=1.3, y=2.2, label=paste0("IC1: GF25 = ",round(coef(mod2518_IC1)[1],2)," + ",round(coef(mod2518_IC1)[2],2)," x GF18"), col=VerCols["IC1"]) +
	geom_point(data=GFs10min[IC2],col=paste0(VerCols["IC2"],50)) +
	geom_abline(intercept=coef(mod2518_IC2)[1],slope=coef(mod2518_IC2)[2],lwd=1.6,col=VerCols["IC2"]) +
	annotate(geom="text", x=1.3, y=2.15, label=paste0("IC2: GF25 = ",round(coef(mod2518_IC2)[1],2)," + ",round(coef(mod2518_IC2)[2],2)," x GF18"), col=VerCols["IC2"]) +
	geom_point(data=GFs10min[-MC],col=paste0(VerCols["IC12"],50)) +
	geom_abline(intercept=coef(mod2518_IC12)[1],slope=coef(mod2518_IC12)[2],lwd=1.6,col=VerCols["IC12"]) +
	annotate(geom="text", x=1.3, y=2.1, label=paste0("IC12: GF25 = ",round(coef(mod2518_IC12)[1],2)," + ",round(coef(mod2518_IC12)[2],2)," x GF18"), col=VerCols["IC12"]) +
	geom_point(data=GFs10min[-MC][-"23.03.2021 to 27.03.2021"][-"21.03.2021 18:30 to 22.03.2021 01:50"],col=paste0(VerCols["IC12E"],50)) +
	geom_abline(intercept=coef(mod2518_IC12E)[1],slope=coef(mod2518_IC12E)[2],lwd=1.6,col=VerCols["IC12E"]) +
	annotate(geom="text", x=1.3, y=2.05, label=paste0("IC12E: GF25 = ",round(coef(mod2518_IC12E)[1],2)," + ",round(coef(mod2518_IC12E)[2],2)," x GF18"), col=VerCols["IC12E"]) +
	geom_point(data=GFs10min[IC2][-"23.03.2021 to 27.03.2021"][-"21.03.2021 18:30 to 22.03.2021 01:50"],col=paste0(VerCols["IC2E"],50)) +
	geom_abline(intercept=coef(mod2518_IC2E)[1],slope=coef(mod2518_IC2E)[2],lwd=1.6,col=VerCols["IC2E"]) +
	annotate(geom="text", x=1.3, y=2.0, label=paste0("IC2E: GF25 = ",round(coef(mod2518_IC2E)[1],2)," + ",round(coef(mod2518_IC2E)[2],2)," x GF18"), col=VerCols["IC2E"]) +
	theme_bw(base_size = 20) +
    theme(legend.position = c(1,1), text = element_text(size=30))

par(mfrow=c(2,2))
	mod2518_2 <- plot(MASS::rlm(GF25 ~ GF18, GFs10min[IC1][iWS], stats = "deming"))
	title(main="IC1",outer=TRUE,line=-2)
par(mfrow=c(2,2))	
	mod2518_3 <- plot(MASS::rlm(GF25 ~ GF18, GFs10min[IC2][iWS], stats = "deming",main="IC2"))
	title(main="IC2",outer=TRUE,line=-2)
par(mfrow=c(2,2))
	mod2518_mMC <- plot(MASS::rlm(GF25 ~ GF18, GFs10min[-MC][iWS], stats = "deming",main="-MC"))
	title(main="-MC",outer=TRUE,line=-2)
par(mfrow=c(2,2))
	mod2518_4 <- plot(MASS::rlm(GF25 ~ GF18, GFs10min[-MC][-"23.03.2021 to 27.03.2021"][iWS], stats = "deming"))
	title(main="IC1 & IC2 -End",outer=TRUE,line=-2)
par(mfrow=c(2,2))
	mod2518_4 <- plot(MASS::rlm(GF25 ~ GF18, GFs10min[IC2][-"23.03.2021 to 27.03.2021"][iWS], stats = "deming"))
	title(main="IC1 & IC2 -End",outer=TRUE,line=-2)
	

# match GF18 to GF26 without the period at the end
par(mfrow=c(3,2))
	mod2618_IC1 <- plot(GF26 ~ GF18, GFs10min[IC1][iWS], stats = "deming",main="IC1")
	mod2618_IC2 <- plot(GF26 ~ GF18, GFs10min[IC2][iWS], stats = "deming",main="IC2")
	# mod2618 <- plot(GF26 ~ GF18, GFs10min[-MC][iWS], stats = "deming",main="-MC")
	mod2618_IC12 <- plot(GF26 ~ GF18, GFs10min[-MC][iWS], stats = "deming",main="IC1 & IC2")
	mod2618_IC12E <- plot(GF26 ~ GF18, GFs10min[-MC][-"23.03.2021 to 27.03.2021"][-"21.03.2021 18:30 to 22.03.2021 01:50"][iWS], stats = "deming",main="IC1 & IC2 -End")
	mod2618_IC2E <- plot(GF26 ~ GF18, GFs10min[IC2][-"23.03.2021 to 27.03.2021"][-"21.03.2021 18:30 to 22.03.2021 01:50"][iWS], stats = "deming",main="IC2 & -End")

ggplot(GFs10min[iWS],aes(x=GF18,y=GF26)) +
	geom_point(data=GFs10min[IC1],col=paste0(VerCols["IC1"],50)) +
	geom_abline(intercept=coef(mod2618_IC1)[1],slope=coef(mod2618_IC1)[2],lwd=1.6,col=VerCols["IC1"]) +
	annotate(geom="text", x=1.3, y=2.2, label=paste0("IC1: GF26 = ",round(coef(mod2618_IC1)[1],2)," + ",round(coef(mod2618_IC1)[2],2)," x GF18"), col=VerCols["IC1"]) +
	geom_point(data=GFs10min[IC2],col=paste0(VerCols["IC2"],50)) +
	geom_abline(intercept=coef(mod2618_IC2)[1],slope=coef(mod2618_IC2)[2],lwd=1.6,col=VerCols["IC2"]) +
	annotate(geom="text", x=1.3, y=2.15, label=paste0("IC2: GF26 = ",round(coef(mod2618_IC2)[1],2)," + ",round(coef(mod2618_IC2)[2],2)," x GF18"), col=VerCols["IC2"]) +
	geom_point(data=GFs10min[-MC],col=paste0(VerCols["IC12"],50)) +
	geom_abline(intercept=coef(mod2618_IC12)[1],slope=coef(mod2618_IC12)[2],lwd=1.6,col=VerCols["IC12"]) +
	annotate(geom="text", x=1.3, y=2.1, label=paste0("IC12: GF26 = ",round(coef(mod2618_IC12)[1],2)," + ",round(coef(mod2618_IC12)[2],2)," x GF18"), col=VerCols["IC12"]) +
	geom_point(data=GFs10min[-MC][-"23.03.2021 to 27.03.2021"][-"21.03.2021 18:30 to 22.03.2021 01:50"],col=paste0(VerCols["IC12E"],50)) +
	geom_abline(intercept=coef(mod2618_IC12E)[1],slope=coef(mod2618_IC12E)[2],lwd=1.6,col=VerCols["IC12E"]) +
	annotate(geom="text", x=1.3, y=2.05, label=paste0("IC12E: GF26 = ",round(coef(mod2618_IC12E)[1],2)," + ",round(coef(mod2618_IC12E)[2],2)," x GF18"), col=VerCols["IC12E"]) +
	geom_point(data=GFs10min[IC2][-"23.03.2021 to 27.03.2021"][-"21.03.2021 18:30 to 22.03.2021 01:50"],col=paste0(VerCols["IC2E"],50)) +
	geom_abline(intercept=coef(mod2618_IC2E)[1],slope=coef(mod2618_IC2E)[2],lwd=1.6,col=VerCols["IC2E"]) +
	annotate(geom="text", x=1.3, y=2.0, label=paste0("IC2E: GF26 = ",round(coef(mod2618_IC2E)[1],2)," + ",round(coef(mod2618_IC2E)[2],2)," x GF18"), col=VerCols["IC2E"]) +
	theme_bw(base_size = 20) +
    theme(legend.position = c(1,1), text = element_text(size=30))


par(mfrow=c(2,2))
	mod2618_2 <- plot(MASS::rlm(GF26 ~ GF18, GFs10min[IC1][iWS], stats = "deming"))
	title(main="IC1",outer=TRUE,line=-2)
par(mfrow=c(2,2))	
	mod2618_3 <- plot(MASS::rlm(GF26 ~ GF18, GFs10min[IC2][iWS], stats = "deming",main="IC2"))
	title(main="IC2",outer=TRUE,line=-2)
par(mfrow=c(2,2))
	mod2618_mMC <- plot(MASS::rlm(GF26 ~ GF18, GFs10min[-MC][iWS], stats = "deming",main="-MC"))
	title(main="-MC",outer=TRUE,line=-2)
par(mfrow=c(2,2))
	mod2618_4 <- plot(MASS::rlm(GF26 ~ GF18, GFs10min[-MC][-"23.03.2021 to 27.03.2021"][-"21.03.2021 18:30 to 22.03.2021 01:50"][iWS], stats = "deming"))
	title(main="IC1 & IC2 -End",outer=TRUE,line=-2)
par(mfrow=c(2,2))
	mod2618_4 <- plot(MASS::rlm(GF26 ~ GF18, GFs10min[IC2][-"23.03.2021 to 27.03.2021"][-"21.03.2021 18:30 to 22.03.2021 01:50"][iWS], stats = "deming"))
	title(main="IC1 & IC2 -End",outer=TRUE,line=-2)
	


	mod2518 <- plot(GF25 ~ GF18, GFs10min[-MC][-"23.03.2021 to 27.03.2021"][-"21.03.2021 18:30 to 22.03.2021 01:50"][iWS], stats = "deming",main="IC2 & -End")
	# mod2518 <- plot(GF25 ~ GF18, GFs10min[IC2][-"23.03.2021 to 27.03.2021"][iWS], stats = "deming",main="IC2 & -End")
	cfs2518 <- coef(mod2518)

# apply calibration
	GFs10min_revised_18 <- GFs10min
	GFs10min_revised_18[,"GF18"] <- cfs2518[1] + cfs2518[2] * GFs10min[,"GF18"]

### match everything with the corrected GF18
	# mod1816 <- plot(GF18 ~ GF16, GFs10min_revised_18[-MC][iWS], stats = "deming")
	mod1816 <- plot(GF18 ~ GF16, GFs10min_revised_18[-MC][-"21.03.2021 18:30 to 22.03.2021 01:50"][iWS], stats = "deming")
	# mod1816 <- plot(GF18 ~ GF16, GFs10min_revised_18[IC1][iWS], stats = "deming")
	# mod1816 <- plot(GF18 ~ GF16, GFs10min_revised_18[IC2][iWS], stats = "deming")
	# mod1816 <- plot(GF18 ~ GF16, GFs10min_revised_18[IC2][-"23.03.2021 to 27.03.2021"][iWS], stats = "deming")
	# mod1816 <- plot(GF18 ~ GF16, GFs10min_revised_18[-MC][-"23.03.2021 to 27.03.2021"][iWS], stats = "deming")
	cfs1816 <- coef(mod1816)
	abline(0,1,col="red",lwd=2)

	# mod1817 <- plot(GF18 ~ GF17, GFs10min_revised_18[-MC][iWS], stats = "deming") ### die ganz tiefen Werte könnten nicht korrekt sein
	mod1817 <- plot(GF18 ~ GF17, GFs10min_revised_18[-MC][-"21.03.2021 18:30 to 22.03.2021 01:50"][iWS], stats = "deming") ### die ganz tiefen Werte könnten nicht korrekt sein
	# mod1817 <- plot(GF18 ~ GF17, GFs10min_revised_18[IC1][iWS], stats = "deming")
	# mod1817 <- plot(GF18 ~ GF17, GFs10min_revised_18[IC2][iWS], stats = "deming")
	cfs1817 <- coef(mod1817)
	abline(0,1,col="red",lwd=2)

	# mod1825 <- plot(GF18 ~ GF25, GFs10min_revised_18[-MC][iWS], stats = "deming")
	# mod1825 <- plot(GF18 ~ GF25, GFs10min_revised_18[IC1][iWS], stats = "deming")
	# cfs1825 <- coef(mod1825)
	# abline(0,1,col="red",lwd=2)

	# mod1826 <- plot(GF18 ~ GF26, GFs10min_revised_18[-MC][iWS], stats = "deming")
	# mod1826 <- plot(GF18 ~ GF26, GFs10min_revised_18[IC1][iWS], stats = "deming")
	# mod1826 <- plot(GF18 ~ GF26, GFs10min_revised_18[IC2][iWS], stats = "deming")
	# mod1826 <- plot(GF18 ~ GF26, GFs10min_revised_18[-MC][-"23.03.2021 to 27.03.2021"][iWS], stats = "deming") # zur Sicherheit das End ausschliessen
	mod1826 <- plot(GF18 ~ GF26, GFs10min_revised_18[-MC][-"23.03.2021 to 27.03.2021"][-"21.03.2021 18:30 to 22.03.2021 01:50"][iWS], stats = "deming") # zur Sicherheit das End ausschliessen
	cfs1826 <- coef(mod1826)
	abline(0,1,col="red",lwd=2)


	cfsList <- list(
		GF16 = cfs1816
		,GF17 = cfs1817
		,GF18 = c(0,1) # already corrected
		,GF25 = c(0,1) # als reference for GF18
		,GF26 = cfs1826
		)

	
########################
### apply correction ###
########################
	
	GFs10min_revised <- GFs10min_revised_18
	GFs10min_revised[,1:5] <- mapply(function(x,y){
		y[1] + y[2]*x
	}, x = GFs10min_revised_18[,1:5], y = cfsList[names(GFs10min_revised_18[,1:5])], SIMPLIFY = FALSE)


#################
### plot data ###
#################

	graphics.off()
	
	x11(width=14, height = 6)
	# png(file.path(PathFigures,"Conc_corr1_IC1.png"),width=13,height=7,unit="in",res=300)
	plot(GFs10min_revised[IC1,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min_revised[IC1,c('GF16','GF17','GF18','GF25','GF26')],na.rm=TRUE),ylab="CH4 [mg/m3]")
	lines(GFs10min_revised[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min_revised[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min_revised[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min_revised[,"GF26"],col=CH4Cols["GF26"])
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")
	title("GasFinder-Intercomparison (10-min averages, calibrated) - IC1",outer=TRUE,line=-1)
	---> in the beginning concentrations are really off
	abline(v=parse_date_time3(c("05.03.2021","06.03.2021 05:00"),tz="Etc/GMT-1"))
	abline(v=parse_date_time3(c("08.03.2021 12:00","08.03.2021 18:00"),tz="Etc/GMT-1"))
	---> GF26 often lower then the others. With low concentrations GF16 is lower
	# dev.off()

	x11(width=14, height = 6)
	# png(file.path(PathFigures,"Conc_corr1_MC.png"),width=13,height=7,unit="in",res=300)
	plot(GFs10min_revised[MC,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min_revised[MC,c('GF16','GF17','GF18','GF25','GF26')],na.rm=TRUE),ylab="CH4 [mg/m3]")
	lines(GFs10min_revised[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min_revised[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min_revised[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min_revised[,"GF26"],col=CH4Cols["GF26"])
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")
	title("GasFinder-Intercomparison (10-min averages, calibrated) - MC",outer=TRUE,line=-1)
	--> GF25 rather high. Between the releases and after the second there is some offset
	# dev.off()

	x11(width=14, height = 6)
	# png(file.path(PathFigures,"Conc_corr1_IC2.png"),width=13,height=7,unit="in",res=300)
	plot(GFs10min_revised[IC2,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min_revised[IC2,c('GF16','GF17','GF18','GF25','GF26')],na.rm=TRUE),ylab="CH4 [mg/m3]")
	lines(GFs10min_revised[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min_revised[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min_revised[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min_revised[,"GF26"],col=CH4Cols["GF26"])
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")
	title("GasFinder-Intercomparison (10-min averages, calibrated) - IC2",outer=TRUE,line=-1)
	---> at the beginning GF17 too low and after the 24.03.2021 12:00 GF26 seems beeing heigher 
	---> after the offset and span correction it does not seem that there is outgasing. But GF25 is quite high
	# dev.off()


	---> There still seems to be an offset between the concentrations. Further, GF25 seems to be a little high during the experiment
		 (the dC should decrease more in my opinion). The offset is larger between releases and after the release than before. Could be due to the power failure.
		 With IC2, GF25 is suddenly lower. But I still use it as an adjustment because of the span.

dt_GF_rev <- cbind(as.data.table(GFs10min_revised[]),st=st(GFs10min_revised[]))

dt_GF_rev[MC %in% c("IC1","MC","IC2"),{
	ggplot(.SD, aes(x=st)) +
	geom_smooth(aes(y=GF17-GF16),col=CH4Cols["GF17"]) +
	geom_point(aes(y=GF17-GF16),col=CH4Cols["GF17"]) +
	geom_smooth(aes(y=GF18-GF16),col=CH4Cols["GF18"]) +
	geom_point(aes(y=GF18-GF16),col=CH4Cols["GF18"]) +
	geom_smooth(aes(y=GF25-GF16),col=CH4Cols["GF25"]) +
	geom_point(aes(y=GF25-GF16),col=CH4Cols["GF25"]) +
	geom_smooth(aes(y=GF26-GF16),col=CH4Cols["GF26"]) +
	geom_point(aes(y=GF26-GF16),col=CH4Cols["GF26"]) +
	ylab("GFX - GF16") +
	ylim(-0.05,0.05) +
	facet_grid(~ MC, scales="free_x", space="free_x") +
	theme_bw(base_size = 20) +
    theme(legend.position = "none", text = element_text(size=30))
}]

dt_GF_rev[MC %in% c("IC1","MC","IC2"),{
	ggplot(.SD, aes(x=st)) +
	geom_smooth(aes(y=GF16-GF17),col=CH4Cols["GF16"]) +
	geom_point(aes(y=GF16-GF17),col=CH4Cols["GF16"]) +
	geom_smooth(aes(y=GF18-GF17),col=CH4Cols["GF18"]) +
	geom_point(aes(y=GF18-GF17),col=CH4Cols["GF18"]) +
	geom_smooth(aes(y=GF25-GF17),col=CH4Cols["GF25"]) +
	geom_point(aes(y=GF25-GF17),col=CH4Cols["GF25"]) +
	geom_smooth(aes(y=GF26-GF17),col=CH4Cols["GF26"]) +
	geom_point(aes(y=GF26-GF17),col=CH4Cols["GF26"]) +
	ylab("GFX - GF17") +
	ylim(-0.05,0.05) +
	facet_grid(~ MC, scales="free_x", space="free_x") +
	theme_bw(base_size = 20) +
    theme(legend.position = "none", text = element_text(size=30))
}]

dt_GF_rev[MC %in% c("IC1","MC","IC2"),{
	ggplot(.SD, aes(x=st)) +
	geom_smooth(aes(y=GF17-GF18),col=CH4Cols["GF17"]) +
	geom_point(aes(y=GF17-GF18),col=CH4Cols["GF17"]) +
	geom_smooth(aes(y=GF16-GF18),col=CH4Cols["GF16"]) +
	geom_point(aes(y=GF16-GF18),col=CH4Cols["GF16"]) +
	geom_smooth(aes(y=GF25-GF18),col=CH4Cols["GF25"]) +
	geom_point(aes(y=GF25-GF18),col=CH4Cols["GF25"]) +
	geom_smooth(aes(y=GF26-GF18),col=CH4Cols["GF26"]) +
	geom_point(aes(y=GF26-GF18),col=CH4Cols["GF26"]) +
	ylab("GFX - GF18") +
	ylim(-0.05,0.05) +
	facet_grid(~ MC, scales="free_x", space="free_x") +
	theme_bw(base_size = 20) +
    theme(legend.position = "none", text = element_text(size=30))
}]

dt_GF_rev[MC %in% c("IC1","MC","IC2"),{
	ggplot(.SD, aes(x=st)) +
	geom_point(aes(y=GF16-GF25),col=CH4Cols["GF16"]) +
	geom_smooth(aes(y=GF16-GF25),col=CH4Cols["GF16"]) +
	geom_point(aes(y=GF17-GF25),col=CH4Cols["GF17"]) +
	geom_smooth(aes(y=GF17-GF25),col=CH4Cols["GF17"]) +
	geom_point(aes(y=GF18-GF25),col=CH4Cols["GF18"]) +
	geom_smooth(aes(y=GF18-GF25),col=CH4Cols["GF18"]) +
	geom_point(aes(y=GF26-GF25),col=CH4Cols["GF26"]) +
	geom_smooth(aes(y=GF26-GF25),col=CH4Cols["GF26"]) +
	ylab("GFX - GF25") +
	ylim(-0.05,0.05) +
	facet_grid(~ MC, scales="free_x", space="free_x") +
	theme_bw(base_size = 20) +
    theme(legend.position = "none", text = element_text(size=30))
}]

dt_GF_rev[MC %in% c("IC1","MC","IC2"),{
	ggplot(.SD, aes(x=st)) +
	geom_point(aes(y=GF16-GF26),col=CH4Cols["GF16"]) +
	geom_smooth(aes(y=GF16-GF26),col=CH4Cols["GF16"]) +
	geom_point(aes(y=GF17-GF26),col=CH4Cols["GF17"]) +
	geom_smooth(aes(y=GF17-GF26),col=CH4Cols["GF17"]) +
	geom_point(aes(y=GF18-GF26),col=CH4Cols["GF18"]) +
	geom_smooth(aes(y=GF18-GF26),col=CH4Cols["GF18"]) +
	geom_point(aes(y=GF25-GF26),col=CH4Cols["GF25"]) +
	geom_smooth(aes(y=GF25-GF26),col=CH4Cols["GF25"]) +
	ylab("GFX - GF26") +
	ylim(-0.05,0.05) +
	facet_grid(~ MC, scales="free_x", space="free_x") +
	theme_bw(base_size = 20) +
    theme(legend.position = "none", text = element_text(size=30))
}]


	graphics.off()


---> an offset correction with the data before (and after) the release might be needed. But first make background concentration in IC2 (saving issues).


#########################################
### Interpolate concentrations in IC2 ###
#########################################

GFs10min_revised[,"GFall"] <- NA_real_
GFs10min_revised[IC2,"GFall"] <- rowMeans(GFs10min_revised[IC2,1:5],na.rm=TRUE)

	plot(GFs10min_revised[IC2,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min_revised[IC2,1:5],na.rm=TRUE),ylab="CH4 [mg/m3]")
	lines(GFs10min_revised[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min_revised[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min_revised[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min_revised[,"GF26"],col=CH4Cols["GF26"])
	lines(GFs10min_revised[,"GFall"],col="black",lwd=2)
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")
	title("GasFinder-Intercomparison (10-min averages, calibrated) - IC1",outer=TRUE,line=-1)
	abline(v=parse_date_time3(c("22.03.2021 09:20","22.03.2021 14:20"),tz="Etc/GMT-1"),col="orange")
	abline(v=parse_date_time3(c("23.03.2021 09:20","23.03.2021 14:20"),tz="Etc/GMT-1"),col="orange",lty=2)

---> there might be outgasing but it is impossible to say that with the data so I just use the time right after the concentration drop

dC_t <- (GFs10min_revised["22.03.2021 14:20"]$GFall - GFs10min_revised["22.03.2021 09:20"]$GFall)/length(GFs10min_revised["22.03.2021 09:20 - 22.03.2021 14:20"]$GFall)
A_t <- GFs10min_revised["22.03.2021 09:20"]$GFall
GFs10min_revised["22.03.2021 09:20 - 22.03.2021 14:20","GFall"] <- A_t + dC_t * seq(1:30)

	plot(GFs10min_revised[IC2,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min_revised[IC2,1:5],na.rm=TRUE),ylab="CH4 [mg/m3]")
	lines(GFs10min_revised[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min_revised[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min_revised[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min_revised[,"GF26"],col=CH4Cols["GF26"])
	lines(GFs10min_revised[,"GFall"],col="black",lwd=2)
	lines(GFs10min_revised["22.03.2021 09:20 - 22.03.2021 14:20","GFall"],col="magenta",lwd=2)

	plot(GFs10min_revised["22.03.2021 - 23.03.2021","GF16"],col=CH4Cols["GF16"],ylim=c(1.2,1.5),ylab="CH4 [mg/m3]")
	lines(GFs10min_revised[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min_revised[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min_revised[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min_revised[,"GF26"],col=CH4Cols["GF26"])
	lines(GFs10min_revised[,"GFall"],col="black",lwd=2)
	lines(GFs10min_revised["22.03.2021 09:20 - 22.03.2021 14:20","GFall"],col="magenta",lwd=2)

graphics.off()


######################################
### Correction of offset during MC ###
######################################

### use data before, between the releases and a few hours after to correct all GasFinders with a average offset to the background

## define times that should be used as reference for correction
	plot(GFs10min_revised[MC,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min_revised[MC,1:5],na.rm=TRUE),ylab="CH4 [mg/m3]")
	lines(GFs10min_revised[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min_revised[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min_revised[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min_revised[,"GF26"],col=CH4Cols["GF26"])
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")
	title("GasFinder-Intercomparison (10-min averages, calibrated) - MC",outer=TRUE,line=-1)
	abline(v=parse_date_time3(c("19.03.2021 07:00","19.03.2021 10:00"),tz="Etc/GMT-1"))
	abline(v=parse_date_time3(c("19.03.2021 17:10","19.03.2021 21:30"),tz="Etc/GMT-1"),col="red")
	abline(v=parse_date_time3(c("20.03.2021 07:30","20.03.2021 12:00"),tz="Etc/GMT-1"),col="blue")
	abline(v=parse_date_time3(c("20.03.2021 18:00","20.03.2021 23:00"),tz="Etc/GMT-1"),col="magenta")
	abline(v=parse_date_time3(c("20.03.2021 18:00","21.03.2021 12:00"),tz="Etc/GMT-1"),col="orange",lty=3)

P1 <- c("19.03.2021 08:50 - 19.03.2021 10:00")
P2 <- c("19.03.2021 17:10 - 19.03.2021 21:30")
P3 <- c("20.03.2021 07:30 - 20.03.2021 11:50")
P4 <- c("20.03.2021 18:00 - 20.03.2021 23:00")
P5 <- c("20.03.2021 18:00 - 21.03.2021 12:00")
P6 <- c("20.03.2021 08:00 - 21.03.2021 12:00")
P123 <- c(P1,P2,P3)
P23 <- c(P2,P3)
P14 <- c(P1,P4)

# wind speed
iWS <- which(GFs10min[,"WS"] > 0.5)

	# GF16
	mod2616P1 <- plot(I(GF26 - GF16) ~ 1, GFs10min_revised[P1][iWS], stats = "rlm")
	mod2616P2 <- plot(I(GF26 - GF16) ~ 1, GFs10min_revised[P2][iWS], stats = "rlm")
	mod2616P3 <- plot(I(GF26 - GF16) ~ 1, GFs10min_revised[P3][iWS], stats = "rlm")
	mod2616P4 <- plot(I(GF26 - GF16) ~ 1, GFs10min_revised[P4][iWS], stats = "rlm")
	mod2616P123 <- plot(I(GF26 - GF16) ~ 1, GFs10min_revised[P123][iWS], stats = "rlm")
	mod2616P23 <- plot(I(GF26 - GF16) ~ 1, GFs10min_revised[P23][iWS], stats = "rlm")
	mod2616P14 <- plot(I(GF26 - GF16) ~ 1, GFs10min_revised[P14][iWS], stats = "rlm")
	mod2616P5 <- plot(I(GF26 - GF16) ~ 1, GFs10min_revised[P5][iWS], stats = "rlm")
	mod2616P6 <- plot(I(GF26 - GF16) ~ 1, GFs10min_revised[P6][iWS], stats = "rlm")
	cfs2616P1 <- coef(mod2616P1)
	cfs2616P2 <- coef(mod2616P2)
	cfs2616P3 <- coef(mod2616P3)
	cfs2616P4 <- coef(mod2616P4)
	cfs2616P123 <- coef(mod2616P123)
	cfs2616P23 <- coef(mod2616P23)
	cfs2616P14 <- coef(mod2616P14)
	cfs2616P5 <- coef(mod2616P5)
	cfs2616P6 <- coef(mod2616P6)
	cfs16 <- data.frame(GF16 = c(P1 = cfs2616P1, P2 = cfs2616P2, P3 = cfs2616P3, P4 = cfs2616P4,
				 P123 = cfs2616P123, P23 = cfs2616P23, P14 = cfs2616P14, P5 = cfs2616P5, P6 = cfs2616P6))

	# GF17
	mod2617P1 <- plot(I(GF26 - GF17) ~ 1, GFs10min_revised[P1][iWS], stats = "rlm")
	mod2617P2 <- plot(I(GF26 - GF17) ~ 1, GFs10min_revised[P2][iWS], stats = "rlm")
	mod2617P3 <- plot(I(GF26 - GF17) ~ 1, GFs10min_revised[P3][iWS], stats = "rlm")
	mod2617P4 <- plot(I(GF26 - GF17) ~ 1, GFs10min_revised[P4][iWS], stats = "rlm")
	mod2617P123 <- plot(I(GF26 - GF17) ~ 1, GFs10min_revised[P123][iWS], stats = "rlm")
	mod2617P23 <- plot(I(GF26 - GF17) ~ 1, GFs10min_revised[P23][iWS], stats = "rlm")
	mod2617P14 <- plot(I(GF26 - GF17) ~ 1, GFs10min_revised[P14][iWS], stats = "rlm")
	mod2617P5 <- plot(I(GF26 - GF17) ~ 1, GFs10min_revised[P5][iWS], stats = "rlm")
	mod2617P6 <- plot(I(GF26 - GF17) ~ 1, GFs10min_revised[P6][iWS], stats = "rlm")
	cfs2617P1 <- coef(mod2617P1)
	cfs2617P2 <- coef(mod2617P2)
	cfs2617P3 <- coef(mod2617P3)
	cfs2617P4 <- coef(mod2617P4)
	cfs2617P123 <- coef(mod2617P123)
	cfs2617P23 <- coef(mod2617P23)
	cfs2617P14 <- coef(mod2617P14)
	cfs2617P5 <- coef(mod2617P5)
	cfs2617P6 <- coef(mod2617P6)
	cfs17 <- data.frame(GF17 = c(P1 = cfs2617P1, P2 = cfs2617P2, P3 = cfs2617P3, P4 = cfs2617P4,
				 P123 = cfs2617P123, P23 = cfs2617P23, P14 = cfs2617P14, P5 = cfs2617P5, P6 = cfs2617P6))

	# GF18
	mod2618P1 <- plot(I(GF26 - GF18) ~ 1, GFs10min_revised[P1][iWS], stats = "rlm")
	mod2618P2 <- plot(I(GF26 - GF18) ~ 1, GFs10min_revised[P2][iWS], stats = "rlm")
	mod2618P3 <- plot(I(GF26 - GF18) ~ 1, GFs10min_revised[P3][iWS], stats = "rlm")
	mod2618P4 <- plot(I(GF26 - GF18) ~ 1, GFs10min_revised[P4][iWS], stats = "rlm")
	mod2618P123 <- plot(I(GF26 - GF18) ~ 1, GFs10min_revised[P123][iWS], stats = "rlm")
	mod2618P23 <- plot(I(GF26 - GF18) ~ 1, GFs10min_revised[P23][iWS], stats = "rlm")
	mod2618P14 <- plot(I(GF26 - GF18) ~ 1, GFs10min_revised[P14][iWS], stats = "rlm")
	mod2618P5 <- plot(I(GF26 - GF18) ~ 1, GFs10min_revised[P5][iWS], stats = "rlm")
	mod2618P6 <- plot(I(GF26 - GF18) ~ 1, GFs10min_revised[P6][iWS], stats = "rlm")
	cfs2618P1 <- coef(mod2618P1)
	cfs2618P2 <- coef(mod2618P2)
	cfs2618P3 <- coef(mod2618P3)
	cfs2618P4 <- coef(mod2618P4)
	cfs2618P123 <- coef(mod2618P123)
	cfs2618P23 <- coef(mod2618P23)
	cfs2618P14 <- coef(mod2618P14)
	cfs2618P5 <- coef(mod2618P5)
	cfs2618P6 <- coef(mod2618P6)
	cfs18 <- data.frame(GF18 = c(P1 = cfs2618P1, P2 = cfs2618P2, P3 = cfs2618P3, P4 = cfs2618P4,
				 P123 = cfs2618P123, P23 = cfs2618P23, P14 = cfs2618P14, P5 = cfs2618P5, P6 = cfs2618P6))

	# GF25
	mod2625P1 <- plot(I(GF26 - GF25) ~ 1, GFs10min_revised[P1][iWS], stats = "rlm")
	mod2625P2 <- plot(I(GF26 - GF25) ~ 1, GFs10min_revised[P2][iWS], stats = "rlm")
	mod2625P3 <- plot(I(GF26 - GF25) ~ 1, GFs10min_revised[P3][iWS], stats = "rlm")
	mod2625P4 <- plot(I(GF26 - GF25) ~ 1, GFs10min_revised[P4][iWS], stats = "rlm")
	mod2625P123 <- plot(I(GF26 - GF25) ~ 1, GFs10min_revised[P123][iWS], stats = "rlm")
	mod2625P23 <- plot(I(GF26 - GF25) ~ 1, GFs10min_revised[P23][iWS], stats = "rlm")
	mod2625P14 <- plot(I(GF26 - GF25) ~ 1, GFs10min_revised[P14][iWS], stats = "rlm")
	mod2625P5 <- plot(I(GF26 - GF25) ~ 1, GFs10min_revised[P5][iWS], stats = "rlm")
	mod2625P6 <- plot(I(GF26 - GF25) ~ 1, GFs10min_revised[P6][iWS], stats = "rlm")
	cfs2625P1 <- coef(mod2625P1)
	cfs2625P2 <- coef(mod2625P2)
	cfs2625P3 <- coef(mod2625P3)
	cfs2625P4 <- coef(mod2625P4)
	cfs2625P123 <- coef(mod2625P123)
	cfs2625P23 <- coef(mod2625P23)
	cfs2625P14 <- coef(mod2625P14)
	cfs2625P5 <- coef(mod2625P5)
	cfs2625P6 <- coef(mod2625P6)
	cfs25 <- data.frame(GF25 = c(P1 = cfs2625P1, P2 = cfs2625P2, P3 = cfs2625P3, P4 = cfs2625P4,
				 P123 = cfs2625P123, P23 = cfs2625P23, P14 = cfs2625P14, P5 = cfs2625P5, P6 = cfs2625P6))

	cfsall <- cbind(cfs16,cfs17,cfs18,cfs25)	
	# with GF25 the differences between the version are largest. Probably us P14
	
	cfsListP1 <- list(GF16 = cfs2616P1, GF17 = cfs2617P1, GF18 = cfs2618P1, GF25 = cfs2625P1, GF26 = c(0,1), GFall =c(0,1))
	cfsListP2 <- list(GF16 = cfs2616P2, GF17 = cfs2617P2, GF18 = cfs2618P2, GF25 = cfs2625P2, GF26 = c(0,1), GFall =c(0,1))
	cfsListP3 <- list(GF16 = cfs2616P3, GF17 = cfs2617P3, GF18 = cfs2618P3, GF25 = cfs2625P3, GF26 = c(0,1), GFall =c(0,1))
	cfsListP4 <- list(GF16 = cfs2616P4, GF17 = cfs2617P4, GF18 = cfs2618P4, GF25 = cfs2625P4, GF26 = c(0,1), GFall =c(0,1))
	cfsListP123 <- list(GF16 = cfs2616P123, GF17 = cfs2617P123, GF18 = cfs2618P123, GF25 = cfs2625P123, GF26 = c(0,1), GFall =c(0,1))
	cfsListP23 <- list(GF16 = cfs2616P23, GF17 = cfs2617P23, GF18 = cfs2618P23, GF25 = cfs2625P23, GF26 = c(0,1), GFall =c(0,1))
	cfsListP14 <- list(GF16 = cfs2616P14, GF17 = cfs2617P14, GF18 = cfs2618P14, GF25 = cfs2625P14, GF26 = c(0,1), GFall =c(0,1))
	cfsListP5 <- list(GF16 = cfs2616P5, GF17 = cfs2617P5, GF18 = cfs2618P5, GF25 = cfs2625P5, GF26 = c(0,1), GFall =c(0,1))
	cfsListP6 <- list(GF16 = cfs2616P6, GF17 = cfs2617P6, GF18 = cfs2618P6, GF25 = cfs2625P6, GF26 = c(0,1), GFall =c(0,1))
		

####################################
### apply correction/calibration ###
####################################

	GFs10min_revised_P1<-GFs10min_revised_P2<-GFs10min_revised_P3<-GFs10min_revised_P4<-GFs10min_revised_P123<-GFs10min_revised_P23<-GFs10min_revised_P14<-GFs10min_revised_P5<-GFs10min_revised_P6<-GFs10min_revised
	## apply only on MC
	GFs10min_revised_P1[MC,1:5] <- mapply(function(x,y){
			y[1] + x
		}, x = GFs10min_revised[MC,1:5], y = cfsListP1[names(GFs10min_revised[,1:5])], SIMPLIFY = FALSE)
	GFs10min_revised_P2[MC,1:5] <- mapply(function(x,y){
			y[1] + x
		}, x = GFs10min_revised[MC,1:5], y = cfsListP2[names(GFs10min_revised[,1:5])], SIMPLIFY = FALSE)
	GFs10min_revised_P3[MC,1:5] <- mapply(function(x,y){
			y[1] + x
		}, x = GFs10min_revised[MC,1:5], y = cfsListP3[names(GFs10min_revised[,1:5])], SIMPLIFY = FALSE)
	GFs10min_revised_P4[MC,1:5] <- mapply(function(x,y){
			y[1] + x
		}, x = GFs10min_revised[MC,1:5], y = cfsListP4[names(GFs10min_revised[,1:5])], SIMPLIFY = FALSE)
	GFs10min_revised_P123[MC,1:5] <- mapply(function(x,y){
			y[1] + x
		}, x = GFs10min_revised[MC,1:5], y = cfsListP123[names(GFs10min_revised[,1:5])], SIMPLIFY = FALSE)
	GFs10min_revised_P23[MC,1:5] <- mapply(function(x,y){
			y[1] + x
		}, x = GFs10min_revised[MC,1:5], y = cfsListP23[names(GFs10min_revised[,1:5])], SIMPLIFY = FALSE)
	GFs10min_revised_P14[MC,1:5] <- mapply(function(x,y){
			y[1] + x
		}, x = GFs10min_revised[MC,1:5], y = cfsListP14[names(GFs10min_revised[,1:5])], SIMPLIFY = FALSE)
	GFs10min_revised_P5[MC,1:5] <- mapply(function(x,y){
			y[1] + x
		}, x = GFs10min_revised[MC,1:5], y = cfsListP5[names(GFs10min_revised[,1:5])], SIMPLIFY = FALSE)
	GFs10min_revised_P6[MC,1:5] <- mapply(function(x,y){
			y[1] + x
		}, x = GFs10min_revised[MC,1:5], y = cfsListP6[names(GFs10min_revised[,1:5])], SIMPLIFY = FALSE)

	
	## P1
	plot(GFs10min_revised_P1[MC,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min_revised_P1[MC,1:5],na.rm=TRUE),ylab="CH4 [mg/m3]",main="P1")
	lines(GFs10min_revised_P1[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min_revised_P1[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min_revised_P1[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min_revised_P1[,"GF26"],col=CH4Cols["GF26"])
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")

	## P2
	plot(GFs10min_revised_P2[MC,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min_revised_P2[MC,1:5],na.rm=TRUE),ylab="CH4 [mg/m3]",main="P2")
	lines(GFs10min_revised_P2[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min_revised_P2[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min_revised_P2[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min_revised_P2[,"GF26"],col=CH4Cols["GF26"])
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")
	
	## P3
	plot(GFs10min_revised_P3[MC,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min_revised_P3[MC,1:5],na.rm=TRUE),ylab="CH4 [mg/m3]",main="P3")
	lines(GFs10min_revised_P3[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min_revised_P3[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min_revised_P3[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min_revised_P3[,"GF26"],col=CH4Cols["GF26"])
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")
		
	## P4
	plot(GFs10min_revised_P4[MC,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min_revised_P4[MC,1:5],na.rm=TRUE),ylab="CH4 [mg/m3]",main="P4")
	lines(GFs10min_revised_P4[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min_revised_P4[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min_revised_P4[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min_revised_P4[,"GF26"],col=CH4Cols["GF26"])
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")
	
	## P123
	plot(GFs10min_revised_P123[MC,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min_revised_P123[MC,1:5],na.rm=TRUE),ylab="CH4 [mg/m3]",main="P123")
	lines(GFs10min_revised_P123[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min_revised_P123[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min_revised_P123[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min_revised_P123[,"GF26"],col=CH4Cols["GF26"])
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")
	
	## P23
	plot(GFs10min_revised_P23[MC,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min_revised_P23[MC,1:5],na.rm=TRUE),ylab="CH4 [mg/m3]",main="P23")
	lines(GFs10min_revised_P23[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min_revised_P23[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min_revised_P23[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min_revised_P23[,"GF26"],col=CH4Cols["GF26"])
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")
	
	## P14
	plot(GFs10min_revised_P14[MC,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min_revised_P14[MC,1:5],na.rm=TRUE),ylab="CH4 [mg/m3]",main="P14")
	lines(GFs10min_revised_P14[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min_revised_P14[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min_revised_P14[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min_revised_P14[,"GF26"],col=CH4Cols["GF26"])
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")
		
	## P5
	plot(GFs10min_revised_P5[MC,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min_revised_P5[MC,1:5],na.rm=TRUE),ylab="CH4 [mg/m3]",main="P5")
	lines(GFs10min_revised_P5[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min_revised_P5[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min_revised_P5[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min_revised_P5[,"GF26"],col=CH4Cols["GF26"])
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")
		
	## P6
	plot(GFs10min_revised_P6[MC,"GF16"],col=CH4Cols["GF16"],ylim=range(GFs10min_revised_P6[MC,1:5],na.rm=TRUE),ylab="CH4 [mg/m3]",main="P6")
	lines(GFs10min_revised_P6[,"GF17"],col=CH4Cols["GF17"])
	lines(GFs10min_revised_P6[,"GF18"],col=CH4Cols["GF18"])
	lines(GFs10min_revised_P6[,"GF25"],col=CH4Cols["GF25"])
	lines(GFs10min_revised_P6[,"GF26"],col=CH4Cols["GF26"])
	legend("topleft",legend=names(CH4Cols),col=CH4Cols,lty=1,bty="n")
	
---> in terms of corrected concentrations P23 seems most reasonable. However, this might include outgasing. On the otherhand P1, P4 or P14 lead to
--> higher background concentrations between the releases. This could be actually the case due to unfavourable wind directions but right before
--> the second release there is still this negative offset.


	GFs10min_revised_MC <- GFs10min_revised_P123
	# GFs10min_revised_MC <- GFs10min_revised_P23
	# GFs10min_revised_MC <- GFs10min_revised_P6
	
	# GFs10min_revised_MC <- GFs10min_revised # without further correction


dt_GF_rev_MC <- cbind(as.data.table(GFs10min_revised_MC[]),st=st(GFs10min_revised_MC[]))

dt_GF_rev_MC[MC %in% c("IC1","MC","IC2"),{
	ggplot(.SD, aes(x=st)) +
	geom_smooth(aes(y=GF17-GF16),col=CH4Cols["GF17"]) +
	geom_point(aes(y=GF17-GF16),col=CH4Cols["GF17"]) +
	geom_smooth(aes(y=GF18-GF16),col=CH4Cols["GF18"]) +
	geom_point(aes(y=GF18-GF16),col=CH4Cols["GF18"]) +
	geom_smooth(aes(y=GF25-GF16),col=CH4Cols["GF25"]) +
	geom_point(aes(y=GF25-GF16),col=CH4Cols["GF25"]) +
	geom_smooth(aes(y=GF26-GF16),col=CH4Cols["GF26"]) +
	geom_point(aes(y=GF26-GF16),col=CH4Cols["GF26"]) +
	ylab("GFX - GF16") +
	ylim(-0.05,0.05) +
	facet_grid(~ MC, scales="free_x", space="free_x") +
	theme_bw(base_size = 20) +
    theme(legend.position = "none", text = element_text(size=30))
}]

dt_GF_rev_MC[MC %in% c("IC1","MC","IC2"),{
	ggplot(.SD, aes(x=st)) +
	geom_smooth(aes(y=GF16-GF17),col=CH4Cols["GF16"]) +
	geom_point(aes(y=GF16-GF17),col=CH4Cols["GF16"]) +
	geom_smooth(aes(y=GF18-GF17),col=CH4Cols["GF18"]) +
	geom_point(aes(y=GF18-GF17),col=CH4Cols["GF18"]) +
	geom_smooth(aes(y=GF25-GF17),col=CH4Cols["GF25"]) +
	geom_point(aes(y=GF25-GF17),col=CH4Cols["GF25"]) +
	geom_smooth(aes(y=GF26-GF17),col=CH4Cols["GF26"]) +
	geom_point(aes(y=GF26-GF17),col=CH4Cols["GF26"]) +
	ylab("GFX - GF17") +
	ylim(-0.05,0.05) +
	facet_grid(~ MC, scales="free_x", space="free_x") +
	theme_bw(base_size = 20) +
    theme(legend.position = "none", text = element_text(size=30))
}]

dt_GF_rev_MC[MC %in% c("IC1","MC","IC2"),{
	ggplot(.SD, aes(x=st)) +
	geom_smooth(aes(y=GF17-GF18),col=CH4Cols["GF17"]) +
	geom_point(aes(y=GF17-GF18),col=CH4Cols["GF17"]) +
	geom_smooth(aes(y=GF16-GF18),col=CH4Cols["GF16"]) +
	geom_point(aes(y=GF16-GF18),col=CH4Cols["GF16"]) +
	geom_smooth(aes(y=GF25-GF18),col=CH4Cols["GF25"]) +
	geom_point(aes(y=GF25-GF18),col=CH4Cols["GF25"]) +
	geom_smooth(aes(y=GF26-GF18),col=CH4Cols["GF26"]) +
	geom_point(aes(y=GF26-GF18),col=CH4Cols["GF26"]) +
	ylab("GFX - GF18") +
	ylim(-0.05,0.05) +
	facet_grid(~ MC, scales="free_x", space="free_x") +
	theme_bw(base_size = 20) +
    theme(legend.position = "none", text = element_text(size=30))
}]

dt_GF_rev_MC[MC %in% c("IC1","MC","IC2"),{
	ggplot(.SD, aes(x=st)) +
	geom_point(aes(y=GF16-GF25),col=CH4Cols["GF16"]) +
	geom_smooth(aes(y=GF16-GF25),col=CH4Cols["GF16"]) +
	geom_point(aes(y=GF17-GF25),col=CH4Cols["GF17"]) +
	geom_smooth(aes(y=GF17-GF25),col=CH4Cols["GF17"]) +
	geom_point(aes(y=GF18-GF25),col=CH4Cols["GF18"]) +
	geom_smooth(aes(y=GF18-GF25),col=CH4Cols["GF18"]) +
	geom_point(aes(y=GF26-GF25),col=CH4Cols["GF26"]) +
	geom_smooth(aes(y=GF26-GF25),col=CH4Cols["GF26"]) +
	ylab("GFX - GF25") +
	ylim(-0.05,0.05) +
	facet_grid(~ MC, scales="free_x", space="free_x") +
	theme_bw(base_size = 20) +
    theme(legend.position = "none", text = element_text(size=30))
}]

dt_GF_rev_MC[MC %in% c("IC1","MC","IC2"),{
	ggplot(.SD, aes(x=st)) +
	geom_point(aes(y=GF16-GF26),col=CH4Cols["GF16"]) +
	geom_smooth(aes(y=GF16-GF26),col=CH4Cols["GF16"]) +
	geom_point(aes(y=GF17-GF26),col=CH4Cols["GF17"]) +
	geom_smooth(aes(y=GF17-GF26),col=CH4Cols["GF17"]) +
	geom_point(aes(y=GF18-GF26),col=CH4Cols["GF18"]) +
	geom_smooth(aes(y=GF18-GF26),col=CH4Cols["GF18"]) +
	geom_point(aes(y=GF25-GF26),col=CH4Cols["GF25"]) +
	geom_smooth(aes(y=GF25-GF26),col=CH4Cols["GF25"]) +
	ylab("GFX - GF26") +
	ylim(-0.05,0.05) +
	facet_grid(~ MC, scales="free_x", space="free_x") +
	theme_bw(base_size = 20) +
    theme(legend.position = "none", text = element_text(size=30))
}]

graphics.off()


#################
### save data ###
#################

saveRDS(GFs10min_revised[,paste0("GF",c(16:18,25,26,"all"))],file=file.path(PathRSaves,"Conc_wo_10min.rds")) # without additional offset correction
saveRDS(GFs10min_revised_P23[,paste0("GF",c(16:18,25,26,"all"))],file=file.path(PathRSaves,"Conc_P23_10min.rds")) # with the P23 correction
saveRDS(GFs10min_revised_P6[,paste0("GF",c(16:18,25,26,"all"))],file=file.path(PathRSaves,"Conc_P6_10min.rds")) # with the P6 correction
saveRDS(GFs10min_revised_P4[,paste0("GF",c(16:18,25,26,"all"))],file=file.path(PathRSaves,"Conc_P4_10min.rds")) # with the P4 correction

