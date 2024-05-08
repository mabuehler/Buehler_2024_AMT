

###############################################
###############################################
#####                                     #####
#####    Datenfilterung der Emissionen    #####
#####                                     #####
###############################################
###############################################


#################################
### Header (Pfade, Libraries) ###
#################################

	library(ibts)
	library(RgoogleMaps)
	library(bLSmodelR)
	library(RColorBrewer)
	library(ggplot2)

PathData <- "Path to /data"		
PathRSaves <- "Path to /RSaves"

	PfadFigures <- file.path(dirname(PfadRSaves),"Figures")
	source("~/repos/3_Scripts/gel-scripts/wgs84-ch1903.r")
	Cat.Path <- paste0(PfadDaten,"/Catalogs")

###################################
### Funktionen und Definitionen ###
###################################

lines_sec2xy <- function(xyMK,sensor,node=1,wd,col="lightblue",lwd=2,...){
	GF <- xyMK[xyMK[,1] %in% sensor,]
	sens <- as.numeric(GF[GF[,3] == node,4:5])
	b <- tan((90 - wd)/180*pi)
	x <- if(wd <= 180) 600 else -600
	y <- sens[2] - (sens[1] - x)*b
	lines(c(sens[1],x),c(sens[2],y),col=col,lwd=lwd,...)
}

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


	### MKs
	QV1 <- " to 28.01.2021 12:00"
	QV2 <- "05.03.2021 to 10.03.2021 18:00"
	QV1u2 <- " to 10.03.2021 18:00"
	MK <- "18.03.2021 11:00 - 21.03.2021 14:00"
	MK_rel <- "19.03.2021 10:00 to 20.03.2021 08:00"
	QV3	<- "21.03.2021 14:00 to "
	# QV3_rel	<- 

#######################
### Laden der Daten ###
#######################

	# Pfadlängen & Geometrie
	load(file.path(PfadRSaves,"Geometry_STO.RData"))
	# Google map
	STO_Map <- ReadMapTile(file.path(PfadFigures,"STO_GoogleMaps.png"))
	# Emissions Resultate
	Emiss_Result_raw <- readRDS(file = file.path(PfadRSaves, "STO_Emiss_orig_P23_10min.rds"))
	# Emiss_Result_raw <- readRDS(file = file.path(PfadRSaves, "STO_Emiss_orig_P6_10min.rds"))
	# Emiss_Result_raw <- readRDS(file = file.path(PfadRSaves, "STO_Emiss_short_P6_10min.rds"))
	# Emiss_Result_raw <- readRDS(file = file.path(PfadRSaves, "STO_Emiss_orig_P23_30min.rds"))

#####################################################
### die verschiedenen Emissionsruns zusammenfügen ###
#####################################################

## nicht nötig für diese Version

# vollständiger Zeitdatansatz generieren (wegen plotten) 18.03.2021 12:50 - 22.03.2021 14:20
time_vect <- data.table(st=seq(0,585)*600 + parse_date_time3("18.03.2021 12:50",tz="Etc/GMT-1"))
Emiss_Result <- merge(time_vect,Emiss_Result_raw,by="st",all=TRUE)
# MFC für Darstellung auf 0 setzen vor und nach Release
# Emiss_Result[Campaign=="MK" & is.na(Q_MFC),]

setkey(Emiss_Result,st)

##########################
### create XY Geometry ###
##########################

Sensors_MK_xy <- CH.to.map(STO_Map,Sensors_MK)
Sensors_QV3_xy <- CH.to.map(STO_Map,Sensors_QV3)
Sources_xy <- CH.to.map(STO_Map,Sources)

################################
### mal ne übersicht plotten ###
################################

plot(Q_GF16 ~ st, Emiss_Result[Campaign %in% c("MK")],type="l",col=CH4Cols["GF16"],ylim=c(-10,10))
lines(Q_GF17 ~ st, Emiss_Result[Campaign %in% c("MK")],col=CH4Cols["GF17"])
lines(Q_GF18 ~ st, Emiss_Result[Campaign %in% c("MK")],col=CH4Cols["GF18"])
lines(Q_GF25 ~ st, Emiss_Result[Campaign %in% c("MK")],col=CH4Cols["GF25"])
lines(Q_GF26 ~ st, Emiss_Result[Campaign %in% c("MK")],col=CH4Cols["GF26"])
lines(Q_MFC ~ st, Emiss_Result[Campaign %in% c("MK")],col="black",lwd=2)

####################################
### Nur Daten mit Bise verwenden ###
####################################

Result_raw <- Emiss_Result[Wind_dir == "NE" | is.na(Wind_dir)]


################################################
### add missing times caused by power outage ###
################################################

t_step <- as.numeric(median(Result_raw[,et-st],na.rm=TRUE))
st_fill <- parse_date_time3("2021-03-19 09:00",tz='Etc/GMT-1')
st_dt <- data.table(st=st_fill + seq(0,length.out=250/t_step,by=t_step) * 60, et=st_fill + seq(0,length.out=250/t_step,by=t_step) * 60 + t_step*60)


Result_ls <- lapply(Result_raw[,unique(Sonic)], function(x) {
	# browser()
	out <- merge(Result_raw[Sonic==x],st_dt,by=c('st','et'),all=TRUE)
	out[,Sonic := x]
	out[st %in% st_dt[1],Campaign := 'MK']
	out[,Source := 'Schopf']
})

Result <- rbindlist(Result_ls)






#####################################
### Windrichtungsfilterung für MK ###
#####################################

## Jenachdem verwende ich diese gar nicht. mal schauen, wird sicher ziemlich ein rekursiver Prozess. 

plot(Sources,Sensors_MK)
# Südwest und Südostpunkte des Schopfs bestimmen
# Schopf Ecke NW
ENW <- as.numeric(Sources[Sources[,3] == max(Sources[,3]),2:3])
# Schopf Ecke NE
ENE <- as.numeric(Sources[Sources[,2] == max(Sources[,2]),2:3])
# Schopf Ecke SW
ESW <- as.numeric(Sources[Sources[,2] == min(Sources[,2]),2:3])
# Schopf Ecke SE
ESE <- as.numeric(Sources[Sources[,3] == min(Sources[,3]),2:3])
## Mitte des GasFinderpfads bestimmen
P16M <- as.numeric(colMeans(Sensors_MK[Sensors_MK[,1] == "GF16" ,c(4,5)]))
P17M <- as.numeric(colMeans(Sensors_MK[Sensors_MK[,1] == "GF17" ,c(4,5)]))
P18M <- as.numeric(colMeans(Sensors_MK[Sensors_MK[,1] == "GF18" ,c(4,5)]))
P25M <- as.numeric(colMeans(Sensors_MK[Sensors_MK[,1] == "GF25" ,c(4,5)]))

## dummy Sensor in der Mitte des Pfades erstellen zur Überprüfung der Winkel
PMitte <- genSensors(data.frame(Name=c("GF16","GF17","GF18","GF25")
	,x = c(P16M[1],P17M[1],P18M[1],P25M[1])
	,y = c(P16M[2],P17M[2],P18M[2],P25M[2])
	,z=rep(1,4),rep(NA,4),rep(1,4)
	))

Mitte_xy <- CH.to.map(STO_Map,PMitte)

### Winkel für Variante Ecken (wir sonst über Karte) bestimmen für Messkampagne (MK)
iWD16 <- c(90 - round(atan(abs(Sensors_MK[Sensors_MK[,1] == "GF16" & Sensors_MK[,3] == 1,5]-ENW[2])/
		abs(Sensors_MK[Sensors_MK[,1] == "GF16" & Sensors_MK[,3] == 1,4]-ENW[1]))*180/pi,0),
	0 + round(atan(abs(Sensors_MK[Sensors_MK[,1] == "GF16" & Sensors_MK[,3] == 2,4]-ESE[1])/
		abs(Sensors_MK[Sensors_MK[,1] == "GF16" & Sensors_MK[,3] == 2,5]-ESE[2]))*180/pi,0))
iWD17 <- c(360 - round(atan(abs(Sensors_MK[Sensors_MK[,1] == "GF17" & Sensors_MK[,3] == 1,4]-ESW[1])/
		abs(Sensors_MK[Sensors_MK[,1] == "GF17" & Sensors_MK[,3] == 1,5]-ESW[2]))*180/pi,0),
	90 + round(atan(abs(Sensors_MK[Sensors_MK[,1] == "GF17" & Sensors_MK[,3] == 2,5]-ESE[2])/
		abs(Sensors_MK[Sensors_MK[,1] == "GF17" & Sensors_MK[,3] == 2,4]-ESE[1]))*180/pi,0))
iWD18 <- c(90 - round(atan(abs(Sensors_MK[Sensors_MK[,1] == "GF18" & Sensors_MK[,3] == 1,5]-ESW[2])/
		abs(Sensors_MK[Sensors_MK[,1] == "GF18" & Sensors_MK[,3] == 1,4]-ESW[1]))*180/pi,0),
	0 + round(atan(abs(Sensors_MK[Sensors_MK[,1] == "GF18" & Sensors_MK[,3] == 2,4]-ESE[1])/
		abs(Sensors_MK[Sensors_MK[,1] == "GF18" & Sensors_MK[,3] == 2,5]-ESE[2]))*180/pi,0))
iWD25 <- c(90 - round(atan(abs(Sensors_MK[Sensors_MK[,1] == "GF25" & Sensors_MK[,3] == 1,5]-ENW[2])/
		abs(Sensors_MK[Sensors_MK[,1] == "GF25" & Sensors_MK[,3] == 1,4]-ENW[1]))*180/pi,0),
	0 + round(atan(abs(Sensors_MK[Sensors_MK[,1] == "GF25" & Sensors_MK[,3] == 2,4]-ESE[1])/
		abs(Sensors_MK[Sensors_MK[,1] == "GF25" & Sensors_MK[,3] == 2,5]-ESE[2]))*180/pi,0))

### Winkel für Variante Mitte bestimmen MK
iWD16M <- c(0 + round(atan(abs(P16M[1]-ENW[1])/abs(P16M[2]-ENW[2]))*180/pi,0),
	90 - round(atan(abs(P16M[2]-ESE[2])/abs(P16M[1]-ESE[1]))*180/pi,0))
iWD17M <- c(0 + round(atan(abs(P17M[1]-ENW[1])/abs(P17M[2]-ENW[2]))*180/pi,0),
	90 - round(atan(abs(P17M[2]-ESE[2])/abs(P17M[1]-ESE[1]))*180/pi,0))
iWD18M <- c(0 + round(atan(abs(P18M[1]-ENW[1])/abs(P18M[2]-ENW[2]))*180/pi,0),
	90 - round(atan(abs(P18M[2]-ESE[2])/abs(P18M[1]-ESE[1]))*180/pi,0))
iWD25M <- c(0 + round(atan(abs(P25M[1]-ENW[1])/abs(P25M[2]-ENW[2]))*180/pi,0),
	90 - round(atan(abs(P25M[2]-ESE[2])/abs(P25M[1]-ESE[1]))*180/pi,0))


######################################
### Windrichtungsfilterung für QV3 ###
######################################

## Jenachdem verwende ich diese gar nicht. mal schauen, wird sicher ziemlich ein rekursiver Prozess. 

plot(Sources,Sensors_QV3)
## Daten nur für GF18 rechnen und dann auf alle anwenden, da QV und Unterschiede wahrscheinlich sehr gering. GF18 ist in der Mitte

# Südwest und Südostpunkte des Schopfs bestimmen
# Schopf Ecke NW
ENW <- as.numeric(Sources[Sources[,3] == max(Sources[,3]),2:3])
# Schopf Ecke NE
ENE <- as.numeric(Sources[Sources[,2] == max(Sources[,2]),2:3])
# Schopf Ecke SW
ESW <- as.numeric(Sources[Sources[,2] == min(Sources[,2]),2:3])
# Schopf Ecke SE
ESE <- as.numeric(Sources[Sources[,3] == min(Sources[,3]),2:3])

## Mitte des GasFinderpfads bestimmen
PMQV3 <- as.numeric(colMeans(Sensors_QV3[Sensors_QV3[,1] == "GF18" ,c(4,5)]))
## dummy Sensor in der Mitte des Pfades erstellen zur Überprüfung der Winkel
PMitteQV3 <- genSensors(data.frame(Name=c("GFM")
	,x = PMQV3[1]
	,y = PMQV3[2]
	,z=1,NA,1
	))

MitteQV3_xy <- CH.to.map(STO_Map,PMitteQV3)

### Winkel für Variante Ecken (wir sonst über Karte) bestimmen für Messkampagne (MK)

iWDQV3 <- c(360 - round(atan(abs(Sensors_QV3[Sensors_QV3[,1] == "GF18" & Sensors_QV3[,3] == 2,4]-ESW[1])/
		abs(Sensors_QV3[Sensors_QV3[,1] == "GF18" & Sensors_QV3[,3] == 2,5]-ESW[2]))*180/pi,0),
	90 + round(atan(abs(Sensors_QV3[Sensors_QV3[,1] == "GF18" & Sensors_QV3[,3] == 1,5]-ESE[2])/
		abs(Sensors_QV3[Sensors_QV3[,1] == "GF18" & Sensors_QV3[,3] == 1,4]-ESE[1]))*180/pi,0))

### Winkel für Variante Mitte bestimmen MK
iWDMQV3 <- c(0 + round(atan(abs(PMQV3[1]-ENW[1])/abs(PMQV3[2]-ENW[2]))*180/pi,0),
	90 - round(atan(abs(PMQV3[2]-ESE[2])/abs(PMQV3[1]-ESE[1]))*180/pi,0))


###################################
### Winkel auf Karte überprüfen ###
###################################
graphics.off()

## GF17
# iWD17neu <- c(22,68)
PlotOnStaticMap(STO_Map)
plot(Sensors_MK_xy[Sensors_MK_xy[,1] %in% c("GF17"),],sensor.text.args=list(labels="",cex=2),lines.args=list(lwd=2,col="orange",lty=3),points.args = list(pch = 20, cex = 2, col ="orange"),add=TRUE)
lines_sec2xy(Sensors_MK_xy,"GF17",1,iWD17[1],lty=2)
lines_sec2xy(Sensors_MK_xy,"GF17",2,iWD17[2],lty=2)
lines_sec2xy(Sensors_MK_xy,"GF17",2,iWD17[1])
lines_sec2xy(Sensors_MK_xy,"GF17",1,iWD17[2])
lines_sec2xy(Mitte_xy,"GF17",1,iWD17M[1],col="pink",lty=3)
lines_sec2xy(Mitte_xy,"GF17",1,iWD17M[2],col="pink",lty=3)
lines_sec2xy(Sensors_MK_xy,"GF17",1,iWD17M[2],col="pink")
lines_sec2xy(Sensors_MK_xy,"GF17",1,iWD17M[1],lty=2,col="pink")
lines_sec2xy(Sensors_MK_xy,"GF17",2,iWD17M[1],col="pink")
lines_sec2xy(Sensors_MK_xy,"GF17",2,iWD17M[2],lty=2,col="pink")

## GF18
# iWD18neu <- c(22,68)
PlotOnStaticMap(STO_Map)
plot(Sensors_MK_xy[Sensors_MK_xy[,1] %in% c("GF18"),],sensor.text.args=list(labels="",cex=2),lines.args=list(lwd=2,col="orange",lty=3),points.args = list(pch = 20, cex = 2, col ="orange"),add=TRUE)
lines_sec2xy(Sensors_MK_xy,"GF18",1,iWD18[1],lty=2)
lines_sec2xy(Sensors_MK_xy,"GF18",2,iWD18[2],lty=2)
lines_sec2xy(Sensors_MK_xy,"GF18",2,iWD18[1])
lines_sec2xy(Sensors_MK_xy,"GF18",1,iWD18[2])
lines_sec2xy(Mitte_xy,"GF18",1,iWD18M[1],col="pink",lty=3)
lines_sec2xy(Mitte_xy,"GF18",1,iWD18M[2],col="pink",lty=3)
lines_sec2xy(Sensors_MK_xy,"GF18",1,iWD18M[2],col="pink")
lines_sec2xy(Sensors_MK_xy,"GF18",1,iWD18M[1],lty=2,col="pink")
lines_sec2xy(Sensors_MK_xy,"GF18",2,iWD18M[1],col="pink")
lines_sec2xy(Sensors_MK_xy,"GF18",2,iWD18M[2],lty=2,col="pink")

## GF16
# iWD16neu <- c(22,68)
PlotOnStaticMap(STO_Map)
plot(Sensors_MK_xy[Sensors_MK_xy[,1] %in% c("GF16"),],sensor.text.args=list(labels="",cex=2),lines.args=list(lwd=2,col="orange",lty=3),points.args = list(pch = 20, cex = 2, col ="orange"),add=TRUE)
lines_sec2xy(Sensors_MK_xy,"GF16",1,iWD16[1],lty=2)
lines_sec2xy(Sensors_MK_xy,"GF16",2,iWD16[2],lty=2)
lines_sec2xy(Sensors_MK_xy,"GF16",2,iWD16[1])
lines_sec2xy(Sensors_MK_xy,"GF16",1,iWD16[2])
lines_sec2xy(Mitte_xy,"GF16",1,iWD16M[1],col="pink",lty=3)
lines_sec2xy(Mitte_xy,"GF16",1,iWD16M[2],col="pink",lty=3)
lines_sec2xy(Sensors_MK_xy,"GF16",1,iWD16M[2],col="pink")
lines_sec2xy(Sensors_MK_xy,"GF16",1,iWD16M[1],lty=2,col="pink")
lines_sec2xy(Sensors_MK_xy,"GF16",2,iWD16M[1],col="pink")
lines_sec2xy(Sensors_MK_xy,"GF16",2,iWD16M[2],lty=2,col="pink")

## GF25
# iWD25neu <- c(22,68)
PlotOnStaticMap(STO_Map)
plot(Sensors_MK_xy[Sensors_MK_xy[,1] %in% c("GF25"),],sensor.text.args=list(labels="",cex=2),lines.args=list(lwd=2,col="orange",lty=3),points.args = list(pch = 20, cex = 2, col ="orange"),add=TRUE)
lines_sec2xy(Sensors_MK_xy,"GF25",1,iWD25[1],lty=2)
lines_sec2xy(Sensors_MK_xy,"GF25",2,iWD25[2],lty=2)
lines_sec2xy(Sensors_MK_xy,"GF25",2,iWD25[1])
lines_sec2xy(Sensors_MK_xy,"GF25",1,iWD25[2])
lines_sec2xy(Mitte_xy,"GF25",1,iWD25M[1],col="pink",lty=3)
lines_sec2xy(Mitte_xy,"GF25",1,iWD25M[2],col="pink",lty=3)
lines_sec2xy(Sensors_MK_xy,"GF25",1,iWD25M[2],col="pink")
lines_sec2xy(Sensors_MK_xy,"GF25",1,iWD25M[1],lty=2,col="pink")
lines_sec2xy(Sensors_MK_xy,"GF25",2,iWD25M[1],col="pink")
lines_sec2xy(Sensors_MK_xy,"GF25",2,iWD25M[2],lty=2,col="pink")


########### QV3
PlotOnStaticMap(STO_Map)
plot(Sensors_QV3_xy[Sensors_QV3_xy[,1] %in% c("GF18","GF26"),],sensor.text.args=list(labels="",cex=2),lines.args=list(lwd=2,col="orange",lty=3),points.args = list(pch = 20, cex = 2, col ="orange"),add=TRUE)
lines_sec2xy(Sensors_QV3_xy,"GF18",1,iWDQV3[2],lty=2)
lines_sec2xy(Sensors_QV3_xy,"GF18",2,iWDQV3[1],lty=2)
lines_sec2xy(Sensors_QV3_xy,"GF18",1,iWDQV3[1])
lines_sec2xy(Sensors_QV3_xy,"GF18",2,iWDQV3[2])
lines_sec2xy(MitteQV3_xy,"GFM",1,iWDMQV3[1],col="pink",lty=3)
lines_sec2xy(MitteQV3_xy,"GFM",1,iWDMQV3[2],col="pink",lty=3)
lines_sec2xy(Sensors_QV3_xy,"GF18",1,iWDMQV3[2],lty=2,col="pink")
lines_sec2xy(Sensors_QV3_xy,"GF18",1,iWDMQV3[1],col="pink")
lines_sec2xy(Sensors_QV3_xy,"GF18",2,iWDMQV3[1],lty=2,col="pink")
lines_sec2xy(Sensors_QV3_xy,"GF18",2,iWDMQV3[2],col="pink")

graphics.off()


## nach Karte Filtern
iWD <- list(
	iWD_GF16 = iWD16
	,iWD_GF17 = iWD17
	,iWD_GF18 = iWD18
	,iWD_GF25 = iWD25
)

## nach Mittelpunkt Pfad filtern (ARA Methode)
iWDM <- list(
	iWD_GF16 = iWD16M
	,iWD_GF17 = iWD17M
	,iWD_GF18 = iWD18M
	,iWD_GF25 = iWD25M
)

##################################
### Emiss (unfiltered) vs. WD  ###
##################################

par(mfrow=c(2,5),mar=c(4,4,0,0))
plot(Q_GF17 ~ I((WD_WS2 + 180) %% 360 - 180),Result[Campaign == "MK"],xlab="WD WS", col= CH4Cols["GF17"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD17 + 180) %% 360 - 180)
for(i in names(SonicCols)){
plot(Q_GF17 ~ I((WD_WS2 + 180) %% 360 - 180),Result[Campaign == "MK" & Sonic == i],xlab="WD WS", col= CH4Cols["GF17"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD17 + 180) %% 360 - 180)
}
plot(Q_GF17 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MK"],xlab="WD WS", col= CH4Cols["GF17"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD17 + 180) %% 360 - 180)
for(i in names(SonicCols)){
plot(Q_GF17 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MK" & Sonic == i],xlab=paste0("WD ", i), col= CH4Cols["GF17"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD17 + 180) %% 360 - 180)
}

par(mfrow=c(2,5),mar=c(4,4,0,0))
plot(Q_GF18 ~ I((WD_WS2 + 180) %% 360 - 180),Result[Campaign == "MK"],xlab="WD WS", col= CH4Cols["GF18"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD18 + 180) %% 360 - 180)
for(i in names(SonicCols)){
plot(Q_GF18 ~ I((WD_WS2 + 180) %% 360 - 180),Result[Campaign == "MK" & Sonic == i],xlab="WD WS", col= CH4Cols["GF18"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD18 + 180) %% 360 - 180)
}
plot(Q_GF18 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MK"],xlab="WD WS", col= CH4Cols["GF18"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD18 + 180) %% 360 - 180)
for(i in names(SonicCols)){
plot(Q_GF18 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MK" & Sonic == i],xlab=paste0("WD ", i), col= CH4Cols["GF18"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD18 + 180) %% 360 - 180)
}

par(mfrow=c(2,5),mar=c(4,4,0,0))
plot(Q_GF16 ~ I((WD_WS2 + 180) %% 360 - 180),Result[Campaign == "MK"],xlab="WD WS", col= CH4Cols["GF16"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD16 + 180) %% 360 - 180)
for(i in names(SonicCols)){
plot(Q_GF16 ~ I((WD_WS2 + 180) %% 360 - 180),Result[Campaign == "MK" & Sonic == i],xlab="WD WS", col= CH4Cols["GF16"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD16 + 180) %% 360 - 180)
}
plot(Q_GF16 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MK"],xlab="WD WS", col= CH4Cols["GF16"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD16 + 180) %% 360 - 180)
for(i in names(SonicCols)){
plot(Q_GF16 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MK" & Sonic == i],xlab=paste0("WD ", i), col= CH4Cols["GF16"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD16 + 180) %% 360 - 180)
}

par(mfrow=c(2,5),mar=c(4,4,0,0))
plot(Q_GF25 ~ I((WD_WS2 + 180) %% 360 - 180),Result[Campaign == "MK"],xlab="WD WS", col= CH4Cols["GF25"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD25 + 180) %% 360 - 180)
for(i in names(SonicCols)){
plot(Q_GF25 ~ I((WD_WS2 + 180) %% 360 - 180),Result[Campaign == "MK" & Sonic == i],xlab="WD WS", col= CH4Cols["GF25"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD25 + 180) %% 360 - 180)
}
plot(Q_GF25 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MK"],xlab="WD WS", col= CH4Cols["GF25"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD25 + 180) %% 360 - 180)
for(i in names(SonicCols)){
plot(Q_GF25 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MK" & Sonic == i],xlab=paste0("WD ", i), col= CH4Cols["GF25"], ylim = c(-20,10),xlim=c(-10,110))
abline(v=(iWD25 + 180) %% 360 - 180)
}


---> WD filtering seems to work better if the wind direction of the individual sonics are used instead of the weather staion.


### Filterung nach Karte und einzelne Sonics als Referenz für die Windrichtung
# MK
for(x in c("GF16","GF17","GF18","GF25")){
	if(x != "GF17"){
	Result[Campaign == "MK" & (WD < iWD[[paste0("iWD_",x)]][1] | WD > iWD[[paste0("iWD_",x)]][2]),
		grep(paste0("Q_",x),names(Result),value=TRUE)	 := NA_real_ ]
		} else {
	Result[Campaign == "MK" & (WD < iWD[[paste0("iWD_",x)]][1] & WD > iWD[[paste0("iWD_",x)]][2]),
		grep(paste0("Q_",x),names(Result),value=TRUE)	 := NA_real_ ]
		}
}
# QV3
for(x in c("GF16","GF17","GF18","GF25","GF26")){
	Result[Campaign == "QV3" & (WD < iWDQV3[1] - 5 & WD > iWDQV3[2] + 5),
		grep(paste0("Q_",x),names(Result),value=TRUE)	 := NA_real_ ]
}


###########################
### bLS model Filterung ###
###########################

y_lim <- c(-10,30)
y_lim2 <- c(-10,200)
y_lim3 <- range(Result[Q_MFC > 0,.(Q_GF16,Q_GF17,Q_GF18,Q_GF25)],na.rm=TRUE)

### Ustar
# MK
par(mfrow=c(2,2))
for(i in names(SonicCols)){
plot(Q_GF17 ~ Ustar,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("u* ",i),ylim=y_lim,xlim=c(0,0.5))
points(Q_GF18 ~ Ustar,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ Ustar,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ Ustar,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF17","GF18","GF16","GF25")])
abline(h=0,v=c(0.05,0.1,0.15),col="blue")
}
# QV3
for(i in names(SonicCols)){
plot(Q_GF26 ~ Ustar,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF26"],ylab="Q [kg/h]",xlab=paste0("u* ",i),ylim=y_lim,xlim=c(0,0.5))
points(Q_GF17 ~ Ustar,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF17"])
points(Q_GF18 ~ Ustar,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ Ustar,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ Ustar,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(0.05,0.1,0.15),col="blue")
}
# nur wenn MFC an
for(i in names(SonicCols)){
plot(Q_GF17 ~ Ustar,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("u* ",i),ylim=y_lim,xlim=c(0,0.5))
points(Q_GF18 ~ Ustar,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ Ustar,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ Ustar,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(0.05,0.1,0.15),col="blue")
}
# ---> keine u* Filterung nötig. Und wenn, dann ohne Probleme bei 0.1 oder sogar 0.15
---> Just take 0.15 and nothing else

### sUu
par(mfrow=c(2,2))
for(i in names(SonicCols)){
plot(Q_GF17 ~ sUu,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("sUu ",i),ylim=y_lim,xlim=c(0,16))
points(Q_GF18 ~ sUu,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sUu,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sUu,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF17","GF18","GF16","GF25")])
abline(h=0,v=c(4,5,6),col="blue")
}
# QV3
for(i in names(SonicCols)){
plot(Q_GF26 ~ sUu,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF26"],ylab="Q [kg/h]",xlab=paste0("sUu",i),ylim=y_lim,xlim=c(0,16))
points(Q_GF17 ~ sUu,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF17"])
points(Q_GF18 ~ sUu,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sUu,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sUu,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(4,5,6),col="blue")
}
# nur wenn MFC an
for(i in names(SonicCols)){
plot(Q_GF17 ~ sUu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("sUu",i),ylim=y_lim,xlim=c(0,16))
points(Q_GF18 ~ sUu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sUu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sUu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(4,5,6),col="blue")
}

# ---> 4 oder 6. 4 schient okay zu sein, da nur bei SonicA und SonicC Werte darüber
---> do not filter

### sVu
par(mfrow=c(2,2))
for(i in names(SonicCols)){
plot(Q_GF17 ~ sVu,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("sVu ",i),ylim=y_lim,xlim=c(0,16))
points(Q_GF18 ~ sVu,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sVu,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sVu,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF17","GF18","GF16","GF25")])
abline(h=0,v=c(4,5,6),col="blue")
}
# QV3
for(i in names(SonicCols)){
plot(Q_GF26 ~ sVu,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF26"],ylab="Q [kg/h]",xlab=paste0("sVu",i),ylim=y_lim,xlim=c(0,16))
points(Q_GF17 ~ sVu,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF17"])
points(Q_GF18 ~ sVu,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sVu,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sVu,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(4,5,6),col="blue")
}
# nur wenn MFC an
for(i in names(SonicCols)){
plot(Q_GF17 ~ sVu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("sVu",i),ylim=y_lim,xlim=c(0,16))
points(Q_GF18 ~ sVu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sVu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sVu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(4,5,6),col="blue")
}

# ---> 4 oder 5. 4 schient okay zu sein, da nur bei SonicA und SonicC Werte darüber
---> do not filter

### sWu
par(mfrow=c(2,2))
for(i in names(SonicCols)){
plot(Q_GF17 ~ sWu,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("sWu ",i),ylim=y_lim,xlim=c(0,10))
points(Q_GF18 ~ sWu,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sWu,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sWu,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF17","GF18","GF16","GF25")])
abline(h=0,v=c(2,3,4),col="blue")
}
# QV3
for(i in names(SonicCols)){
plot(Q_GF26 ~ sWu,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF26"],ylab="Q [kg/h]",xlab=paste0("sWu",i),ylim=y_lim,xlim=c(0,10))
points(Q_GF17 ~ sWu,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF17"])
points(Q_GF18 ~ sWu,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sWu,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sWu,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(2,3,4),col="blue")
}
# nur wenn MFC an
for(i in names(SonicCols)){
plot(Q_GF17 ~ sWu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("sWu",i),ylim=y_lim,xlim=c(0,10))
points(Q_GF18 ~ sWu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ sWu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ sWu,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(2,3,4),col="blue")
}

# ---> 2 oder 4. Besser mal 4 oder nicht filtern, da bei SonicA viele anscheinend sinnvolle Werte über 2.
---> do not filter

### C0
par(mfrow=c(2,2))
for(i in names(SonicCols)){
plot(Q_GF17 ~ C0,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("C0 ",i),ylim=y_lim,xlim=c(0,50))
points(Q_GF18 ~ C0,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ C0,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ C0,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF17","GF18","GF16","GF25")])
abline(h=0,v=c(7,10,30),col="blue")
}
# QV3
for(i in names(SonicCols)){
plot(Q_GF26 ~ C0,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF26"],ylab="Q [kg/h]",xlab=paste0("C0",i),ylim=y_lim,xlim=c(0,50))
points(Q_GF17 ~ C0,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF17"])
points(Q_GF18 ~ C0,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ C0,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ C0,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(7,10,30),col="blue")
}
# nur wenn MFC an
for(i in names(SonicCols)){
plot(Q_GF17 ~ C0,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("C0",i),ylim=y_lim,xlim=c(0,50))
points(Q_GF18 ~ C0,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ C0,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ C0,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(7,10,30),col="blue")
}

# ---> z.T über 10, darum vielleicht 30 obwohl nicht plausibel. SonicA und SonicC betroffen
----> do not filter

### L
# MK
for(i in names(SonicCols)){
plot(Q_GF17 ~ I(1/L),Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF17"],xlab=paste0("1/L ",i),ylab="Q [kg/h]",ylim=y_lim,xlim=c(-2,2))
points(Q_GF18 ~ I(1/L),Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ I(1/L),Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ I(1/L),Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF17","GF18","GF16","GF25")])
abline(h=0,v=c(-0.5,0.5),col="blue")
}
# QV3
for(i in names(SonicCols)){
plot(Q_GF26 ~ I(1/L),Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF26"],xlab=paste0("1/L ",i),ylab="Q [kg/h]",ylim=y_lim,xlim=c(-2,2))
points(Q_GF17 ~ I(1/L),Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF17"])
points(Q_GF18 ~ I(1/L),Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ I(1/L),Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ I(1/L),Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(-0.5,0.5),col="blue")
}
# nur wenn MFC an
for(i in names(SonicCols)){
plot(Q_GF17 ~ I(1/L),Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF17"],xlab=paste0("1/L ",i),ylab="Q [kg/h]",ylim=y_lim,xlim=c(-2,2))
points(Q_GF18 ~ I(1/L),Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ I(1/L),Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ I(1/L),Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(-0.5,0.5),col="blue")
}

--> abs(0.5) or could even take abs(0.2)


### Zo
# MK
par(mfrow=c(2,2))
for(i in names(SonicCols)){
plot(Q_GF17 ~ Zo,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("Zo ",i),ylim=y_lim,xlim=c(0,1))
points(Q_GF18 ~ Zo,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ Zo,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ Zo,Result[Campaign == "MK" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF17","GF18","GF16","GF25")])
abline(h=0,v=c(0.05,0.07,0.1),col="blue")
}
# QV3
for(i in names(SonicCols)){
plot(Q_GF26 ~ Zo,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF26"],ylab="Q [kg/h]",xlab=paste0("Zo",i),ylim=y_lim,xlim=c(0,1))
points(Q_GF17 ~ Zo,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF17"])
points(Q_GF18 ~ Zo,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ Zo,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ Zo,Result[Campaign == "QV3" & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(0.05,0.07,0.1),col="blue")
}
# nur wenn MFC an
for(i in names(SonicCols)){
plot(Q_GF17 ~ Zo,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF17"],ylab="Q [kg/h]",xlab=paste0("Zo",i),ylim=y_lim,xlim=c(0,0.2))
points(Q_GF18 ~ Zo,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF18"])
points(Q_GF16 ~ Zo,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF16"])
points(Q_GF25 ~ Zo,Result[Q_MFC > 0 & Sonic == i],pch=19,col=CH4Cols["GF25"])
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])
abline(h=0,v=c(0.05,0.07,0.1),col="blue")
}

# ---> 0.07 oder 0.1 bzw gar nicht nötig zu filtern
---> do not filter

# für Plot object kopieren, damit ich ungefilterted Daten plotten kann

Res_unfil <- copy(Result)

# erstmals sehr grob und nur einzelne Parameter filtern
# Result[sUu >= 4 | sVu >= 4 | sWu >= 2 | abs(1/L) >= 0.5 | Zo >= 0.07 | C0 >= 10, grep("Q_GF",names(Result),value = TRUE) := NA_real_]
# Result[ abs(1/L) >= 0.5 |  C0 >= 10, grep("Q_GF",names(Result),value = TRUE) := NA_real_] # all the other parameters are actually not necessary and a bit arbitrary.
# above was some old filter. I rather exclude C0 and take the others.
Result[Ustar < 0.15, grep("Q_GF",names(Result),value = TRUE) := NA_real_]


## alle Sonics und GFs
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


#############################
### Emissionen betrachten ###
#############################

### alle Emissionen mit MFC
par(mfrow=c(3,1))
# Emiss
plot(Q_GF16 ~ st, Result[Campaign %in% c("MK","QV3")],type="l",col=CH4Cols["GF16"],ylim=c(-10,10),ylab="Q [kg/h]")
lines(Q_GF17 ~ st, Result[Campaign %in% c("MK","QV3")],col=CH4Cols["GF17"])
lines(Q_GF18 ~ st, Result[Campaign %in% c("MK","QV3")],col=CH4Cols["GF18"])
lines(Q_GF25 ~ st, Result[Campaign %in% c("MK","QV3")],col=CH4Cols["GF25"])
lines(Q_GF26 ~ st, Result[Campaign %in% c("MK","QV3")],col=CH4Cols["GF26"])
lines(Q_MFC ~ st, Result[Campaign %in% c("MK","QV3")],col="black",lwd=2)
# wind speed
plot(U_sonic ~ st, Result[Campaign %in% c("MK","QV3") & Sonic == "SonicB"],type="l",col="blue",ylab="U_sonic [m/s]")
# wind direction
plot(WD_WS2 ~ st, Result[Campaign %in% c("MK","QV3") & Sonic == "SonicB"],type="l",col="red",ylab="wind direction")

### alle Emissionen mit MFC
par(mfrow=c(3,1))
# Emiss
plot(Q_GF16 ~ st, Result[Campaign %in% c("MK")],type="l",col=CH4Cols["GF16"],ylim=c(-10,10),ylab="Q [kg/h]")
lines(Q_GF17 ~ st, Result[Campaign %in% c("MK")],col=CH4Cols["GF17"])
lines(Q_GF18 ~ st, Result[Campaign %in% c("MK")],col=CH4Cols["GF18"])
lines(Q_GF25 ~ st, Result[Campaign %in% c("MK")],col=CH4Cols["GF25"])
lines(Q_GF26 ~ st, Result[Campaign %in% c("MK")],col=CH4Cols["GF26"])
lines(Q_MFC ~ st, Result[Campaign %in% c("MK")],col="black",lwd=2)
legend("topright",legend=c("GF26","GF17","GF18","GF16","GF25"), fill=CH4Cols[c("GF26","GF17","GF18","GF16","GF25")])

# wind speed
plot(U_sonic ~ st, Result[Campaign %in% c("MK") & Sonic == "SonicB"],type="l",col="blue",ylab="U_sonic [m/s]")
# wind direction
plot(WD_WS2 ~ st, Result[Campaign %in% c("MK") & Sonic == "SonicB"],type="l",col="red",ylab="wind direction")

graphics.off()
###########################
### nochmals WD filtern ###
###########################

par(mfrow=c(2,5),mar=c(4,4,0,0))
plot(Q_GF17 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MK"],xlab="WD all Sonics", col= CH4Cols["GF17"], ylim = c(-20,10),xlim=c(-10,110))
abline(h=0, v=(iWD17 + 180) %% 360 - 180)
abline(v=c(21),lty=3,col="red")
for(i in names(SonicCols)){
plot(Q_GF17 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MK" & Sonic == i],xlab=paste0("WD ", i), col= CH4Cols["GF17"], ylim = c(-20,10),xlim=c(-10,110))
abline(h=0, v=(iWD17 + 180) %% 360 - 180)
abline(v=c(21),lty=3,col="red")
}
plot(Q_GF18 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MK"],xlab="WD all Sonics", col= CH4Cols["GF18"], ylim = c(-20,10),xlim=c(-10,110))
abline(h=0, v=(iWD18 + 180) %% 360 - 180)
abline(v=c(21,71),lty=3,col="red")
for(i in names(SonicCols)){
plot(Q_GF18 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MK" & Sonic == i],xlab=paste0("WD ", i), col= CH4Cols["GF18"], ylim = c(-20,10),xlim=c(-10,110))
abline(h=0, v=(iWD18 + 180) %% 360 - 180)
abline(v=c(21,71),lty=3,col="red")
}

par(mfrow=c(2,5),mar=c(4,4,0,0))
plot(Q_GF16 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MK"],xlab="WD all Sonics", col= CH4Cols["GF16"], ylim = c(-20,10),xlim=c(-10,110))
abline(h=0, v=(iWD16 + 180) %% 360 - 180)
abline(v=c(30,64),lty=3,col="red")
for(i in names(SonicCols)){
plot(Q_GF16 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MK" & Sonic == i],xlab=paste0("WD ", i), col= CH4Cols["GF16"], ylim = c(-20,10),xlim=c(-10,110))
abline(h=0, v=(iWD16 + 180) %% 360 - 180)
abline(v=c(30,64),lty=3,col="red")
}
plot(Q_GF25 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MK"],xlab="WD all Sonics", col= CH4Cols["GF25"], ylim = c(-20,10),xlim=c(-10,110))
abline(h=0, v=(iWD25 + 180) %% 360 - 180)
abline(v=c(30,61),lty=3,col="red")
for(i in names(SonicCols)){
plot(Q_GF25 ~ I((WD + 180) %% 360 - 180),Result[Campaign == "MK" & Sonic == i],xlab=paste0("WD ", i), col= CH4Cols["GF25"], ylim = c(-20,10),xlim=c(-10,110))
abline(h=0, v=(iWD25 + 180) %% 360 - 180)
abline(v=c(30,61),lty=3,col="red")
}

## nach Karte Filtern
iWD_refined <- list(
	iWD_GF16 = c(30,64)
	,iWD_GF17 = c(21,iWD17[2])
	,iWD_GF18 = c(21,71)
	,iWD_GF25 = c(30,61)
)



### Filterung nach Karte und einzelne Sonics als Referenz für die Windrichtung
# MK
for(x in c("GF16","GF17","GF18","GF25")){
	Result[Campaign == "MK" & (WD < iWD_refined[[paste0("iWD_",x)]][1] | WD > iWD_refined[[paste0("iWD_",x)]][2]),
		grep(paste0("Q_",x),names(Result),value=TRUE)	 := NA_real_ ]
}

# # QV3 # this is not done and probably not necessary
# for(x in c("GF16","GF17","GF18","GF25","GF26")){
# 	Result[Campaign == "QV3" & (WD < iWDQV3[1] - 5 & WD > iWDQV3[2] + 5),
# 		grep(paste0("Q_",x),names(Result),value=TRUE)	 := NA_real_ ]
# }


---> vielleicht noch eine CQ Filterung machen

##########################
### CQ Werte anschauen ###
##########################

CQ_thresh <- 8E-5

par(mfrow=c(2,5),mar=c(4,4,0,0))
plot(Q_GF17 ~ CQ_GF17,Result[Campaign == "MK"],xlab=paste0("CQ_GF17 ",i), col= CH4Cols["GF17"], ylim = c(-20,10),xlim=c(0,0.0015))
abline(h=0, v=CQ_thresh)
for(i in names(SonicCols)){
plot(Q_GF17 ~ CQ_GF17,Result[Campaign == "MK" & Sonic == i],xlab=paste0("CQ_GF17 ", i), col= CH4Cols["GF17"], ylim = c(-20,10),xlim=c(0,0.0015))
abline(h=0, v=CQ_thresh)
}
plot(Q_GF18 ~ CQ_GF18,Result[Campaign == "MK"],xlab=paste0("CQ_GF18 ",i), col= CH4Cols["GF18"], ylim = c(-20,10),xlim=c(0,0.0015))
abline(h=0, v=CQ_thresh)
for(i in names(SonicCols)){
plot(Q_GF18 ~ CQ_GF18,Result[Campaign == "MK" & Sonic == i],xlab=paste0("CQ_GF18 ", i), col= CH4Cols["GF18"], ylim = c(-20,10),xlim=c(0,0.0015))
abline(h=0, v=CQ_thresh)
}

par(mfrow=c(2,5),mar=c(4,4,0,0))
plot(Q_GF16 ~ CQ_GF16,Result[Campaign == "MK"],xlab=paste0("CQ_GF16 ",i), col= CH4Cols["GF16"], ylim = c(-20,10),xlim=c(0,0.0015))
abline(h=0, v=CQ_thresh)
for(i in names(SonicCols)){
plot(Q_GF16 ~ CQ_GF16,Result[Campaign == "MK" & Sonic == i],xlab=paste0("CQ_GF16 ", i), col= CH4Cols["GF16"], ylim = c(-20,10),xlim=c(0,0.0015))
abline(h=0, v=CQ_thresh)
}
plot(Q_GF25 ~ CQ_GF25,Result[Campaign == "MK"],xlab=paste0("CQ_GF25 ",i), col= CH4Cols["GF25"], ylim = c(-20,10),xlim=c(0,0.0015))
abline(h=0, v=CQ_thresh)
for(i in names(SonicCols)){
plot(Q_GF25 ~ CQ_GF25,Result[Campaign == "MK" & Sonic == i],xlab=paste0("CQ_GF25 ", i), col= CH4Cols["GF25"], ylim = c(-20,10),xlim=c(0,0.0015))
abline(h=0, v=CQ_thresh)
}

---> CQ Filterung nicht sinnvoll. Da für jeden GF und Sonic anders gewählt werden müsste.


graphics.off()
###########################

saveRDS(Result,file=file.path(PfadRSaves,"STO_Results_orig_P23.rds"))
