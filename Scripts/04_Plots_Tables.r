
##############################################
##############################################
#####                                    #####
#####    Creation of plots and tables    #####
#####                                    #####
##############################################
##############################################

# Author: Marcel Bühler
# Date: August 5, 2024
# Contact: mb@bce.au.dk or Christoph Häni christoph.haeni@bfh.ch
# Description: With this script one can recreate all the plots and values in the tables of the publication, the supplement, and the initial submission.
#
# Note: This code was written by Marcel Bühler and is intended to follow the publication 'Applicability of the inverse dispersion method to measure emissions from animal housings' in AMT. 
# Please feel free to use and modify it, but attribution is appreciated.


#################
### Libraries ###
#################

## there might be some libraries that are not necessary.
library(ibts)
library(RgoogleMaps)
library(bLSmodelR)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggpointdensity) # this might not be necessary
library(shape) # this might not be necessary
library(plotly) # this might not be necessary


#############
### Paths ###
#############

PathData <- "Path to /data"		
PathRSaves <- "Path to /RSaves"
PathFigures <- 'Path to /Figures'


#################
### Functions ###
#################

source("https://raw.githubusercontent.com/hafl-gel/gel-scripts/main/windrose.r")
source("https://raw.githubusercontent.com/hafl-gel/gel-scripts/main/wgs84-ch1903.r")
source(file.path(file.path(dirname(PathRSaves),"Other/contourXY.r")))
source(file.path(file.path(dirname(PathRSaves),"Other/contourXZ.r")))


#################
### Campaigns ###
#################

	IC1 <- "05.03.2021 to 10.03.2021 18:00"
	MC <- "18.03.2021 11:00 - 21.03.2021 14:00"
	IC2	<- "21.03.2021 14:00 to "


#################
### load data ###
#################

	## Geometry
	load(file.path(PathRSaves,"Geometry.RData"))
	STO_Map <- ReadMapTile(file.path(PathFigures,"STO_GoogleMaps.png"))
	## weather station
	WS700_2 <- readRDS(file.path(PathRSaves, "WS700.rds"))
	## Emissions
	Result <- readRDS(file=file.path(PathRSaves,"Results_orig_P23.rds"))
	## Sonics
	load(file.path(PathRSaves,"Sonics_10min.RData"))
	## MFC
	MFC <- readRDS(file.path(PathRSaves,"MFC.rds"))
	## KELLER pressure sensor
	CC_Sensor <- readRDS(file.path(PathRSaves,"CC_Sensor_raw.rds"))
	# Concentration data
	GFs10min <- readRDS(file.path(PathRSaves,"Conc_P23_10min.rds"))

	## contour plots (needed for a plot in the initial submission)
	# XY_SonicA <- readRDS(file.path(PathRSaves,"XY_SonicA_MC.rds"))
	# XY_SonicB <- readRDS(file.path(PathRSaves,"XY_SonicB_MC.rds"))
	# XY_SonicC <- readRDS(file.path(PathRSaves,"XY_SonicC_MC.rds"))
	# XY_Sonic2 <- readRDS(file.path(PathRSaves,"XY_Sonic2_MC.rds"))
	## Emissions of WD variation.
	# Result_var <- readRDS(file=file.path(PathRSaves,"Result_variation.rds")) # this is not provided as not used in the final publication. But you can calculate it with the provided data.
	

##########################
### create XY Geometry ###
##########################

Sensors_MC_xy <- ch_to_map(STO_Map,Sensors_MC)
Sources_xy <- ch_to_map(STO_Map,Source)
Tree_xy <- ch_to_map(STO_Map,Tree)

centre <- colMeans(Source[,c(2,3)])
Sonics_MC <- Sensors_MC[Sensors_MC[,1] %in% c('SonicA','SonicB','SonicC','Sonic2'),]


#######################
### Weather station ###
#######################

WS2 <- shift_time(WS700_2[[1]]$data,"-5mins") # Apply 5min offset correction
WS2[which(WS2$WS_0460 == 0),c("WS_0460","WD_corr")] <- NA_real_ 	# remove no wind data
WS_10min <- pool(WS2,granularity="10mins",tz="Etc/GMT-1",st="2021-03-18 15:10")
WS_1min <- pool(WS2[,c("WS_0480","WD_corr")],granularity="1mins",tz="Etc/GMT-1",st="2021-03-18 15:10")
WS2_10min <- merge.ibts(pool(WS700_2[[2]]$data[,c("WS_0260", "WS_0360", "WS_0620", "WS_0625", "WS_0820")],granularity="10mins",st.to="07.12.2020 12:00"),WS700_2[[1]]$data[,c("WS_0160","WS_0480","WD_corr")])
names(WS2_10min) <- c("relHum","Press","Precip1","Precip2","Precip3","Temp","WS","WD_WS")

WS2_10min['19.03.2021 11:10 - 19.03.2021 11:30', c('relHum','Press','Temp','WS','WD_WS')] <- NA_real_ # set the values during power outage to NA
WS2_dt <- cbind(as.data.table(WS2_10min),st=st(WS2_10min))


##############
### Result ###
##############

## rename Sonics
Result$Sonic_ord2 <- factor(Result$Sonic, levels=c("SonicC","SonicA","Sonic2","SonicB"),labels=c("UA-UW","UA-2.0h","UA-5.3h","UA-8.6h"))
Result$Sonic_ord <- factor(Result$Sonic, levels=c("SonicC","SonicA","Sonic2","SonicB"),labels=c("UA-UW","UA-50m","UA-100m","UA-150m"))

## calcualte recovery, but only use data with a certain flow through the MFC
Result[Q_MFCex > 4,R_GF16 := Q_GF16 / Q_MFCex]
Result[Q_MFCex > 4,R_GF17 := Q_GF17 / Q_MFCex]
Result[Q_MFCex > 4,R_GF18 := Q_GF18 / Q_MFCex]
Result[Q_MFCex > 4,R_GF25 := Q_GF25 / Q_MFCex]
Result[Q_MFCex > 4,R_GF26 := Q_GF26 / Q_MFCex]


##############
### Sonics ###
##############

SonicB <- SonicB[,names(SonicB)!='MC']
Sonic2[,"MC"] <-  NA_real_
Sonic2[MC,"MC"] <-  "MC"
Sonic2[IC2,"MC"] <-  "IC2"

SonicAll_10min <- merge(merge(
        merge(
            merge(cbind(as.data.table(SonicB),start_interval=st(SonicB),end_interval = et(SonicB)), cbind(as.data.table(SonicA),start_interval=st(SonicA),end_interval = et(SonicA)), by = c('start_interval', 'end_interval'), suffixes = c('', '.a'))
            , cbind(as.data.table(SonicC),start_interval=st(SonicC),end_interval = et(SonicC)), by = c('start_interval', 'end_interval'), suffixes = c('', '.c'))
                , cbind(as.data.table(Sonic2),start_interval=st(Sonic2),end_interval = et(Sonic2)), by = c('start_interval', 'end_interval'), suffixes = c('.b', '.2'))
                    , cbind(as.data.table(WS_10min),start_interval=st(WS_10min),end_interval=et(WS_10min)), by = c('start_interval','end_interval'))


SonicAll_10min[, ':=' (
	dWD.a = (((WD.a - WD.c) %% 360) + 180) %% 360 -180,
	dWD.2 = (((WD.2 - WD.c) %% 360) + 180) %% 360 -180,
	dWD.b = (((WD.b - WD.c) %% 360) + 180) %% 360 -180,
	dUstar.a = (Ustar.a - Ustar.c),
	dUstar.2 = (Ustar.2 - Ustar.c),
	dUstar.b = (Ustar.b - Ustar.c),
	dsVu.a = (sVu.a - sVu.c),
	dsVu.2 = (sVu.2 - sVu.c),
	dsVu.b = (sVu.b - sVu.c),
	dsWu.a = (sWu.a - sWu.c),
	dsWu.2 = (sWu.2 - sWu.c),
	dsWu.b = (sWu.b - sWu.c),
	dU_sonic.a = (U_sonic.a - U_sonic.c),
	dU_sonic.2 = (U_sonic.2 - U_sonic.c),
	dU_sonic.b = (U_sonic.b - U_sonic.c),
	# diff to weather station
	dWDWS.a = (((WD.a - WD_corr) %% 360) + 180) %% 360 -180,
	dWDWS.2 = (((WD.2 - WD_corr) %% 360) + 180) %% 360 -180,
	dWDWS.b = (((WD.b - WD_corr) %% 360) + 180) %% 360 -180,
	dWDWS.c = (((WD.c - WD_corr) %% 360) + 180) %% 360 -180,
	dU_sonicWS.a = (U_sonic.a - WS_0480),
	dU_sonicWS.2 = (U_sonic.2 - WS_0480),
	dU_sonicWS.b = (U_sonic.b - WS_0480),
	dU_sonicWS.c = (U_sonic.c - WS_0480)
	)]


SonicAll_10min_cast <- melt(SonicAll_10min,id.vars = c('start_interval','MC','WD.c','U_sonic.c','Ustar.c','sVu.c','sWu.c','WD_corr','WS_0480'),
	measure.vars = c('dWD.a','dWD.2','dWD.b','dU_sonic.a','dU_sonic.2','dU_sonic.b',
		'dUstar.a','dUstar.2','dUstar.b','dsVu.a','dsVu.2','dsVu.b','dsWu.a','dsWu.2','dsWu.b',
		'dWDWS.a','dWDWS.2','dWDWS.b','dWDWS.c','dU_sonicWS.a','dU_sonicWS.2','dU_sonicWS.b','dU_sonicWS.c'))

SonicAll_10min_cast[grep('.a$', variable), Sonic := 'UA-2.0h']
SonicAll_10min_cast[grep('.2$', variable), Sonic := 'UA-5.3h']
SonicAll_10min_cast[grep('.b$', variable), Sonic := 'UA-8.6h']
SonicAll_10min_cast[grep('.c$', variable), Sonic := 'UA-UW']

SonicAll_10min_cast$Sonic_ordh <- factor(SonicAll_10min_cast$Sonic, levels=c("UA-2.0h","UA-5.3h","UA-8.6h","UA-UW"),
	labels=c("UA-2.0h","UA-5.3h","UA-8.6h","UA-UW"))

SonicAll_10min_cast$Sonic_ord <- factor(SonicAll_10min_cast$Sonic, levels=c("UA-UW","UA-2.0h","UA-5.3h","UA-8.6h"),
	labels=c("UA-UW","UA-50m","UA-100m","UA-150m"))

SonicAll_10min_cast[grep('WD', variable), top_variable := 'WD']
SonicAll_10min_cast[grep('U_sonic', variable), top_variable := 'U_sonic']
SonicAll_10min_cast[grep('Ustar', variable), top_variable := 'Ustar']
SonicAll_10min_cast[grep('sVu', variable), top_variable := 'sVu']
SonicAll_10min_cast[grep('sWu', variable), top_variable := 'sWu']

SonicAll_10min_cast[grep('WS',variable),comparison := 'WS']
SonicAll_10min_cast[!grep('WS',variable),comparison := 'UW']


#############################
### MFC & Pressure sensor ###
#############################

isub <- "19.03.2021 09:00 - 20.03.2021 09:00"
MFC_p <- MFC[isub]
CC_p <- pool(CC_Sensor[isub],granularity="1mins",st.to="19.03.2021 09:00:00")
tvect <- data.table(st=seq(0,1440)*60 + parse_date_time3("19.03.2021 09:00",tz="Etc/GMT-1"))

MFC_dt <- data.table(cbind(MFC_p,st=st(MFC_p)))
CC_dt <- merge(tvect,as.data.table(cbind(CC_p,st=st(CC_p))),by="st",all.x=TRUE)
CC_dt[,Press_CCex := 0]
CC_dt[!is.na(Press_CC),Press_CCex := Press_CC]
coeff_MFC <- MFC_dt[,max(Q_MFC,na.rm=TRUE) / 16]
coeff_CC <- CC_dt[,max(Press_CC,na.rm=TRUE) / 9]


#-----------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------#
					#/////////////////////////////////---    PUBLICATION    ---\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#-----------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------#


#####################################
#####################################
#####                           #####
#####    Figures publication    #####
#####                           #####
#####################################
#####################################

#############
### Figure 1:

# this is a picture


#############
### Figure 2:

# this is a schematic drawn with a design program (Inkscape)


######
### Figure 3: Schematic overview of the measurement setup during the measurement campaign. OP: open-path device. UA: 3D ultrasonic anemometer. UW: upwind. The numbers behind the OP and UA represent the fetch.

{
png(file.path(PathFigures,"Figure_3.png"),width=23,height=15,unit="cm",res=600)
# par(mfrow=c(1,1),mar=c(0.0,0.0,0.0,0.0))
par(mfrow=c(1,1),mar=c(0.1,0.1,0.1,0.1))
plot(type="n",1, xlim=c(centre[1]-255,centre[1]+200),ylim=c(centre[2]-215,centre[2]+97),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Source,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[9,][,4],y0=Sensors_MC[9,][,5],
	x1=Sensors_MC[10,][,4],y1=Sensors_MC[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,code=1)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[9,][,4],y0=Sensors_MC[9,][,5],
	x1=Sensors_MC[10,][,4],y1=Sensors_MC[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)

plot(Tree,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
plot(Sonics_MC,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
plot(Sensors_MC[Sensors_MC[,1] == "WS2",],sensors.text.args=list(labels=""),points.args = list(pch = 17, cex = 2, col ="black"),add=TRUE)
bx <- centre[1]+80
by <- centre[2]-205
dx1 <- 50
dx2 <- 100
dy <- 4
## scale bar
polygon(c(bx,bx+dx1,bx+dx1,bx),c(by,by,by-dy,by-dy),col="black")
polygon(c(bx+dx1,bx+dx2,bx+dx2,bx+dx1),c(by,by,by-dy,by-dy),col="white",border="black")
lines(c(bx,bx),c(by,by+dy*(2/3)),col="black")
lines(c(bx+dx1,bx+dx1),c(by,by+dy*(1/3)),col="black")
lines(c(bx+dx2,bx+dx2),c(by,by+dy*(2/3)),col="black")
text(bx,by+dy*3,"0 m",cex=1.3,adj=c(0.5,0.5),col="black")
text(bx+dx1,by+dy*3,"50 m",cex=1.3,adj=c(0.5,0.5),col="black")
text(bx+dx2,by+dy*3,"100 m",cex=1.3,adj=c(0.5,0.5),col="black")
## North direction
wr.x <- bx+85
wr.y <- by+40  
polygon(c(wr.x,wr.x-5,wr.x),c(wr.y,wr.y-5,wr.y+40),border=1)
polygon(c(wr.x,wr.x+5,wr.x),c(wr.y,wr.y-5,wr.y+40),border=1,col=1)
text(x=wr.x,y=wr.y+50,labels="N",cex=1.5)
# text(x=centre[1]-185,y=centre[2]+120,"MC",cex=1.5)
## label GasFinder
text(x=Sensors_MC[1,][,4]+21,y=Sensors_MC[1,][,5]-4,"OP-8.6h")
text(x=Sensors_MC[3,][,4]+21,y=Sensors_MC[3,][,5]-4,"OP-2.0h")
text(x=Sensors_MC[5,][,4]+21,y=Sensors_MC[5,][,5]-4,"OP-5.3h")
text(x=Sensors_MC[7,][,4]+21,y=Sensors_MC[7,][,5]-4,"OP-12h")
text(x=Sensors_MC[10,][,4]-21,y=Sensors_MC[10,][,5],"OP-UW")
## label Sonics
text(x=Sonics_MC[1,][,4]-16,y=Sonics_MC[1,][,5]-15,"UA-5.3h")
text(x=Sonics_MC[2,][,4]-15,y=Sonics_MC[2,][,5]-15,"UA-2.0h")
text(x=Sonics_MC[3,][,4]-15,y=Sonics_MC[3,][,5]-15,"UA-8.6h")
text(x=Sonics_MC[4,][,4],y=Sonics_MC[4,][,5]-15,"UA-UW")

## Legend
xl <- centre[1] - 270
yl <- centre[2] + 97
dyl <- 16
dxl <- 8
Arrows(x0=xl-dxl,x1=xl+dxl,y0=yl,y1=yl,arr.type="triangle",arr.length=0.1,arr.adj=1)
Arrows(x0=xl+dxl,x1=xl-dxl,y0=yl,y1=yl,arr.type="T",arr.length=0.1,arr.adj=1)
points(x=c(xl,xl),y=c(yl-dyl,yl-dyl*2),pch=c(8,17),cex=1)
points(x=xl,y=yl-dyl*3,pch=21,cex=2,col='grey60',bg='grey95')
rect(xleft=xl-dxl,xright=xl+dxl,ybottom=yl-dyl*4-4,ytop=yl-dyl*4+4,col="grey95")
text(x=rep(xl+dxl*2,5),y=c(yl,yl-dyl,yl-dyl*2,yl-dyl*3,yl-dyl*4),labels=c("OP","UA","Weather station","Tree","Barn"),pos=4)
# lines(x=c(xl-dxl*2,xl-dxl*2,xl+200),y=c(yl-100,yl+dyl,yl+dyl))
lines(x=c(xl-20,xl+100,xl+100),y=c(yl-80,yl-80,yl+20))

dev.off()
}


#############
### Figure 4: Weather conditions as 10 min averages measured with the on-site weather station (temperature) and the UA-UW (wind direction and wind speed) during the measurement campaign. The grey shaded areas indicate the times during which CH4 was released.

limits_WSUA  <- parse_date_time3(c("19.03.2021","21.03.2021"),tz="Etc/GMT-1")

## Temperature from weather station
WS_TempUA <- WS2_dt[st >= parse_date_time3("18.03.2021 20:00",tz='Etc/GMT-1') & st <= parse_date_time3('21.03.2021 04:00',tz='Etc/GMT-1'),{
 	ggplot(.SD,aes(x=st,y=Temp)) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 11:19",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 11:40",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 01:08",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("20.03.2021 01:21",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_hline(yintercept = 0, colour="grey") +
 	geom_line() +
 	scale_x_datetime(limits=limits_WSUA,expand=c(0.02,0.02)) +
 	ylim(-5,10) +
 	ylab("Temperature [°C]") +
 	xlab(NULL) +
 	theme_bw(base_size=18) +
 	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3),
	 axis.text.x = element_blank(), axis.ticks.x = element_blank(),plot.margin = unit(c(0.2, 1.05, -2, 0.05), "cm"))
 }]

## Wind speed and wind direction from UA-UW
Fig_WD_UAUW <- Result[Sonic == 'SonicC' & st >= parse_date_time3("19.03.2021",tz='Etc/GMT-1') & st <= parse_date_time3('21.03.2021',tz='Etc/GMT-1'),{
	ggplot(.SD,aes(x=st)) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 11:19",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 11:40",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 01:08",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("20.03.2021 01:21",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_line(aes(y=(WD - 180) %% 360 - 180)) +
 	geom_line(aes(y=(WD + 180) %% 360 + 180)) +
 	scale_x_datetime(limits=limits_WSUA,expand=c(0.02,0.02)) +
	scale_y_continuous(breaks = seq(0, 360, by=90), limits=c(-90,450),expand=c(-0.15,-0.15)) +
 	ylab("Wind direction [°]") +
 	xlab(NULL) +
 	theme_bw(base_size=18) +
 	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3), 
 		axis.text.x = element_blank(), axis.ticks.x = element_blank(),plot.margin = unit(c(-3, 1.05, -3, 0.05), "cm"))
}]

Fig_WS_UAUW <- Result[Sonic == 'SonicC' & st >= parse_date_time3("19.03.2021",tz='Etc/GMT-1') & st <= parse_date_time3('21.03.2021',tz='Etc/GMT-1'),{
	ggplot(.SD,aes(x=st,y=U_sonic)) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 11:19",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 11:40",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 01:08",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("20.03.2021 01:21",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_line(size=0.5) +
 	scale_x_datetime(limits=limits_WSUA,expand=c(0.02,0.02)) +
	scale_y_continuous(breaks = seq(0, 6, by=2), limits=c(0,6),expand=c(0.02,0.02)) +
 	ylab(expression("Wind speed [m s"^-1*"]")) +
 	xlab(NULL) +
 	theme_bw(base_size=18) +
 	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3),plot.margin = unit(c(-3, 1.05, 0.05, 0.05), "cm"))
}]

Meteo_MK_WSUA <- ggarrange(WS_TempUA,Fig_WD_UAUW,Fig_WS_UAUW,ncol=1,align='hv')

ggsave(file.path(PathFigures, "Figure_4.png"), Meteo_MK_WSUA, width = 30, height = 30/1.8, units="cm")


#############
### Figure 5: (a) Recovery rate for the measurement campaign. The colours indicate the OP used to calculate the recovery rate. (b) Atmospheric stability recorded with the UA-UW. Grey shaded areas are the times during which CH4 was released. The time zone is UTCC1.

limits_recovery_MC <- c(parse_date_time3("2021-03-19 08:00",tz='Etc/GMT-1'), parse_date_time3("2021-03-20 09:00",tz='Etc/GMT-1'))

Fig_recovery_UWSonic_IDM <- Result[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & Sonic_ord2=='UA-UW',{
	ggplot(.SD, aes(x=st)) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 11:19",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 11:40",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 01:08",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("20.03.2021 01:21",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_hline(yintercept=1,linetype=3) +
	geom_line(aes(y=R_GF17,colour="OP-2.0h"),size=0.6) +
	geom_point(aes(y=R_GF17,colour="OP-2.0h"),size=1) +
	geom_line(aes(y=R_GF18,colour="OP-5.3h"),size=0.6) +
	geom_point(aes(y=R_GF18,colour="OP-5.3h"),size=1) +
	geom_line(aes(y=R_GF16,colour="OP-8.6h"),size=0.6) +
	geom_point(aes(y=R_GF16,colour="OP-8.6h"),size=1) +
	geom_line(aes(y=R_GF25,colour="OP-12h"),size=0.6) +
	geom_point(aes(y=R_GF25,colour="OP-12h"),size=1) +
	xlab(NULL) +
	scale_x_datetime(limits=limits_recovery_MC) +
	ylab('IDM recovery rate [%]') +
	scale_colour_manual(name=NULL,values=c('OP-2.0h'='#F8766D','OP-5.3h'='#7CAE00','OP-8.6h'='#00BFC4','OP-12h'='#C77CFF'),breaks=c('OP-2.0h','OP-5.3h','OP-8.6h','OP-12h')) +
	theme_bw(base_size = 18) +
	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3),
	strip.background =element_rect(fill="white"), legend.title=element_blank(), legend.position='top',
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),plot.margin = unit(c(0, 0.1, 0, 0.1), "cm"))
}]


Fig_stability_UWSonic <- Result[st > "2021-03-19 08:00" & st < "2021-03-20 09:00" & Sonic_ord2=='UA-UW',{
	ggplot(.SD, aes(x=st)) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 11:19",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 11:40",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 01:08",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("20.03.2021 01:21",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_hline(yintercept=0,linetype=3) +
	geom_line(aes(y=1/L),size=0.6) +
	geom_point(aes(y=1/L),size=1) +
	scale_y_continuous(breaks=seq(-0.2,0.4,by=0.05),limits=c(-0.4,0.4),expand=c(-0.108,-0.208)) +
	xlab(NULL) +
	scale_x_datetime(limits=limits_recovery_MC) +
	ylab(expression('L'^-1*' [m'^-1*']')) +
	scale_colour_manual(name=NULL,values=c('UA-UW'='black')) +
	theme_bw(base_size = 18) +
	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3),
	legend.title=element_blank(),plot.margin = unit(c(0, 0.1, 0.1, 0.1), "cm"))
}]


Fig_recovery_stability_UWSonic_IDM <- ggarrange(Fig_recovery_UWSonic_IDM,Fig_stability_UWSonic,align='v',ncol=1,heights=c(5,2))


Fig_5 <- Fig_recovery_stability_UWSonic_IDM +
annotate("text", x = 0.983, y = 0.80, label = 'atop(bold("(a)"))', color = "black", parse = TRUE, cex = 5) +
annotate("text", x = 0.983, y = 0.23, label = 'atop(bold("(b)"))', color = "black", parse = TRUE, cex = 5) +
annotate("text", x = 0.25, y = 0.898, label = 'Daytime release', size = 5.1, hjust = 0) +
annotate('segment', xend = 0.206, x = 0.242, y = 0.898, yend = 0.898, arrow = arrow(length = unit(0.1, 'cm'), type = 'closed', angle = 90), size = 0.35) +
annotate('segment', x = 0.382, xend = 0.418, y = 0.898, yend = 0.898, arrow = arrow(length = unit(0.1, 'cm'), type = 'closed', angle = 90), size = 0.35) +
annotate("text", x = 0.6675, y = 0.898, label = 'Nighttime release', size = 5.1, hjust = 0) +
annotate('segment', xend = 0.58, x = 0.6593, y = 0.898, yend = 0.898, arrow = arrow(length = unit(0.1, 'cm'), type = 'closed', angle = 90), size = 0.35) +
annotate('segment', x = 0.8088, xend = 0.888, y = 0.898, yend = 0.898, arrow = arrow(length = unit(0.1, 'cm'), type = 'closed', angle = 90), size = 0.35)

ggsave(file.path(PathFigures,"Figure_5.png"),Fig_5,width = 30, height = 30/2, units="cm")


#############
### Figure 6: Absolute difference in the wind direction between the three downwind UAs (UA-DW) and the upwind UA (UA-UW) recorded during the entire measurement campaign, given as 10 min data. The exact locations of the UAs are given in Fig. 3.

Fig_dWD_10min_woT <- SonicAll_10min_cast[MC == "MC" & top_variable == 'WD' & comparison == 'UW', {
    ggplot(.SD, aes(x=(WD.c + 180) %% 360 - 180, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    scale_x_continuous(limits=c(-10,120),breaks=seq(0,120,by=20)) +
    scale_y_continuous(limits=c(-30,30), breaks=seq(-30,30,by=10)) +
    xlab(expression("WD"["UA-UW"]*" [°]")) +
    facet_grid(Sonic_ordh ~.) +
    ylab(expression("WD"["UA-DW"]*"  -  WD"["UA-UW"]*"  [°]")) + 
    theme_bw(base_size = 18) +
    theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_blank(), strip.background = element_rect(fill="white"), 
	plot.margin = unit(c(0.08, 0.08, 0.08, 0.08), "cm"))
}]

ggsave(file.path(PathFigures, "Figure_6.png"), Fig_dWD_10min_woT, width = 20, height = 20,units="cm")


####################################
####################################
#####                          #####
#####    Tables publication    #####
#####                          #####
####################################
####################################

############
### Table 1: Precision of the OP determined according to Häni et al. (2021). N is the number of intervals used to determine precision.

GFs10min[,'Campaign'] <- NA
GFs10min[IC1,'Campaign'] <- 'IC1'
GFs10min[MC,'Campaign'] <- 'MC'
GFs10min[IC2,'Campaign'] <- 'IC2'

dt_GF <- as.data.table(cbind(st=st(GFs10min),GFs10min))
Conc <- merge(merge(dt_GF[!is.na(Campaign)],Result[Sonic=='SonicC',.(st,Q_MFCex)],all=TRUE),as.data.table(cbind(st=st(WS2_10min),WS2_10min[,c('Press','Temp','WS')])))
Conc[,R := 8.31446261815324] # m3 Pa K-1 mol-1
Conc[,MM := 16.043] # g/mol

Conc[,':=' (dCppm_OP2h = ((GF17-GF26) * R * (Temp+273.15)* 10) / (MM * Press), dCppm_OP5h = ((GF18-GF26) * R * (Temp+273.15)* 10) / (MM * Press),
	 dCppm_OP9h = ((GF16-GF26) * R * (Temp+273.15)* 10) / (MM * Press), dCppm_OP12h = ((GF25-GF26) * R * (Temp+273.15)* 10) / (MM * Press))] 

Conc[,MFC := 1]
Conc[, help_MFC := shift(Q_MFCex, type = "lead")]
Conc[, help2_MFC := shift(Q_MFCex,6, type = "lag")] # remove 1 h
Conc[Q_MFCex > 0 | help_MFC > 0 | help2_MFC > 0, MFC := NA]
Conc[,help_MFC := NULL]
Conc[,help2_MFC := NULL]

Conc_melt <- melt(Conc,id.vars=c('st','Campaign','MFC','Q_MFCex'),measure.vars=c('dCppm_OP2h','dCppm_OP5h','dCppm_OP9h','dCppm_OP12h'),value.name='dCppm')
Conc_melt[!is.finite(dCppm),dCppm := NA_real_] # remove infinite values

OP_Precision <- sapply(Conc_melt[,unique(Campaign)], function(j){
	OPvar <- sapply(Conc_melt[,unique(variable)], function(k){
		out <- Conc_melt[Campaign == j & variable == k, round(2.9*mad(get('dCppm'),na.rm=TRUE) * 110 /sqrt(2),1)]
		names(out) <- k
		return(out)
	})	
})
row.names(OP_Precision) <- c('prec_OP-2.0h','prec_OP-5.3h','prec_OP-8.6h','prec_OP-12h')


OP_Precision_all <- sapply(Conc_melt[,unique(variable)], function(k){
		out <- Conc_melt[variable == k, round(2.9*mad(get('dCppm'),na.rm=TRUE) * 110 /sqrt(2),1)]
		names(out) <- k
		return(out)
	})
names(OP_Precision_all) <- c('prec_OP-2.0h','prec_OP-5.3h','prec_OP-8.6h','prec_OP-12h')


OP_Precision # precision in ppm-m
Conc_melt[!is.na(MFC) & !is.na(dCppm),.(N=.N),by=.(Campaign,variable)] # Number of measurements

OP_Precision_all # precision in ppm-m
Conc_melt[!is.na(MFC) & !is.na(dCppm),.(N=.N),by=.(variable)] # Number of measurements


############
### Table 2: Median recovery rates with standard deviation, median concentration enhancements (1C) with standard deviation, and number of 10 min intervals (N) for all OPs using the data from the UA-UW for the two releases in the MC. Daytime release (unstable atmospheric conditions) and nighttime release (stable atmospheric conditions).

## recovery rates
Result_R <- data.table(rbind(
 	cbind(melt(Result[Sonic_ord2 == "UA-UW"], id.vars=c('st','L','U_sonic','Ustar','WD'),measure.vars=c("R_GF16","R_GF17","R_GF18","R_GF25","R_GF26"), variable.name="R_GF", value.name="Recovery"),Sonic="UA-UW")
 	,cbind(melt(Result[Sonic_ord2 == "UA-2.0h"], id.vars=c('st','L','U_sonic','Ustar','WD'),measure.vars=c("R_GF16","R_GF17","R_GF18","R_GF25","R_GF26"), variable.name="R_GF", value.name="Recovery"),Sonic="UA-50m")
 	,cbind(melt(Result[Sonic_ord2 == "UA-5.3h"], id.vars=c('st','L','U_sonic','Ustar','WD'),measure.vars=c("R_GF16","R_GF17","R_GF18","R_GF25","R_GF26"), variable.name="R_GF", value.name="Recovery"),Sonic="UA-100m")
 	,cbind(melt(Result[Sonic_ord2 == "UA-8.6h"], id.vars=c('st','L','U_sonic','Ustar','WD'),measure.vars=c("R_GF16","R_GF17","R_GF18","R_GF25","R_GF26"), variable.name="R_GF", value.name="Recovery"),Sonic="UA-150m")
	))
Result_R$OP <- factor(Result_R$R_GF, levels=c("R_GF17","R_GF18","R_GF16","R_GF25"),labels=c("OP-2.0h","OP-5.3h","OP-8.6h","OP-12h"))

Result_R[L > 0,stability := 'stable']
Result_R[L < 0,stability := 'unstable']

# revovery rate entire MC
Result_R[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & Sonic=='UA-UW' & !is.na(stability) & !is.na(OP),
.(mean=round(mean(Recovery,na.rm=TRUE),2),median=round(median(Recovery,na.rm=TRUE),2),sd=round(sd(Recovery,na.rm=TRUE),2)),by=.(OP)][order(OP)]

## recovery daytime and nighttime release
Result_R[st > "2021-03-19 08:30" & st < "2021-03-20 08:00" & !is.na(stability) & !is.na(OP),
.(mean=round(mean(Recovery,na.rm=TRUE),2),median=round(median(Recovery,na.rm=TRUE),2),sd=round(sd(Recovery,na.rm=TRUE),2)),by=.(stability,OP)][order(stability,OP)]

## delta C during MC
# only wanna know times when we determined emissions
dt_help <- Result[Sonic == 'SonicC']
dt_help[Q_MFCex > 4,R_OP9h := Q_GF16 / Q_MFCex]
dt_help[Q_MFCex > 4,R_OP2h := Q_GF17 / Q_MFCex]
dt_help[Q_MFCex > 4,R_OP5h := Q_GF18 / Q_MFCex]
dt_help[Q_MFCex > 4,R_OP12h := Q_GF25 / Q_MFCex]

V1 <- merge(Conc_melt[Campaign == 'MC' & variable == 'dCppm_OP2h'],dt_help[,.(st,Rate = R_OP2h,L=1/L)])
V2 <- merge(Conc_melt[Campaign == 'MC' & variable == 'dCppm_OP5h'],dt_help[,.(st,Rate = R_OP5h,L=1/L)])
V3 <- merge(Conc_melt[Campaign == 'MC' & variable == 'dCppm_OP9h'],dt_help[,.(st,Rate = R_OP9h,L=1/L)])
V4 <- merge(Conc_melt[Campaign == 'MC' & variable == 'dCppm_OP12h'],dt_help[,.(st,Rate = R_OP12h,L=1/L)])

dt_Conc <- rbind(Conc_melt[Campaign != 'MC'],V1, V2, V3, V4,fill=TRUE)
dt_Conc[L < 0,stability := 'instable']
dt_Conc[L > 0,stability := 'stable']

# dC and N entire MC
dt_Conc[!is.na(Rate) & Campaign == 'MC',.(mean_ppmm=round(mean(dCppm*100),1),median_ppmm=round(median(dCppm*100),1),sd=round(sd(dCppm*100),1),N=.N),by=variable]
# dC and N nighttime and daytime release
dt_Conc[!is.na(Rate) & Campaign == 'MC',.(mean_ppmm=round(mean(dCppm*100),1),median_ppmm=round(median(dCppm*100),1),sd=round(sd(dCppm*100),1),N=.N),by=.(variable,stability)]


############
### Table 3: Mean wind direction (WD), mean wind speed (WS), mean friction velocity (u), and the mean of the inverse of the Obukhov length (L) recorded by the UA during the two release phases in the MC.

Result[L > 0, stability := 'stable']
Result[L < 0, stability := 'unstable']

Result[!is.na(Sonic_ord) & Campaign == "MC" & Q_MFCex > 4, .(meanWS=round(mean(U_sonic,na.rm=TRUE),1),
	meanWD=round((360 + atan2(mean(U_sonic * sin(WD * pi/180),na.rm=TRUE),mean(U_sonic * cos(WD * pi/180),na.rm=TRUE)) * 180/pi) %% 360,1),
	meanUstar=round(mean(Ustar,na.rm=TRUE),2),meanL=round(mean(1/L,na.rm=TRUE),2)),by=.(stability,Sonic_ord)][order(-stability,Sonic_ord)]


#-----------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------#
					#/////////////////////////////////---    SUPPLEMENT    ---\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#-----------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------#


####################################
####################################
#####                          #####
#####    Figures supplement    #####
#####                          #####
####################################
####################################

#########################
### Supplement - Fig. S1:

limits_recovery_IC2 <- c(parse_date_time3("2021-03-22 09:00",tz='Etc/GMT-1'), parse_date_time3("2021-03-22 15:00",tz='Etc/GMT-1'))

Fig_recovery_IC2_IDM <- Result[st > "2021-03-22 09:30" & st < "2021-03-22 14:15" & Sonic_ord2=='UA-2.0h',{
	ggplot(.SD, aes(x=st)) +
	geom_rect(aes(xmin=parse_date_time3("22.03.2021 09:40",tz="Etc/GMT-1"),xmax=parse_date_time3("22.03.2021 12:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_rect(aes(xmin=parse_date_time3("22.03.2021 12:37",tz="Etc/GMT-1"),xmax=parse_date_time3("22.03.2021 14:11",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_hline(yintercept=1,linetype=3) +
	geom_line(aes(y=R_GF17,colour="OP-2.0h"),size=0.6) +
	geom_point(aes(y=R_GF17,colour="OP-2.0h"),size=1) +
	geom_line(aes(y=R_GF18,colour="OP-5.3h"),size=0.6) +
	geom_point(aes(y=R_GF18,colour="OP-5.3h"),size=1) +
	geom_line(aes(y=R_GF16,colour="OP-8.6h"),size=0.6) +
	geom_point(aes(y=R_GF16,colour="OP-8.6h"),size=1) +
	geom_line(aes(y=R_GF25,colour="OP-12h"),size=0.6) +
	geom_point(aes(y=R_GF25,colour="OP-12h"),size=1) +
	geom_line(aes(y=R_GF26,colour="OP-UW"),size=0.6) +
	geom_point(aes(y=R_GF26,colour="OP-UW"),size=1) +
	xlab(NULL) +
	scale_x_datetime(limits=limits_recovery_IC2,date_labels = "%b %d %H:%M") +
	scale_y_continuous(limits=c(0.2,1),breaks=seq(0.25,1,0.25)) +
	ylab('IDM recovery rate [%]') +
	scale_colour_manual(name=NULL,values=c('OP-2.0h'='#F8766D','OP-5.3h'='#7CAE00','OP-8.6h'='#00BFC4','OP-12h'='#C77CFF','OP-UW'='black'),breaks=c('OP-2.0h','OP-5.3h','OP-8.6h','OP-12h','OP-UW')) +
	theme_bw(base_size = 18) +
	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3),
	strip.background =element_rect(fill="white"), legend.title=element_blank(), legend.position='top')
}]

ggsave(file.path(PathFigures,"Supplement/Supplement - Fig. S1.png"),Fig_recovery_IC2_IDM,width = 30, height = 30/2, units="cm")


#########################
### Supplement - Fig. S2:

Fig_MFC_paper <- MFC_dt[,{
	ggplot(.SD, aes(x=st)) +
	geom_line(aes(y=Q_MFCex,col='Q')) +
	geom_line(aes(y=Temp_MFC*coeff_MFC,col='Temp')) +
	xlab(NULL) +
	geom_hline(yintercept=0,linetype=3) +
	scale_y_continuous(name=expression("Release rate MFC [kg h"^-1*"]"),sec.axis = sec_axis(~./coeff_MFC, name="Temperature MFC [°C]")) +
	scale_colour_manual(name=NULL,values=c('black','grey70')) +
	theme_classic(base_size = 18) +
	theme(legend.position = "none", axis.title.y.right = element_text(color = "grey70"))
}]

Fig_CC_paper <- CC_dt[,{
	ggplot(.SD, aes(x=st)) +
	geom_line(aes(y=Press_CCex,colour="Q")) +
	geom_line(aes(y=Temp_CC*coeff_CC,colour="Temp")) +
	xlab(NULL) +
	geom_hline(yintercept=0,linetype=3) +
	scale_y_continuous(name=expression("Pressure gas bundle [hPa]"),sec.axis = sec_axis(~./coeff_CC, name="Temperature pressure sensor [°C]")) +
	scale_colour_manual(name=NULL,values=c('black','grey70')) +
	theme_classic(base_size = 18) +
	theme(legend.position = "none", axis.title.y.right = element_text(color = "grey70"))
}]


Fig_MFC_Press_paper <- ggarrange(Fig_CC_paper,Fig_MFC_paper,ncol=1,align="hv")

ggsave(Fig_MFC_Press_paper,file=file.path(PathFigures,"/Supplement/Supplement - Fig. S2.png"),width=35,height=35/1.57, units="cm")


#########################
### Supplement - Fig. S3:

{
png(file.path(PathFigures,"/Supplement/Supplement - Fig. S3.png"),width=30,height=10,unit="cm",res=600)
par(mfrow=c(1,3),mar=c(0,0,0,0),oma=c(0.2, 0.2, 0.2, 0.2))
# IC1
plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Source,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
Arrows(x0=Sensors_IC1[2,][,4],y0=Sensors_IC1[2,][,5],
	x1=Sensors_IC1[1,][,4],y1=Sensors_IC1[1,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_IC1[4,][,4],y0=Sensors_IC1[4,][,5],
	x1=Sensors_IC1[3,][,4],y1=Sensors_IC1[3,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_IC1[6,][,4],y0=Sensors_IC1[6,][,5],
	x1=Sensors_IC1[5,][,4],y1=Sensors_IC1[5,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_IC1[8,][,4],y0=Sensors_IC1[8,][,5],
	x1=Sensors_IC1[7,][,4],y1=Sensors_IC1[7,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_IC1[10,][,4],y0=Sensors_IC1[10,][,5],
	x1=Sensors_IC1[9,][,4],y1=Sensors_IC1[9,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_IC1[2,][,4],y0=Sensors_IC1[2,][,5],
	x1=Sensors_IC1[1,][,4],y1=Sensors_IC1[1,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_IC1[4,][,4],y0=Sensors_IC1[4,][,5],
	x1=Sensors_IC1[3,][,4],y1=Sensors_IC1[3,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_IC1[6,][,4],y0=Sensors_IC1[6,][,5],
	x1=Sensors_IC1[5,][,4],y1=Sensors_IC1[5,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_IC1[8,][,4],y0=Sensors_IC1[8,][,5],
	x1=Sensors_IC1[7,][,4],y1=Sensors_IC1[7,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_IC1[10,][,4],y0=Sensors_IC1[10,][,5],
	x1=Sensors_IC1[9,][,4],y1=Sensors_IC1[9,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
plot(Tree,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
plot(Sensors_IC1[Sensors_IC1[,1] == 'Sonic2',],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black", lwd=3),add=TRUE)
plot(Sensors_MC[Sensors_MC[,1] == "WS2",],sensors.text.args=list(labels=""),points.args = list(pch = 17, cex = 2, col ="black"),add=TRUE)
text(x=centre[1]-185,y=centre[2]+120,"IC1",cex=1.5)

# MC
plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
# abline(v=seq(-400,length.out=10,by=100),h=seq(-400,length.out=10,by=100),lty=3,col="lightgrey",lwd=1.5)
plot(Source,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[9,][,4],y0=Sensors_MC[9,][,5],
	x1=Sensors_MC[10,][,4],y1=Sensors_MC[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,code=1)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[9,][,4],y0=Sensors_MC[9,][,5],
	x1=Sensors_MC[10,][,4],y1=Sensors_MC[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)

plot(Tree,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
plot(Sonics_MC,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
plot(Sensors_MC[Sensors_MC[,1] == "WS2",],sensors.text.args=list(labels=""),points.args = list(pch = 17, cex = 2, col ="black"),add=TRUE)
bx <- centre[1]+30
by <- centre[2]-205
dx1 <- 50
dx2 <- 100
dy <- 4
polygon(c(bx,bx+dx1,bx+dx1,bx),c(by,by,by-dy,by-dy),col="black")
polygon(c(bx+dx1,bx+dx2,bx+dx2,bx+dx1),c(by,by,by-dy,by-dy),col="white",border="black")
lines(c(bx,bx),c(by,by+dy*(2/3)),col="black")
lines(c(bx+dx1,bx+dx1),c(by,by+dy*(1/3)),col="black")
lines(c(bx+dx2,bx+dx2),c(by,by+dy*(2/3)),col="black")
text(bx,by+dy*3,"0 m",cex=1.3,adj=c(0.5,0.5),col="black")
text(bx+dx1,by+dy*3,"50 m",cex=1.3,adj=c(0.5,0.5),col="black")
text(bx+dx2,by+dy*3,"100 m",cex=1.3,adj=c(0.5,0.5),col="black")
## Direction north
wr.x <- bx+85
wr.y <- by+40  
polygon(c(wr.x,wr.x-5,wr.x),c(wr.y,wr.y-5,wr.y+40),border=1)
polygon(c(wr.x,wr.x+5,wr.x),c(wr.y,wr.y-5,wr.y+40),border=1,col=1)
text(x=wr.x,y=wr.y+50,labels="N",cex=1.5)
text(x=centre[1]-185,y=centre[2]+120,"MC",cex=1.5)

# IC2
plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Source,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
Arrows(x0=Sensors_IC2[1,][,4],y0=Sensors_IC2[1,][,5],
	x1=Sensors_IC2[2,][,4],y1=Sensors_IC2[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_IC2[3,][,4],y0=Sensors_IC2[3,][,5],
	x1=Sensors_IC2[4,][,4],y1=Sensors_IC2[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_IC2[5,][,4],y0=Sensors_IC2[5,][,5],
	x1=Sensors_IC2[6,][,4],y1=Sensors_IC2[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_IC2[7,][,4],y0=Sensors_IC2[7,][,5],
	x1=Sensors_IC2[8,][,4],y1=Sensors_IC2[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_IC2[9,][,4],y0=Sensors_IC2[9,][,5],
	x1=Sensors_IC2[10,][,4],y1=Sensors_IC2[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_IC2[1,][,4],y0=Sensors_IC2[1,][,5],
	x1=Sensors_IC2[2,][,4],y1=Sensors_IC2[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_IC2[3,][,4],y0=Sensors_IC2[3,][,5],
	x1=Sensors_IC2[4,][,4],y1=Sensors_IC2[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_IC2[5,][,4],y0=Sensors_IC2[5,][,5],
	x1=Sensors_IC2[6,][,4],y1=Sensors_IC2[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_IC2[7,][,4],y0=Sensors_IC2[7,][,5],
	x1=Sensors_IC2[8,][,4],y1=Sensors_IC2[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_IC2[9,][,4],y0=Sensors_IC2[9,][,5],
	x1=Sensors_IC2[10,][,4],y1=Sensors_IC2[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
plot(Tree,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
plot(Sensors_IC2[Sensors_IC2[,1] == 'Sonic2',] ,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black", lwd=3),add=TRUE)
plot(Sensors_MC[Sensors_MC[,1] == "WS2",],sensors.text.args=list(labels=""),points.args = list(pch = 17, cex = 2, col ="black"),add=TRUE)
text(x=centre[1]-185,y=centre[2]+120,"IC2",cex=1.5)

## Legend
xl <- centre[1] + 58
yl <- centre[2] - 152
dyl <- 16
dxl <- 8
Arrows(x0=xl-dxl,x1=xl+dxl,y0=yl,y1=yl,arr.type="triangle",arr.length=0.1,arr.adj=1)
Arrows(x0=xl+dxl,x1=xl-dxl,y0=yl,y1=yl,arr.type="T",arr.length=0.1,arr.adj=1)
points(x=c(xl,xl),y=c(yl-dyl,yl-dyl*2),pch=c(8,17),cex=1)
points(x=xl,y=yl-dyl*3,pch=21,cex=2,col='grey60',bg='grey95')
rect(xleft=xl-dxl,xright=xl+dxl,ybottom=yl-dyl*4-4,ytop=yl-dyl*4+4,col="grey95")
text(x=rep(xl+dxl*2,5),y=c(yl,yl-dyl,yl-dyl*2,yl-dyl*3,yl-dyl*4),labels=c("OP","UA","Weather station","Tree","Barn"),pos=4)
lines(x=c(xl-dxl*2,xl-dxl*2,xl+200),y=c(yl-100,yl+dyl,yl+dyl))
lines(x=c(xl-20,xl+100,xl+100),y=c(yl-80,yl-80,yl+20))

dev.off()
}


#########################
### Supplement - Fig. S4:

Fig_dUsonic_WD_10min <- SonicAll_10min_cast[MC == "MC" & U_sonic.c > 0.0 & top_variable == 'U_sonic' & comparison == 'UW', {
    ggplot(.SD, aes(x=(WD.c + 180) %% 360 - 180, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    scale_x_continuous(limits=c(-10,120),breaks=seq(0,120,by=20)) +
    scale_y_continuous(limits=c(-2,2), breaks=round(seq(-2,2,by=1),0)) +
    xlab(expression("WD"["UA-UW"]*" [°]")) +
    facet_grid(Sonic_ordh ~.) +
    ylab(expression(paste(italic('u')['UA-DW']*'  -  ',italic('u')['UA-UW']*'  [m s'^-1*']'))) + 
    theme_bw(base_size = 18) +
    theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_blank(), strip.background = element_rect(fill="white"))
}]


Fig_dUstar_WD_10min <- SonicAll_10min_cast[MC == "MC" & U_sonic.c > 0.0 & top_variable == 'Ustar' & comparison == 'UW', {
    ggplot(.SD, aes(x=(WD.c + 180) %% 360 - 180, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    scale_x_continuous(limits=c(-10,120),breaks=seq(0,120,by=20)) +
    scale_y_continuous(limits=c(-0.3,0.3), breaks=round(seq(-0.3,0.3,by=0.1),1)) +
    xlab(expression("WD"["UA-UW"]*" [°]")) +
    facet_grid(Sonic_ordh ~.) +
    ylab(expression(paste(italic('u'['*'])['UA-DW']*'  -  ',italic('u'['*'])['UA-UW']*'  [m s'^-1*']'))) + 
    theme_bw(base_size = 18) +
    theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_blank(), strip.background = element_rect(fill="white"))
}]


Fig_dsVu_WD_10min <- SonicAll_10min_cast[MC == "MC" & U_sonic.c > 0.0 & top_variable == 'sVu' & comparison == 'UW', {
    ggplot(.SD, aes(x=(WD.c + 180) %% 360 - 180, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    scale_x_continuous(limits=c(-10,120),breaks=seq(0,120,by=20)) +
    scale_y_continuous(limits=c(-6,6), breaks=round(seq(-6,6,by=2),0)) +
    xlab(expression("WD"["UA-UW"]*" [°]")) +
    facet_grid(Sonic_ordh ~.) +
  	ylab(expression(italic(sigma[v])*italic(u['*'])^-1*phantom()['UA-DW']*'  -  '*italic(sigma[v])*italic(u['*'])^-1*phantom()['UA-UW'])) +
    theme_bw(base_size = 18) +
    theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_blank(), strip.background = element_rect(fill="white"))
}]


Fig_dsWu_WD_10min <- SonicAll_10min_cast[MC == "MC" & U_sonic.c > 0.0 & top_variable == 'sWu' & comparison == 'UW', {
    ggplot(.SD, aes(x=(WD.c + 180) %% 360 - 180, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    scale_x_continuous(limits=c(-10,120),breaks=seq(0,120,by=20)) +
    scale_y_continuous(limits=c(-4,4), breaks=round(seq(-6,6,by=2),0)) +
    xlab(expression("WD"["UA-UW"]*" [°]")) +
    facet_grid(Sonic_ordh ~.) +
   	ylab(expression(italic(sigma[w])*italic(u['*'])^-1*phantom()['UA-DW']*'  -  '*italic(sigma[w])*italic(u['*'])^-1*phantom()['UA-UW'])) +
    theme_bw(base_size = 18) +
    theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_blank(), strip.background = element_rect(fill="white"))
}]


Fig_turbulence <- ggarrange(Fig_dUsonic_WD_10min,Fig_dUstar_WD_10min,Fig_dsVu_WD_10min,Fig_dsWu_WD_10min,ncol=2,nrow=2,align="hv")
ggsave(file.path(PathFigures,'/Supplement/Supplement - Fig. S4.png'),Fig_turbulence,width=30, height=30, units='cm')


#########################
### Supplement - Fig. S5:

Fig_recovery_disturbedSonic_IDM <- Result[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & !is.na(Sonic_ord2) & Sonic_ord2 != 'UA-UW',{
	ggplot(.SD, aes(x=st)) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 11:19",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 11:40",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 01:08",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("20.03.2021 01:21",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_hline(yintercept=1,linetype=3) +
	geom_line(aes(y=R_GF17,colour="OP-2.0h"),size=0.6) +
	geom_point(aes(y=R_GF17,colour="OP-2.0h"),size=1) +
	geom_line(aes(y=R_GF18,colour="OP-5.3h"),size=0.6) +
	geom_point(aes(y=R_GF18,colour="OP-5.3h"),size=1) +
	geom_line(aes(y=R_GF16,colour="OP-8.6h"),size=0.6) +
	geom_point(aes(y=R_GF16,colour="OP-8.6h"),size=1) +
	geom_line(aes(y=R_GF25,colour="OP-12h"),size=0.6) +
	geom_point(aes(y=R_GF25,colour="OP-12h"),size=1) +
	xlab(NULL) +
	ylab('IDM recovery rate [%]') +
	facet_grid(Sonic_ord2~.) +
	scale_colour_manual(name=NULL,values=c('OP-2.0h'='#F8766D','OP-5.3h'='#7CAE00','OP-8.6h'='#00BFC4','OP-12h'='#C77CFF'),breaks=c('OP-2.0h','OP-5.3h','OP-8.6h','OP-12h')) +
	theme_bw(base_size = 18) +
	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3), strip.background =element_rect(fill="white"),
	legend.title=element_blank(), legend.position='top')
    # theme(legend.position = "right", text = element_text(size=30))
}]

ggsave(file.path(PathFigures,"/Supplement/Supplement - Fig. S5.png"),Fig_recovery_disturbedSonic_IDM,width = 30, height = 30/1.5, units="cm")


#########################
### Supplement - Fig. S6:

{
png(file.path(PathFigures,"/Supplement/Supplement - Fig. 6.png"),width=30,height=30,unit="cm",res=600)

indMK_MFC_1 <- parse_date_time3("19.03.2021 10:30",tz="Etc/GMT-1") + c(0:37)*600
indMK_MFC_2 <- parse_date_time3("19.03.2021 21:50",tz="Etc/GMT-1") + c(0:54)*600
indMK_MFC_3 <- c(indMK_MFC_1,indMK_MFC_2)

plot_range <- data.frame(x=c(583735,583735+520),y=c(210050,210050+520))
bx <- plot_range[2,1]-220
by <- plot_range[1,2]+10
dx1 <- 100
dx2 <- 200
dy <- 4
wr.x <- bx+185
wr.y <- by+50

par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(0.2,0.2,0.2,0.2))
## SonicC
plot(type="n",plot_range,xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
windrose(SonicC[indMK_MFC_3], wd = "WD", ws = "U_sonic", max_freq = 10, draw.grid = FALSE,
    ws_breaks = c(0,2,3,4,7),delta_wd = 6, add = TRUE, scale = 0.5,show.legend=FALSE, start=0, center = c(583952.8,210336.9))
plot(Source,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Tree,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[9,][,4],y0=Sensors_MC[9,][,5],
	x1=Sensors_MC[10,][,4],y1=Sensors_MC[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,code=1,col='grey60')
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[9,][,4],y0=Sensors_MC[9,][,5],
	x1=Sensors_MC[10,][,4],y1=Sensors_MC[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
plot(Sensors_MC[Sensors_MC[1]=='SonicC',],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
text(x=plot_range[1,1]-20,y=plot_range[2,2],"UA-UW",cex=2,pos=4)

## SonicA
plot(type="n",plot_range,xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
windrose(SonicA[indMK_MFC_3], wd = "WD", ws = "U_sonic", max_freq = 10, draw.grid = FALSE,
    ws_breaks = c(0,2,3,4,7),delta_wd = 6, add = TRUE, scale = 0.5,show.legend=FALSE, start=0, center = c(583862.8,210199.9))
plot(Source,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Tree,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[9,][,4],y0=Sensors_MC[9,][,5],
	x1=Sensors_MC[10,][,4],y1=Sensors_MC[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60',code=1)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[9,][,4],y0=Sensors_MC[9,][,5],
	x1=Sensors_MC[10,][,4],y1=Sensors_MC[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
plot(Sensors_MC[Sensors_MC[1] == 'SonicA',],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
text(x=plot_range[1,1]-20,y=plot_range[2,2],"UA-2.0h",cex=2,pos=4)
## Legende
dyl <- 30
dxl <- 15
xl <- plot_range[2,1] - 90
yl <- plot_range[1,2] + 130 - dyl
Arrows(x0=xl-dxl,x1=xl+dxl,y0=yl,y1=yl,arr.type="triangle",arr.length=0.15,lwd=1.5,arr.adj=1)
Arrows(x0=xl+dxl,x1=xl-dxl,y0=yl,y1=yl,arr.type="T",arr.length=0.15,lwd=1.5,arr.adj=1)
points(x=xl,y=yl-dyl,pch=8,cex=1.5)
points(x=xl,y=yl-dyl,pch=1,cex=1.7)
points(x=xl,y=yl-dyl*2,pch=21,cex=2.3,col='grey60',bg='grey95')
rect(xleft=xl-dxl,xright=xl+dxl,ybottom=yl-dyl*3-8,ytop=yl-dyl*3+8,col="grey95")
text(x=rep(xl+dxl*2,4),y=c(yl,yl-dyl,yl-dyl*2,yl-dyl*3),labels=c("OP","UA","Tree","Barn"),pos=4,cex=1.5)
lines(x=c(xl-dxl*2,xl-dxl*2,xl+200),y=c(yl-200,yl+dyl,yl+dyl))

## Sonic2 (SonicD)
plot(type="n",plot_range,xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
windrose(Sonic2[indMK_MFC_3], wd = "WD", ws = "U_sonic", max_freq = 10, draw.grid = FALSE,
    ws_breaks = c(0,2,3,4,7),delta_wd = 6, add = TRUE, scale = 0.5,
    unit=expression("u [m s"^-1*"]"), show.legend = FALSE, start=0, center = c(583828.7,210164.5))
plot(Source,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Tree,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[9,][,4],y0=Sensors_MC[9,][,5],
	x1=Sensors_MC[10,][,4],y1=Sensors_MC[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60',code=1)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[9,][,4],y0=Sensors_MC[9,][,5],
	x1=Sensors_MC[10,][,4],y1=Sensors_MC[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
plot(Sensors_MC[Sensors_MC[1]=='Sonic2',],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
text(x=plot_range[1,1]-20,y=plot_range[2,2],"UA-5.3h",cex=2,pos=4)
## Massstab
polygon(c(bx,bx+dx1,bx+dx1,bx),c(by,by,by-dy,by-dy),col="black")
polygon(c(bx+dx1,bx+dx2,bx+dx2,bx+dx1),c(by,by,by-dy,by-dy),col="white",border="black")
lines(c(bx,bx),c(by,by+dy*(2/3)),col="black")
lines(c(bx+dx1,bx+dx1),c(by,by+dy*(1/3)),col="black")
lines(c(bx+dx2,bx+dx2),c(by,by+dy*(2/3)),col="black")
text(bx,by+dy*3,"0 m",cex=1.6,adj=c(0.5,0.4),col="black")
text(bx+dx1,by+dy*3,"100 m",cex=1.6,adj=c(0.5,0.4),col="black")
text(bx+dx2,by+dy*3,"200 m",cex=1.6,adj=c(0.5,0.4),col="black")
## Richtung Nord
polygon(c(wr.x,wr.x-5,wr.x),c(wr.y,wr.y-5,wr.y+40),border=1)
polygon(c(wr.x,wr.x+5,wr.x),c(wr.y,wr.y-5,wr.y+40),border=1,col=1)
text(x=wr.x,y=wr.y+60,labels="N",cex=2)

## SonicB
plot(type="n",plot_range,xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
windrose(SonicB[indMK_MFC_3], wd = "WD", ws = "U_sonic", max_freq = 10, draw.grid = FALSE,
    ws_breaks = c(0,2,3,4,7),delta_wd = 6, add = TRUE, scale = 0.5, legend.args = list(title = expression("u [m s"^-1*"]"), x = "bottomright", cex=1.5, bty="n"), start=0, center=c( 583792.4,210128.0))
plot(Source,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Tree,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[9,][,4],y0=Sensors_MC[9,][,5],
	x1=Sensors_MC[10,][,4],y1=Sensors_MC[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60',code=1)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MC[9,][,4],y0=Sensors_MC[9,][,5],
	x1=Sensors_MC[10,][,4],y1=Sensors_MC[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
plot(Sensors_MC[Sensors_MC[1]=='SonicB',],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
lines(x=c(bx+140,bx+140,bx+400),y=c(by-300,by+130,by+130))
text(x=plot_range[1,1]-20,y=plot_range[2,2],"UA-8.6h",cex=2,pos=4)

dev.off()
}


###################################
###################################
#####                         #####
#####    Tables supplement    #####
#####                         #####
###################################
###################################

# no tables that need code. For Table S1 see Apply filter script


#-----------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------#
					#/////////////////////////////////---    INITIAL SUBMISSION    ---\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#-----------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------#


########################################
########################################
#####                              #####
#####    Figures initial submit    #####
#####                              #####
########################################
########################################

##############################
### Figure 1 - initial submit:

# photo


##############################
### Figure 2 - initial submit:

# schematic


##############################
### Figure 3 - initial submit:

{
png(file.path(PathFigures,"Initial submit/Initial submit - Fig. 3.png"),width=23,height=15,unit="cm",res=600)
par(mfrow=c(1,1),mar=c(0.1,0.1,0.1,0.1))
# MK
plot(type="n",1, xlim=c(centre[1]-255,centre[1]+200),ylim=c(centre[2]-215,centre[2]+97),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Source,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[9,][,4],y0=Sensors_MC[9,][,5],
	x1=Sensors_MC[10,][,4],y1=Sensors_MC[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,code=1)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[9,][,4],y0=Sensors_MC[9,][,5],
	x1=Sensors_MC[10,][,4],y1=Sensors_MC[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
plot(Tree,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
plot(Sonics_MC,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
plot(Sensors_MC[Sensors_MC[,1] == "WS2",],sensors.text.args=list(labels=""),points.args = list(pch = 17, cex = 2, col ="black"),add=TRUE)
bx <- centre[1]+80
by <- centre[2]-205
dx1 <- 50
dx2 <- 100
dy <- 4
## scale bar
polygon(c(bx,bx+dx1,bx+dx1,bx),c(by,by,by-dy,by-dy),col="black")
polygon(c(bx+dx1,bx+dx2,bx+dx2,bx+dx1),c(by,by,by-dy,by-dy),col="white",border="black")
lines(c(bx,bx),c(by,by+dy*(2/3)),col="black")
lines(c(bx+dx1,bx+dx1),c(by,by+dy*(1/3)),col="black")
lines(c(bx+dx2,bx+dx2),c(by,by+dy*(2/3)),col="black")
text(bx,by+dy*3,"0 m",cex=1.3,adj=c(0.5,0.5),col="black")
text(bx+dx1,by+dy*3,"50 m",cex=1.3,adj=c(0.5,0.5),col="black")
text(bx+dx2,by+dy*3,"100 m",cex=1.3,adj=c(0.5,0.5),col="black")
## North direction
wr.x <- bx+85
wr.y <- by+40  
polygon(c(wr.x,wr.x-5,wr.x),c(wr.y,wr.y-5,wr.y+40),border=1)
polygon(c(wr.x,wr.x+5,wr.x),c(wr.y,wr.y-5,wr.y+40),border=1,col=1)
text(x=wr.x,y=wr.y+50,labels="N",cex=1.5)
# text(x=centre[1]-185,y=centre[2]+120,"MC",cex=1.5)
## label GFs
text(x=Sensors_MC[1,][,4]+21,y=Sensors_MC[1,][,5]-4,"OP-150m")
text(x=Sensors_MC[3,][,4]+21,y=Sensors_MC[3,][,5]-4,"OP-50m")
text(x=Sensors_MC[5,][,4]+21,y=Sensors_MC[5,][,5]-4,"OP-100m")
text(x=Sensors_MC[7,][,4]+21,y=Sensors_MC[7,][,5]-4,"OP-200m")
text(x=Sensors_MC[10,][,4]-21,y=Sensors_MC[10,][,5],"OP-UW")
## label Sonics
text(x=Sonics_MC[1,][,4]-16,y=Sonics_MC[1,][,5]-15,"UA-100m")
text(x=Sonics_MC[2,][,4]-15,y=Sonics_MC[2,][,5]-15,"UA-50m")
text(x=Sonics_MC[3,][,4]-15,y=Sonics_MC[3,][,5]-15,"UA-150m")
text(x=Sonics_MC[4,][,4],y=Sonics_MC[4,][,5]-15,"UA-UW")

## Legend
xl <- centre[1] - 270
yl <- centre[2] + 97
dyl <- 16
dxl <- 8
Arrows(x0=xl-dxl,x1=xl+dxl,y0=yl,y1=yl,arr.type="triangle",arr.length=0.1,arr.adj=1)
Arrows(x0=xl+dxl,x1=xl-dxl,y0=yl,y1=yl,arr.type="T",arr.length=0.1,arr.adj=1)
points(x=c(xl,xl),y=c(yl-dyl,yl-dyl*2),pch=c(8,17),cex=1)
points(x=xl,y=yl-dyl*3,pch=21,cex=2,col='grey60',bg='grey95')
rect(xleft=xl-dxl,xright=xl+dxl,ybottom=yl-dyl*4-4,ytop=yl-dyl*4+4,col="grey95")
text(x=rep(xl+dxl*2,5),y=c(yl,yl-dyl,yl-dyl*2,yl-dyl*3,yl-dyl*4),labels=c("OP","UA","Weather station","Tree","Barn"),pos=4)
lines(x=c(xl-20,xl+100,xl+100),y=c(yl-80,yl-80,yl+20))

dev.off()
}


##############################
### Figure 4 - initial submit:

limits_WSMK  <- parse_date_time3(c("18.03.2021 11:00","21.03.2021 14:00"),tz="Etc/GMT-1")

WS_TempMK <- WS2_dt[st >= parse_date_time3("18.03.2021 11:00",tz='Etc/GMT-1') & st <= parse_date_time3('21.03.2021 14:00',tz='Etc/GMT-1'),{
 	ggplot(.SD,aes(x=st,y=Temp)) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_hline(yintercept = 0, colour="grey") +
 	geom_line() +
 	scale_x_datetime(limits=limits_WSMK,expand=c(0.02,0.02)) +
 	ylim(-5,10) +
 	ylab("Temperature [°C]") +
 	xlab(NULL) +
 	theme_bw(base_size=18) +
 	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3))
 }]

WS_WDMK <- WS2_dt[st >= parse_date_time3("18.03.2021 11:00",tz='Etc/GMT-1') & st <= parse_date_time3('21.03.2021 14:00',tz='Etc/GMT-1'),{
 	ggplot(.SD,aes(x=st,y=WD_WS)) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_point(size=0.5) +
 	scale_x_datetime(limits=limits_WSMK,expand=c(0.02,0.02)) +
	scale_y_continuous(breaks = seq(0, 360, by=90), limits=c(0,360)) +
 	ylab("Wind direction [°]") +
 	xlab(NULL) +
 	theme_bw(base_size=18) +
 	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3))
 }]

WS_WSMK <- WS2_dt[st >= parse_date_time3("18.03.2021 11:00",tz='Etc/GMT-1') & st <= parse_date_time3('21.03.2021 14:00',tz='Etc/GMT-1'),{
 	ggplot(.SD,aes(x=st,y=WS)) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_line() +
 	scale_x_datetime(limits=limits_WSMK,expand=c(0.02,0.02)) +
 	ylab(expression("Wind speed [m s"^-1*"]")) +
 	xlab(NULL) +
 	theme_bw(base_size=18) +
 	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3))
 }]

Meteo_MK <- ggarrange(WS_TempMK,WS_WDMK,WS_WSMK,ncol=1,align="hv")
ggsave(file.path(PathFigures, "Initial submit/Initial submit - Fig. 4.png"), Meteo_MK, width = 30, height = 30/1.57, units="cm")


#############
### Figure 5 - initial submit:

Fig_5_inSub <- SonicAll_10min_cast[MC == "MC" & top_variable == 'WD' & comparison == 'UW', {
    ggplot(.SD, aes(x=(WD.c + 180) %% 360 - 180, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    scale_x_continuous(limits=c(-10,120),breaks=seq(0,120,by=20)) +
    scale_y_continuous(limits=c(-30,30), breaks=seq(-30,30,by=10)) +
    xlab(expression("WD"["UA-UW"]*" [°]")) +
    facet_grid(Sonic_ord ~.) +
    ylab(expression("WD"["UA-DW"]*"  -  WD"["UA-UW"]*"  [°]")) + 
    theme_bw(base_size = 18) +
    theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_blank(), strip.background = element_rect(fill="white"), 
	plot.margin = unit(c(0.08, 0.08, 0.08, 0.08), "cm"))
}]

ggsave(file.path(PathFigures, "Initial submit/Initial submit - Fig. 5.png"), Fig_5_inSub, width = 20, height = 20,units="cm")


## calculate the mean offset for sVu and sWu for the 10min intervals
SonicAll_10min_cast[MC == "MC" & top_variable == 'sVu' & comparison == 'UW', mean(value,na.rm=TRUE), by = Sonic]


##############################
### Figure 6 - initial submit:

Fig_6_inSub <- Result[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & !is.na(Sonic_ord),{
	ggplot(.SD, aes(x=st)) +
	geom_hline(yintercept=1,linetype=3) +
	geom_line(aes(y=R_GF16,colour="OP-150m"),size=0.6) +
	geom_line(aes(y=R_GF17,colour="OP-50m"),size=0.6) +
	geom_line(aes(y=R_GF18,colour="OP-100m"),size=0.6) +
	geom_line(aes(y=R_GF25,colour="OP-200m"),size=0.6) +
	geom_point(aes(y=R_GF16,colour="OP-150m"),size=1) +
	geom_point(aes(y=R_GF17,colour="OP-50m"),size=1) +
	geom_point(aes(y=R_GF18,colour="OP-100m"),size=1) +
	geom_point(aes(y=R_GF25,colour="OP-200m"),size=1) +
	xlab(NULL) +
	ylab('IDM recovery rate [%]') +
	facet_grid(Sonic_ord ~.) +
	scale_colour_manual(name=NULL,values=c('OP-50m'='#F8766D','OP-100m'='#7CAE00','OP-150m'='#00BFC4','OP-200m'='#C77CFF'),breaks=c('OP-50m','OP-100m','OP-150m','OP-200m')) +
	theme_bw(base_size = 18) + theme(strip.background =element_rect(fill="white"))
}]

ggsave(file.path(PathFigures,"Initial submit/Initial submit - Fig. 6.png"),Fig_6_inSub,width = 30, height = 30/1.1, units="cm")


##############################
### Figure 7 - initial submit:

mgs_to_kgh <- (3600)/1E6

XY_SonicA_Conc <- XY_SonicA
XY_SonicA_Conc[[1]]$z <- XY_SonicA_Conc[[1]]$z * 6/mgs_to_kgh/681.43 # Area of the source (681.43 m2)
XY_SonicB_Conc <- XY_SonicB
XY_SonicB_Conc[[1]]$z <- XY_SonicB_Conc[[1]]$z * 6/mgs_to_kgh/681.43
XY_SonicC_Conc <- XY_SonicC
XY_SonicC_Conc[[1]]$z <- XY_SonicC_Conc[[1]]$z * 6/mgs_to_kgh/681.43
XY_Sonic2_Conc <- XY_Sonic2
XY_Sonic2_Conc[[1]]$z <- XY_Sonic2_Conc[[1]]$z * 6/mgs_to_kgh/681.43


{
png(file.path(PathFigures,"/Initial submit/Initial submit - Fig. 7.png"),width=30,height=30,unit="cm",res=600)

bx <- centre[1]-55
by <- centre[2]-230
dx1 <- 50
dx2 <- 100
dy <- 4
wr.x <- bx+85
wr.y <- by+40
dyl <- 16
dxl <- 8
xl <- centre[1] - 10
yl <- centre[2] - 172 - dyl

par(mfrow=c(2,2),mar=c(0,3.7,3.7,0),oma=c(0,0,0,0))
## SonicC
plot(type="n",1, xlim=c(centre[1]-220,centre[1]+60),ylim=c(centre[2]-220,centre[2]+20),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Source,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Tree,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
plot(Sonics_MC[Sonics_MC[,1] == "SonicC",],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
# contourXY(XY_SonicC_Conc, fill = TRUE,add=TRUE,showLegend=TRUE,lpos="bottomright")
contourXY(XY_SonicC_Conc, fill = TRUE,add=TRUE,showLegend=FALSE)
# text(x=centre[1]-180,y=centre[2]+35,"3DUA-UW",cex=2)
text(x=centre[1]-190,y=centre[2]+35,"UA-UW",cex=2)

## SonicA
par(mar=c(0,0,3.7,3.7))
plot(type="n",1, xlim=c(centre[1]-220,centre[1]+60),ylim=c(centre[2]-220,centre[2]+20),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Source,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Tree,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
plot(Sonics_MC[Sonics_MC[,1] == "SonicA",],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
contourXY(XY_SonicA_Conc, fill = TRUE,add=TRUE,showLegend=FALSE)
# text(x=centre[1]-180,y=centre[2]+35,"3DUA-50m",cex=2)
text(x=centre[1]-190,y=centre[2]+35,"UA-50m",cex=2)
## Legende
Arrows(x0=xl-dxl,x1=xl+dxl,y0=yl,y1=yl,arr.type="triangle",arr.length=0.15,lwd=1.5,arr.adj=1)
Arrows(x0=xl+dxl,x1=xl-dxl,y0=yl,y1=yl,arr.type="T",arr.length=0.15,lwd=1.5,arr.adj=1)
points(x=xl,y=yl-dyl,pch=8,cex=1.5)
points(x=xl,y=yl-dyl*2,pch=21,cex=2.3,col='grey60',bg='grey95')
rect(xleft=xl-dxl,xright=xl+dxl,ybottom=yl-dyl*3-4,ytop=yl-dyl*3+4,col="grey95")
text(x=rep(xl+dxl*2,4),y=c(yl,yl-dyl,yl-dyl*2,yl-dyl*3),labels=c("OP","UA","Tree","Barn"),pos=4,cex=1.5)
lines(x=c(xl-dxl*2,xl-dxl*2,xl+200),y=c(yl-100,yl+dyl,yl+dyl))

## Sonic2
par(mar=c(3.7,3.7,0,0))
plot(type="n",1, xlim=c(centre[1]-220,centre[1]+60),ylim=c(centre[2]-220,centre[2]+20),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Source,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Tree,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
plot(Sonics_MC[Sonics_MC[,1] == "Sonic2",],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
contourXY(XY_Sonic2_Conc, fill = TRUE,add=TRUE,showLegend=FALSE)
# text(x=centre[1]-180,y=centre[2]+35,"3DUA-100m",cex=2)
text(x=centre[1]-190,y=centre[2]+35,"UA-100m",cex=2)
# Massstab
polygon(c(bx,bx+dx1,bx+dx1,bx),c(by,by,by-dy,by-dy),col="black")
polygon(c(bx+dx1,bx+dx2,bx+dx2,bx+dx1),c(by,by,by-dy,by-dy),col="white",border="black")
lines(c(bx,bx),c(by,by+dy*(2/3)),col="black")
lines(c(bx+dx1,bx+dx1),c(by,by+dy*(1/3)),col="black")
lines(c(bx+dx2,bx+dx2),c(by,by+dy*(2/3)),col="black")
text(bx,by+dy*3,"0 m",cex=1.6,adj=c(0.5,0.4),col="black")
text(bx+dx1,by+dy*3,"50 m",cex=1.6,adj=c(0.5,0.4),col="black")
text(bx+dx2,by+dy*3,"100 m",cex=1.6,adj=c(0.5,0.4),col="black")
## Richtung Nord
polygon(c(wr.x,wr.x-5,wr.x),c(wr.y,wr.y-5,wr.y+40),border=1)
polygon(c(wr.x,wr.x+5,wr.x),c(wr.y,wr.y-5,wr.y+40),border=1,col=1)
text(x=wr.x,y=wr.y+50,labels="N",cex=2)

## SonicB
par(mar=c(3.7,0,0,3.7))
plot(type="n",1, xlim=c(centre[1]-220,centre[1]+60),ylim=c(centre[2]-220,centre[2]+20),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Source,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Tree,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MC[1,][,4],y0=Sensors_MC[1,][,5],
	x1=Sensors_MC[2,][,4],y1=Sensors_MC[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[3,][,4],y0=Sensors_MC[3,][,5],
	x1=Sensors_MC[4,][,4],y1=Sensors_MC[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[5,][,4],y0=Sensors_MC[5,][,5],
	x1=Sensors_MC[6,][,4],y1=Sensors_MC[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MC[7,][,4],y0=Sensors_MC[7,][,5],
	x1=Sensors_MC[8,][,4],y1=Sensors_MC[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
plot(Sonics_MC[Sonics_MC[,1] == "SonicB",],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
contourXY(XY_SonicB_Conc, fill = TRUE,add=TRUE,showLegend=FALSE)
# text(x=centre[1]-180,y=centre[2]+35,"3DUA-150m",cex=2)
text(x=centre[1]-190,y=centre[2]+35,"UA-150m",cex=2)
# legend
legend("bottomright",legend=c("(0.06,0.3]","(0.3,0.6]","(0.6,3]","(3,5]",">5"),fill=paste0(ConcPalette(5),"4c"),title=expression("Conc. [mg m"^-3*"]"),cex=1.5,bty="n")
lines(x=c(bx+35,bx+35,bx+200),y=c(by-100,by+85,by+85))

dev.off()
}


##############################
### Figure 8 - initial submit:

WDmCols_10 <- colorRampPalette(c("#00ff00", "#003300"))(10)
WDpCols_10 <- colorRampPalette(c("#00ffff", "#0033ff"))(10)
WDvarCols_10 <- c(WDmCols_10,"red",WDpCols_10)
names(WDvarCols_10) <- paste0("SonicB",c(paste0("_m",1:10),"_0",paste0("_p",1:10)))
	
Result_var[Q_MFCex > 4,R_GF16 := Q_GF16 / Q_MFCex]
Result_var[Q_MFCex > 4,R_GF17 := Q_GF17 / Q_MFCex]
Result_var[Q_MFCex > 4,R_GF18 := Q_GF18 / Q_MFCex]
Result_var[Q_MFCex > 4,R_GF25 := Q_GF25 / Q_MFCex]
Result_var[Q_MFCex > 4,R_GF26 := Q_GF26 / Q_MFCex]

sub_Result_var <- Result_var[st > "2021-03-19 09:00" & st < "2021-03-20 09:00" & Conc_corr == "P23" 
				& Sonic %in% names(WDvarCols_10),]

Result_Qvar <- data.table(
 	melt(sub_Result_var, id.vars=c("st","Sonic"),measure.vars=c("R_GF17","R_GF18","R_GF16","R_GF25"), variable.name="Q_GF", value.name="Recovery")
	)

Result_Qvar$OP <- factor(Result_Qvar$Q_GF, levels=c("R_GF17","R_GF18","R_GF16","R_GF25"),labels=c("OP-50m","OP-100m","OP-150m","OP-200m"))

Fig_recovery_WDvar10 <- Result_Qvar[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & GF %in% c("GF-50m","GF-100m","GF-150m","GF-200m")
								 & Sonic %in% names(WDvarCols_10),{
	ggplot(.SD, aes(x=st,y=Recovery,colour=Sonic)) +
	geom_hline(yintercept=1,linetype=3) +
	geom_line(linewidth=0.6) +
	geom_point(size=1) +
	xlab(NULL) +
	scale_y_continuous(breaks=c(0,0.5,1,1.5),limits=c(0,1.75)) +
	ylab('IDM recovery rate [%]') +
	facet_grid(OP ~ .) +
	scale_colour_manual(name=NULL,values=WDvarCols_10) +
	theme_bw(base_size = 18) +
	theme(legend.position='none',strip.background =element_rect(fill="white"),panel.grid.major=element_line(size=0.6),panel.grid.minor=element_line(size=0.3))
}]

ggsave(file.path(PathFigures,"Initial submit/Initial submit - Fig. 8.png"),Fig_recovery_WDvar10,width = 30, height = 30/0.97, units="cm")

# the legned of the plot needs to be added manualy.


###########################################
###########################################
#####                                 #####
#####    Tables initial submisison    #####
#####                                 #####
###########################################
###########################################

############
### Table 1:

Result[L > 0, stability := 'stable']
Result[L < 0, stability := 'unstable']

# during entire MC
Result[!is.na(Sonic_ord) & Campaign == "MC", .(meanWS=round(mean(U_sonic,na.rm=TRUE),1),
	meanWD=round((360 + atan2(mean(U_sonic * sin(WD * pi/180),na.rm=TRUE),mean(U_sonic * cos(WD * pi/180),na.rm=TRUE)) * 180/pi) %% 360,1),
	meanUstar=round(mean(Ustar,na.rm=TRUE),2)),by=Sonic_ord]

# during release only
Result[!is.na(Sonic_ord) & Campaign == "MC" & Q_MFCex > 0, .(meanWS=round(mean(U_sonic,na.rm=TRUE),1),
	meanWD=round((360 + atan2(mean(U_sonic * sin(WD * pi/180),na.rm=TRUE),mean(U_sonic * cos(WD * pi/180),na.rm=TRUE)) * 180/pi) %% 360,1),
	meanUstar=round(mean(Ustar,na.rm=TRUE),2)),by=Sonic_ord]
### there is actually an error, as the threshold of Q_MFCex should be 4 and not 0.


############
### Table 2:

sonic_names <- Result[!is.na(Sonic),unique(Sonic)]

ls_16 <- ls_17 <- ls_18 <- ls_25 <- vector(mode="list",4)
for(j in 1:4){
	ls_16[[j]]	<- as.data.table(Result[Campaign == "MC" & Sonic == sonic_names[j] & Q_MFC > 4 & is.na(Q_GF16),.N]/
		Result[Campaign == "MC" & Sonic == sonic_names[j] & Q_MFC > 4 ,.N])
	ls_17[[j]]	<- as.data.table(Result[Campaign == "MC" & Sonic == sonic_names[j] & Q_MFC > 4 & is.na(Q_GF17),.N]/
		Result[Campaign == "MC" & Sonic == sonic_names[j] & Q_MFC > 4 ,.N])
	ls_18[[j]]	<- as.data.table(Result[Campaign == "MC" & Sonic == sonic_names[j] & Q_MFC > 4 & is.na(Q_GF18),.N]/
		Result[Campaign == "MC" & Sonic == sonic_names[j] & Q_MFC > 4 ,.N])
	ls_25[[j]]	<- as.data.table(Result[Campaign == "MC" & Sonic == sonic_names[j] & Q_MFC > 4 & is.na(Q_GF25),.N]/
		Result[Campaign == "MC" & Sonic == sonic_names[j] & Q_MFC > 4 ,.N])
}

data.loss <- round(rbind(rbindlist(ls_17), rbindlist(ls_18),rbindlist(ls_16),rbindlist(ls_25)),2)
colnames(data.loss) <- c("loss")
data.loss[,NAME := c(paste0(rep("GF17_",4),sonic_names),paste0(rep("GF18_",4),
	sonic_names),paste0(rep("GF16_",4),sonic_names),paste0(rep("GF25_",4),sonic_names))]

data.loss[grep('SonicA',NAME), ':=' (Sonic = 'UA-50m', GF = c(paste0('OP',c('50m','100m','150m','200m'))))]
data.loss[grep('SonicB',NAME), ':=' (Sonic = 'UA-150m', GF = c(paste0('OP',c('50m','100m','150m','200m'))))]
data.loss[grep('SonicC',NAME), ':=' (Sonic = 'UA-UW', GF = c(paste0('OP',c('50m','100m','150m','200m'))))]
data.loss[grep('Sonic2',NAME), ':=' (Sonic = 'UA-100m', GF = c(paste0('OP',c('50m','100m','150m','200m'))))]

data.loss

data.loss[Sonic=='UA-50m']
data.loss[Sonic=='UA-100m']
data.loss[Sonic=='UA-150m']
data.loss[Sonic=='UA-UW']

round(mean(data.loss[grep('SonicC',NAME),loss]),2)
round(mean(data.loss[grep('SonicA',NAME),loss]),2)
round(mean(data.loss[grep('Sonic2',NAME),loss]),2) # in the manuscript was 18%, which is probably wrong
round(mean(data.loss[grep('SonicB',NAME),loss]),2) # in the manuscript was 18%, which is probably wrong


############
### Table 3:

## recovery rates
Result_R <- data.table(rbind(
 	cbind(melt(Result[Sonic_ord == "UA-UW"], id.vars=c('st','L','U_sonic','Ustar','WD'),measure.vars=c("R_GF16","R_GF17","R_GF18","R_GF25","R_GF26"), variable.name="R_GF", value.name="Recovery"),Sonic="UA-UW")
 	,cbind(melt(Result[Sonic_ord == "UA-50m"], id.vars=c('st','L','U_sonic','Ustar','WD'),measure.vars=c("R_GF16","R_GF17","R_GF18","R_GF25","R_GF26"), variable.name="R_GF", value.name="Recovery"),Sonic="UA-50m")
 	,cbind(melt(Result[Sonic_ord == "UA-100m"], id.vars=c('st','L','U_sonic','Ustar','WD'),measure.vars=c("R_GF16","R_GF17","R_GF18","R_GF25","R_GF26"), variable.name="R_GF", value.name="Recovery"),Sonic="UA-100m")
 	,cbind(melt(Result[Sonic_ord == "UA-150m"], id.vars=c('st','L','U_sonic','Ustar','WD'),measure.vars=c("R_GF16","R_GF17","R_GF18","R_GF25","R_GF26"), variable.name="R_GF", value.name="Recovery"),Sonic="UA-150m")
	))
Result_R$OPold <- factor(Result_R$R_GF, levels=c("R_GF17","R_GF18","R_GF16","R_GF25"),labels=c("OP-50m","OP-100m","OP-150m","OP-200m"))

# revovery rate entire MC
Result_R[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & !is.na(OPold),
.(median=round(median(Recovery,na.rm=TRUE),2)),by=.(OPold,Sonic)][order(OPold,Sonic)]

Result_R[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & !is.na(OPold),
.(median=round(median(Recovery,na.rm=TRUE),2)),by=.(Sonic)][order(Sonic)]
Result_R[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & !is.na(OPold),
.(median=round(median(Recovery,na.rm=TRUE),2)),by=.(OPold)][order(OPold)]

