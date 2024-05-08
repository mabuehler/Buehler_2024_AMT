

#######################
#######################
#####             #####
#####    Plots    #####
#####             #####
#######################
#######################



#################################
### Header (Pfade, Libraries) ###
#################################

	library(ibts)
	library(RgoogleMaps)
	library(bLSmodelR)
	library(RColorBrewer)
	library(ggplot2)
	library(ggthemes)
	library(ggpointdensity)
	library(ggpubr)
	library(shape)
	library(plotly)
	library(viridis)



	PfadDaten <- "~/LFE/01_Projekte/14_Quellenversuch_Stall/Daten"		
	if(!dir.exists(PfadDaten)){
	}
	if(!dir.exists(PfadDaten)){
		# Marcel's Laptop
		PfadDaten <- "//bfh.ch/data/LFE/HAFL/MST-Hofduenger/01_Projekte/14_Quellenversuch_Stall/Daten"
	}
	PfadRSaves <- "~/repos/4_Projects/8_FerLAS/RSaves"
	PfadFigures <- file.path(dirname(PfadRSaves),"Figures")
	source("~/repos/2_Packages/1_bLSmodelR/contour_plots/contourXY.r")
	source("~/repos/2_Packages/1_bLSmodelR/contour_plots/contourXZ.r")
	source("~/repos/3_Scripts/gel-scripts/sonic-turbulence.r")
	source("~/repos/3_Scripts/gel-scripts/windrose.r")
	source("~/repos/3_Scripts/gel-scripts/wgs84-ch1903.R")

	Cat.Path <- paste0(PfadDaten,"/Catalogs")

###################################
### Funktionen und Definitionen ###
###################################

lines_sec2xy <- function(xyMK,sensor,node=1,wd,col="lightblue",lwd=2,...){
	# browser()
	GF <- xyMK[xyMK[,1] %in% sensor,]
	sens <- as.numeric(GF[GF[,3] == node,4:5])
	b <- tan((90 - wd)/180*pi)
	# x <- if(wd <= 180) 600 else -600
	x <- if(wd <= 180) 1E6 else -1E6
	y <- sens[2] - (sens[1] - x)*b
	lines(c(sens[1],x),c(sens[2],y),col=col,lwd=lwd,...)
}

	Col_CH4 <- "#1e88e5"
	# display.brewer.all()
	CH4Cols_GF <- CH4Cols_OP <- CH4Cols <- brewer.pal(5,"Dark2")
	names(CH4Cols) <- paste0("GF",c(16:18,25:26))
	names(CH4Cols_OP) <- paste0("OP-",c(1:5))
	names(CH4Cols_GF) <- paste0("GF-",c('UW','50m','100m','150m','200m'))

	DUACols <- SonicCols <- brewer.pal(4,"BrBG")
	names(SonicCols) <- paste0("Sonic",c("A","D","B","C"))
	names(DUACols) <- paste0("3DUA-",c('UW','50m','100m','150m'))

	# colSonics <- brewer.pal(4,"Set1")
	# colSonics <- c('#D00000','#FFBA08','#3F88C5','#032B43')
	colSonics <- c('#F8766D','#7CAE00','#00BFC4','#C77CFF')
	names(colSonics) <- c('UW','50m','100m','150m')

	WDmCols <- colorRampPalette(c("#00ff00", "#003300"))(15)
	WDpCols <- colorRampPalette(c("#00ffff", "#0033ff"))(15)
	# WDmCols <- colorRampPalette(c("#0000ff", "#00ffff"))(15)
	# WDpCols <- colorRampPalette(c("#ff0000", "#ffff00"))(15)
	WDvarCols <- c(WDmCols,"red",WDpCols)
	names(WDvarCols) <- paste0("SonicB",c(paste0("_m",1:15),"_0",paste0("_p",1:15)))
	
	mgs_to_kgd <- (24*3600)/1E6
	mgs_to_kgh <- (3600)/1E6
	mgm2s_to_gm2h <- 3600/1E3


	### MKs
	MKStart <- "19.01.2021"
	MKStop <- "26.03.2021 12:00"
	QV1 <- " to 28.01.2021 12:00"
	QV2 <- "05.03.2021 to 10.03.2021 18:00"
	QV1u2 <- " to 10.03.2021 18:00"
	MK <- "18.03.2021 11:00 - 21.03.2021 14:00"
	QV3	<- "21.03.2021 14:00 to "

#######################
### Laden der Daten ###
#######################

	# Pfadlängen & Geometrie & Map
	load(file.path(PfadRSaves,"Geometry_STO.RData"))
	STO_Map <- ReadMapTile(file.path(PfadFigures,"STO_GoogleMaps.png"))
	## bLS Resultate
	# bLS_Result <- readRDS(file = file.path(PfadRSaves, "STO_Result_orig_P6_10min.rds"))
	## Konzentrations Raw
	# GF_16 <- readRDS(file.path(PfadRSaves,"STO_GF_16.rds"))
	# GF_17 <- readRDS(file.path(PfadRSaves,"STO_GF_17.rds"))
	# GF_18 <- readRDS(file.path(PfadRSaves,"STO_GF_18.rds"))
	# GF_25 <- readRDS(file.path(PfadRSaves,"STO_GF_25.rds"))
	# GF_26 <- readRDS(file.path(PfadRSaves,"STO_GF_26.rds"))
	GFs10min <- readRDS(file.path(PfadRSaves,"STO_Conc_P23_10min.rds"))
	GFs10min_uncorr <- readRDS(file=file.path(PfadRSaves,"STO_Conc_uncorr_10min.rds"))
	GFs10min_corr1 <- readRDS(file=file.path(PfadRSaves,"STO_Conc_corr1_10min.rds"))

	## Sonic data
	load(file.path(PfadRSaves,"STO_Sonics_30min.RData"))
	load(file.path(PfadRSaves,"STO_Sonics_10min.RData"))
	load(file.path(PfadRSaves,"STO_Sonics_1min.RData"))
	load(file.path(PfadRSaves,"STO_Sonics_1sec.RData"))
	load(file=file.path(PfadRSaves,"STO_SonicC_planar.RData"))

	## weather station
	WS700_2 <- readRDS(file.path(PfadRSaves, "WS700_2_STO.rds"))
	WS2 <- shift_time(WS700_2[[1]]$data,"-5mins") # 5min Offset Korrektur (noch nicht angewendet)
	WS2[which(WS2$WS_0460 == 0),c("WS_0460","WD_corr")] <- NA_real_ 	# remove no wind data
	WS_10min <- pool(WS2,granularity="10mins",tz="Etc/GMT-1",st="2021-03-18 15:10")
	WS_1min <- pool(WS2[,c("WS_0480","WD_corr")],granularity="1mins",tz="Etc/GMT-1",st="2021-03-18 15:10")
	WS2_10min <- merge.ibts(pool(WS700_2[[2]]$data[,c("WS_0260", "WS_0360", "WS_0620", "WS_0625", "WS_0820")],granularity="10mins",st.to="07.12.2020 12:00"),WS700_2[[1]]$data[,c("WS_0160","WS_0480","WD_corr")])
	names(WS2_10min) <- c("relHum","Press","Precip1","Precip2","Precip3","Temp","WS","WD_WS")
	## contour plots
	XY_SonicA <- readRDS(file.path(PfadRSaves,"XY_SonicA_MK.rds"))
	XY_SonicB <- readRDS(file.path(PfadRSaves,"XY_SonicB_MK.rds"))
	XY_SonicC <- readRDS(file.path(PfadRSaves,"XY_SonicC_MK.rds"))
	XY_Sonic2 <- readRDS(file.path(PfadRSaves,"XY_Sonic2_MK.rds"))
	## WindTrax forwadrmodelling data
	# WT <- fread(file.path(dirname(PfadRSaves),"WindTrax_Forward/WindTrax_Output_Run1.txt"),sep=";")
	WT <- fread(file.path(dirname(PfadRSaves),"WindTrax_Forward/WindTrax_Output_Run2.txt"),sep=";")
	# WT_header <- fread(file.path(dirname(PfadRSaves),"WindTrax_Forward/WindTrax_Output_Run1.txt"),skip=7,sep=":",nrows=31,header=FALSE)
	# WT_header <- fread(file.path(dirname(PfadRSaves),"WindTrax_Forward/WindTrax_Output_Run1.txt"),skip=7,sep="",nrows=31,header=FALSE)
	WT_header <- fread(file.path(dirname(PfadRSaves),"WindTrax_Forward/WindTrax_Output_Run2.txt"),skip=7,sep="",nrows=31,header=FALSE)
	# colnames(WT) <- WT_header$V1
	## MFC und CC Sensor
	MFC <- readRDS(file.path(PfadRSaves,"STO_MFC.rds"))
	CC_Sensor <- readRDS(file.path(PfadRSaves,"STO_CC_Sensor_raw.rds"))
	## Emissionen
	Result <- readRDS(file=file.path(PfadRSaves,"STO_Results_orig_P23.rds"))
	Result_var <- readRDS(file=file.path(PfadRSaves,"STO_Result_variation.rds"))

 	centre <- colMeans(Sources[,c(2,3)])

##########################
### create XY Geometry ###
##########################

# Sensors_QV1u2_xy <- CH.to.xy(STO_Map,Sensors_QV1u2)
Sensors_MK_xy <- CH.to.map(STO_Map,Sensors_MK)
# Sensors_QV3_xy <- CH.to.map(STO_Map,Sensors_QV3)
Sources_xy <- CH.to.map(STO_Map,Sources)
# Sonics_QV1u2_xy <- CH.to.map(STO_Map,Sonics_QV1u2)
Sonics_MK_xy <- CH.to.map(STO_Map,Sonics_MK)
# Sonics_QV3_xy <- CH.to.map(STO_Map,Sonics_QV3)
# WS_xy <- CH.to.map(STO_Map,WeatherStation)
Baum_xy <- CH.to.map(STO_Map,Baum)

########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################


#########################################
#########################################
#####                               #####
#####    Scematic overview plots    #####
#####                               #####
#########################################
#########################################


# ### QV1u2
# png(file.path(PfadFigures,"Overview_QV1u2_schematic.png"),width=12,height=12,unit="in",res=300)
# plot(Sensors_QV1u2_xy, Sources_xy,sensors.text.args=list(labels=""),lines.args=list(lwd=2,col="black",lty=3)
# 	,points.args = list(pch = 20, cex = 2, col ="black"),polygon.args=list(col="grey95"),sources.text.args=list(labels=""))
# plot(Baum_xy,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 4, col ="grey60",bg="grey95",lwd=2),add=TRUE)
# plot(Sonics_QV1u2_xy,sensors.text.args=list(labels=""),points.args = list(pch = 17, cex = 2, col ="black"),add=TRUE)
# plot(WS_xy[WS_xy[,1] == "WS2",],sensors.text.args=list(labels=""),points.args = list(pch = "*", cex = 3, col ="black"),add=TRUE)
# dev.off()


### MK only
# pdf(file.path(PfadFigures,"Overview_MKuQV3_schematic.pdf"),width=16,height=8)
png(file.path(PfadFigures,"Overview_MK_schematic.png"),width=23,height=15,unit="cm",res=600)
par(mfrow=c(1,1),mar=c(0,0,0,0))
# MK
plot(type="n",1, xlim=c(centre[1]-255,centre[1]+200),ylim=c(centre[2]-215,centre[2]+97),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
# plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
# abline(v=seq(-400,length.out=10,by=100),h=seq(-400,length.out=10,by=100),lty=3,col="lightgrey",lwd=1.5)
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,code=1)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)

plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
plot(Sonics_MK,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
plot(WeatherStation[WeatherStation[,1] == "WS2",],sensors.text.args=list(labels=""),points.args = list(pch = 17, cex = 2, col ="black"),add=TRUE)
bx <- centre[1]+80
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
## Richtung Nord
wr.x <- bx+85
wr.y <- by+40  
polygon(c(wr.x,wr.x-5,wr.x),c(wr.y,wr.y-5,wr.y+40),border=1)
polygon(c(wr.x,wr.x+5,wr.x),c(wr.y,wr.y-5,wr.y+40),border=1,col=1)
text(x=wr.x,y=wr.y+50,labels="N",cex=1.5)
# text(x=centre[1]-185,y=centre[2]+120,"MC",cex=1.5)
## GasFinde anschreiben
text(x=Sensors_MK[1,][,4]+21,y=Sensors_MK[1,][,5]-4,"OP-150m")
text(x=Sensors_MK[3,][,4]+21,y=Sensors_MK[3,][,5]-4,"OP-50m")
text(x=Sensors_MK[5,][,4]+21,y=Sensors_MK[5,][,5]-4,"OP-100m")
text(x=Sensors_MK[7,][,4]+21,y=Sensors_MK[7,][,5]-4,"OP-200m")
text(x=Sensors_MK[10,][,4]-21,y=Sensors_MK[10,][,5],"OP-UW")
## Sonics anschreiben
text(x=Sonics_MK[1,][,4]-16,y=Sonics_MK[1,][,5]-15,"UA-100m")
text(x=Sonics_MK[2,][,4]-15,y=Sonics_MK[2,][,5]-15,"UA-50m")
text(x=Sonics_MK[3,][,4]-15,y=Sonics_MK[3,][,5]-15,"UA-150m")
text(x=Sonics_MK[4,][,4],y=Sonics_MK[4,][,5]-15,"UA-UW")

## Legende
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


### MK only but with hx names

# pdf(file.path(PfadFigures,"Overview_MKuQV3_schematic.pdf"),width=16,height=8)
png(file.path(PfadFigures,"Overview_MK_schematic_hx.png"),width=23,height=15,unit="cm",res=600)
par(mfrow=c(1,1),mar=c(0,0,0,0))
# MK
plot(type="n",1, xlim=c(centre[1]-255,centre[1]+200),ylim=c(centre[2]-215,centre[2]+97),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
# plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
# abline(v=seq(-400,length.out=10,by=100),h=seq(-400,length.out=10,by=100),lty=3,col="lightgrey",lwd=1.5)
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,code=1)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)

plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
plot(Sonics_MK,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
plot(WeatherStation[WeatherStation[,1] == "WS2",],sensors.text.args=list(labels=""),points.args = list(pch = 17, cex = 2, col ="black"),add=TRUE)
bx <- centre[1]+80
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
## Richtung Nord
wr.x <- bx+85
wr.y <- by+40  
polygon(c(wr.x,wr.x-5,wr.x),c(wr.y,wr.y-5,wr.y+40),border=1)
polygon(c(wr.x,wr.x+5,wr.x),c(wr.y,wr.y-5,wr.y+40),border=1,col=1)
text(x=wr.x,y=wr.y+50,labels="N",cex=1.5)
# text(x=centre[1]-185,y=centre[2]+120,"MC",cex=1.5)
## GasFinde anschreiben
text(x=Sensors_MK[1,][,4]+21,y=Sensors_MK[1,][,5]-4,"OP-8.6h")
text(x=Sensors_MK[3,][,4]+21,y=Sensors_MK[3,][,5]-4,"OP-2.0h")
text(x=Sensors_MK[5,][,4]+21,y=Sensors_MK[5,][,5]-4,"OP-5.3h")
text(x=Sensors_MK[7,][,4]+21,y=Sensors_MK[7,][,5]-4,"OP-12h")
text(x=Sensors_MK[10,][,4]-21,y=Sensors_MK[10,][,5],"OP-UW")
## Sonics anschreiben
text(x=Sonics_MK[1,][,4]-16,y=Sonics_MK[1,][,5]-15,"UA-5.3h")
text(x=Sonics_MK[2,][,4]-15,y=Sonics_MK[2,][,5]-15,"UA-2.0h")
text(x=Sonics_MK[3,][,4]-15,y=Sonics_MK[3,][,5]-15,"UA-8.6h")
text(x=Sonics_MK[4,][,4],y=Sonics_MK[4,][,5]-15,"UA-UW")

## Legende
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


### IC1, MC and IC2
{
png(file.path(PfadFigures,"Overview_IC1aMCaIC2_schematic.png"),width=30,height=10,unit="cm",res=600)
par(mfrow=c(1,3),mar=c(0.2,0,0.2,0))
# IC1
plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
Arrows(x0=Sensors_QV1u2[2,][,4],y0=Sensors_QV1u2[2,][,5],
	x1=Sensors_QV1u2[1,][,4],y1=Sensors_QV1u2[1,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_QV1u2[4,][,4],y0=Sensors_QV1u2[4,][,5],
	x1=Sensors_QV1u2[3,][,4],y1=Sensors_QV1u2[3,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_QV1u2[6,][,4],y0=Sensors_QV1u2[6,][,5],
	x1=Sensors_QV1u2[5,][,4],y1=Sensors_QV1u2[5,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_QV1u2[8,][,4],y0=Sensors_QV1u2[8,][,5],
	x1=Sensors_QV1u2[7,][,4],y1=Sensors_QV1u2[7,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_QV1u2[10,][,4],y0=Sensors_QV1u2[10,][,5],
	x1=Sensors_QV1u2[9,][,4],y1=Sensors_QV1u2[9,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_QV1u2[2,][,4],y0=Sensors_QV1u2[2,][,5],
	x1=Sensors_QV1u2[1,][,4],y1=Sensors_QV1u2[1,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_QV1u2[4,][,4],y0=Sensors_QV1u2[4,][,5],
	x1=Sensors_QV1u2[3,][,4],y1=Sensors_QV1u2[3,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_QV1u2[6,][,4],y0=Sensors_QV1u2[6,][,5],
	x1=Sensors_QV1u2[5,][,4],y1=Sensors_QV1u2[5,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_QV1u2[8,][,4],y0=Sensors_QV1u2[8,][,5],
	x1=Sensors_QV1u2[7,][,4],y1=Sensors_QV1u2[7,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_QV1u2[10,][,4],y0=Sensors_QV1u2[10,][,5],
	x1=Sensors_QV1u2[9,][,4],y1=Sensors_QV1u2[9,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
plot(Sonics_QV1u2,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black", lwd=3),add=TRUE)
plot(WeatherStation[WeatherStation[,1] == "WS2",],sensors.text.args=list(labels=""),points.args = list(pch = 17, cex = 2, col ="black"),add=TRUE)
text(x=centre[1]-185,y=centre[2]+120,"IC1",cex=1.5)

# MC
plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
# abline(v=seq(-400,length.out=10,by=100),h=seq(-400,length.out=10,by=100),lty=3,col="lightgrey",lwd=1.5)
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,code=1)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)

plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
plot(Sonics_MK,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
plot(WeatherStation[WeatherStation[,1] == "WS2",],sensors.text.args=list(labels=""),points.args = list(pch = 17, cex = 2, col ="black"),add=TRUE)
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
## Richtung Nord
wr.x <- bx+85
wr.y <- by+40  
polygon(c(wr.x,wr.x-5,wr.x),c(wr.y,wr.y-5,wr.y+40),border=1)
polygon(c(wr.x,wr.x+5,wr.x),c(wr.y,wr.y-5,wr.y+40),border=1,col=1)
text(x=wr.x,y=wr.y+50,labels="N",cex=1.5)
text(x=centre[1]-185,y=centre[2]+120,"MC",cex=1.5)
# ## GasFinde anschreiben
# text(x=Sensors_MK[1,][,4]+20,y=Sensors_MK[1,][,5]-4,"GF16")
# text(x=Sensors_MK[3,][,4]+20,y=Sensors_MK[3,][,5]-4,"GF17")
# text(x=Sensors_MK[5,][,4]+20,y=Sensors_MK[5,][,5]-4,"GF18")
# text(x=Sensors_MK[7,][,4]+20,y=Sensors_MK[7,][,5]-4,"GF25")
# text(x=Sensors_MK[10,][,4]-21,y=Sensors_MK[10,][,5],"GF26")
# ## Sonics anschreiben
# text(x=Sonics_MK[1,][,4]-13,y=Sonics_MK[1,][,5]-15,"SonicD")
# text(x=Sonics_MK[2,][,4]-13,y=Sonics_MK[2,][,5]-15,"SonicA")
# text(x=Sonics_MK[3,][,4]-13,y=Sonics_MK[3,][,5]-15,"SonicB")
# text(x=Sonics_MK[4,][,4],y=Sonics_MK[4,][,5]-15,"SonicC")

# IC2
plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
Arrows(x0=Sensors_QV3[1,][,4],y0=Sensors_QV3[1,][,5],
	x1=Sensors_QV3[2,][,4],y1=Sensors_QV3[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_QV3[3,][,4],y0=Sensors_QV3[3,][,5],
	x1=Sensors_QV3[4,][,4],y1=Sensors_QV3[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_QV3[5,][,4],y0=Sensors_QV3[5,][,5],
	x1=Sensors_QV3[6,][,4],y1=Sensors_QV3[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_QV3[7,][,4],y0=Sensors_QV3[7,][,5],
	x1=Sensors_QV3[8,][,4],y1=Sensors_QV3[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_QV3[9,][,4],y0=Sensors_QV3[9,][,5],
	x1=Sensors_QV3[10,][,4],y1=Sensors_QV3[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_QV3[1,][,4],y0=Sensors_QV3[1,][,5],
	x1=Sensors_QV3[2,][,4],y1=Sensors_QV3[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_QV3[3,][,4],y0=Sensors_QV3[3,][,5],
	x1=Sensors_QV3[4,][,4],y1=Sensors_QV3[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_QV3[5,][,4],y0=Sensors_QV3[5,][,5],
	x1=Sensors_QV3[6,][,4],y1=Sensors_QV3[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_QV3[7,][,4],y0=Sensors_QV3[7,][,5],
	x1=Sensors_QV3[8,][,4],y1=Sensors_QV3[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_QV3[9,][,4],y0=Sensors_QV3[9,][,5],
	x1=Sensors_QV3[10,][,4],y1=Sensors_QV3[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
plot(Sonics_QV3[Sonics_QV3[,1] == 'Sonic2',] ,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black", lwd=3),add=TRUE)
plot(WeatherStation[WeatherStation[,1] == "WS2",],sensors.text.args=list(labels=""),points.args = list(pch = 17, cex = 2, col ="black"),add=TRUE)
text(x=centre[1]-185,y=centre[2]+120,"IC2",cex=1.5)

## Legende
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


####################################################
####################################################
#####                                          #####
#####    Scematic Overview & Windrose plots    #####
#####                                          #####
####################################################
####################################################

MFC_10min <- pool(MFC,granularity="10mins", st='01.03.2021 10:00')
MFC_1min <- pool(MFC,granularity="1mins")
which(MFC_10min[MK,"Q_MFC"] > 0)
which(MFC_1min[MK,"Q_MFC"] > 0)
# start_1 "18.03.2021 12:57" - 5min
# end_1 "18.03.2021 13:50" + 10min
# indMK_MFC_1 <- parse_date_time3("18.03.2021 12:50",tz="Etc/GMT-1") + c(0:7)*600 
# start_2 "19.03.2021 10:30" - 5min
# end_2 "19.03.2021 16:48" + 10min
indMK_MFC_1 <- parse_date_time3("19.03.2021 10:30",tz="Etc/GMT-1") + c(0:37)*600
# start_3 "19.03.2021 21:52" - 5min
# end_3 "20.03.2021 06:52" + 10min
indMK_MFC_2 <- parse_date_time3("19.03.2021 21:50",tz="Etc/GMT-1") + c(0:54)*600
indMK_MFC_3 <- c(indMK_MFC_1,indMK_MFC_2)

length(st(SonicA[indMK_MFC_3]))
length(st(SonicB[indMK_MFC_3]))
length(st(SonicC[indMK_MFC_3]))
length(st(Sonic2[indMK_MFC_3]))
# Sonic2 hat am wenigsten daten
indMK_MFC <- st(Sonic2[indMK_MFC_3])

max(SonicA[indMK_MFC,"Ustar"],na.rm=TRUE)
max(SonicB[indMK_MFC,"Ustar"],na.rm=TRUE)
max(SonicC[indMK_MFC,"Ustar"],na.rm=TRUE)
max(Sonic2[indMK_MFC,"Ustar"],na.rm=TRUE)

####################3

bx <- centre[1]+30
by <- centre[2]-205
dx1 <- 50
dx2 <- 100
dy <- 4
wr.x <- bx+85
wr.y <- by+40 

############################
### u* friciton velocity ###
############################

png(file.path(PfadFigures,"Windrose_SonicAll_corrected.png"),width=30,height=30,unit="cm",res=600)
par(mfrow=c(2,2),mar=c(0,0,0,0))
## SonicC
plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,code=1)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
plot(Sonics_MK,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
# windrose(SonicC["2021-03-19 09:00 - 2021-03-20 09:00"], wd = "WD", ws = "Ustar", max_freq = 16, draw.grid = FALSE,
windrose(SonicC[indMK_MFC_3], wd = "WD", ws = "Ustar", max_freq = 16, draw.grid = FALSE,
    ws_breaks = c(0,0.15,0.25,0.35,0.46),delta_wd = 12, add = TRUE, scale = 0.6,legend.x=NA,
    # ws_breaks = c(0,0.15,0.25,0.32,0.46),delta_wd = 12, add = TRUE, scale = 0.6,
    # center = c(-150,0),unit=expression("u"["*"]*" [m s"^-1*"]"),legend.x = "bottomright",legend.cex=2,legend.bty="n")
    center = c(mean(Sensors_MK[Sensors_MK[1]=="GF25",4]),mean(Sensors_MK[Sensors_MK[1]=="GF25",5])))
text(x=centre[1]-200,y=centre[2]+115,"3DUA-UW",cex=2,pos=4)

## SonicA
plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,code=1)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
plot(Sonics_MK,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
# windrose(SonicA["2021-03-19 09:00 - 2021-03-20 09:00"], wd = "WD", ws = "Ustar", max_freq = 16, draw.grid = FALSE,
windrose(SonicA[indMK_MFC_3], wd = "WD", ws = "Ustar", max_freq = 16, draw.grid = FALSE,
    ws_breaks = c(0,0.15,0.25,0.35,0.46),delta_wd = 12, add = TRUE, scale = 0.6,
    center = c(Sonics_MK[Sonics_MK[1] == "SonicA",4],Sonics_MK[Sonics_MK[1] == "SonicA",5]),legend.x = "NA")
text(x=centre[1]-200,y=centre[2]+115,"3DUA-50m",cex=2,pos=4)
## Legende
dyl <- 16
dxl <- 8
xl <- centre[1] + 62
yl <- centre[2] - 152 - dyl
Arrows(x0=xl-dxl,x1=xl+dxl,y0=yl,y1=yl,arr.type="triangle",arr.length=0.15,lwd=1.5,arr.adj=1)
Arrows(x0=xl+dxl,x1=xl-dxl,y0=yl,y1=yl,arr.type="T",arr.length=0.15,lwd=1.5,arr.adj=1)
points(x=xl,y=yl-dyl,pch=8,cex=1.5)
points(x=xl,y=yl-dyl*2,pch=21,cex=2.3,col='grey60',bg='grey95')
rect(xleft=xl-dxl,xright=xl+dxl,ybottom=yl-dyl*3-4,ytop=yl-dyl*3+4,col="grey95")
text(x=rep(xl+dxl*2,4),y=c(yl,yl-dyl,yl-dyl*2,yl-dyl*3),labels=c("GasFinder","3DUA","Tree","Barn"),pos=4,cex=1.5)
lines(x=c(xl-dxl*2,xl-dxl*2,xl+200),y=c(yl-100,yl+dyl,yl+dyl))

## Sonic2 (SonicD)
plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,code=1)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
plot(Sonics_MK,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
# windrose(Sonic2["2021-03-19 09:00 - 2021-03-20 09:00"], wd = "WD", ws = "Ustar", max_freq = 16, draw.grid = FALSE,
windrose(Sonic2[indMK_MFC_3], wd = "WD", ws = "Ustar", max_freq = 16, draw.grid = FALSE,
    ws_breaks = c(0,0.15,0.25,0.35,0.46),delta_wd = 12, add = TRUE, scale = 0.6,
    center = c(Sonics_MK[Sonics_MK[1] == "Sonic2",4],Sonics_MK[Sonics_MK[1] == "Sonic2",5]),unit=expression("u"["*"]*" [m s"^-1*"]"), 
    legend.x= NA)
    # legend.x= "bottomright",legend.cex=2,legend.bty="n")
text(x=centre[1]-200,y=centre[2]+115,"3DUA-100m",cex=2,pos=4)
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
plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,code=1)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
plot(Sonics_MK,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
# windrose(SonicB["2021-03-19 09:00 - 2021-03-20 09:00"], wd = "WD", ws = "Ustar", max_freq = 16, draw.grid = FALSE,
windrose(SonicB[indMK_MFC_3], wd = "WD", ws = "Ustar", max_freq = 16, draw.grid = FALSE,
    ws_breaks = c(0,0.15,0.25,0.35,0.46),delta_wd = 12, add = TRUE, scale = 0.6,
    center = c(Sonics_MK[Sonics_MK[1] == "SonicB",4],Sonics_MK[Sonics_MK[1] == "SonicB",5]),unit=expression("u"["*"]*" [m s"^-1*"]"),legend.x = "bottomright",legend.cex=1.5,legend.bty="n")
lines(x=c(bx+28,bx+28,bx+200),y=c(by-100,by+75,by+75))
text(x=centre[1]-200,y=centre[2]+115,"3DUA-150m",cex=2,pos=4)

dev.off()

##################################################################
### friction velocity but removed paths and UW placed on the side:
{
png(file.path(PfadFigures,"Windrose_SonicAll_ustar_paper.png"),width=30,height=30,unit="cm",res=600)
par(mfrow=c(2,2),mar=c(0,0,0,0))
## SonicC
# plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(type="n",1, xlim=c(centre[1]-210,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
windrose(SonicC[indMK_MFC_3], wd = "WD", ws = "Ustar", max_freq = 16, draw.grid = FALSE,
    ws_breaks = c(0,0.15,0.25,0.35,0.46),delta_wd = 12, add = TRUE, scale = 0.6,legend.x=NA, start=6, center = c(583707,210173))
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,code=1,col='grey60')
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
plot(Sonics_MK[Sonics_MK[1]=='SonicC',],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
text(x=centre[1]-225,y=centre[2]+130,"UA-UW",cex=2,pos=4)

## SonicA
plot(type="n",1, xlim=c(centre[1]-210,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
windrose(SonicA[indMK_MFC_3], wd = "WD", ws = "Ustar", max_freq = 16, draw.grid = FALSE,
    ws_breaks = c(0,0.15,0.25,0.35,0.46),delta_wd = 12, add = TRUE, scale = 0.6,legend.x= NA, start=6, center = c(583707,210173))
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60',code=1)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
plot(Sonics_MK[Sonics_MK[1] == 'SonicA',],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
text(x=centre[1]-225,y=centre[2]+130,"UA-2.0h",cex=2,pos=4)
## Legende
dyl <- 16
dxl <- 8
xl <- centre[1] + 62
yl <- centre[2] - 152 - dyl
Arrows(x0=xl-dxl,x1=xl+dxl,y0=yl,y1=yl,arr.type="triangle",arr.length=0.15,lwd=1.5,arr.adj=1)
Arrows(x0=xl+dxl,x1=xl-dxl,y0=yl,y1=yl,arr.type="T",arr.length=0.15,lwd=1.5,arr.adj=1)
points(x=xl,y=yl-dyl,pch=8,cex=1.5)
points(x=xl,y=yl-dyl*2,pch=21,cex=2.3,col='grey60',bg='grey95')
rect(xleft=xl-dxl,xright=xl+dxl,ybottom=yl-dyl*3-4,ytop=yl-dyl*3+4,col="grey95")
text(x=rep(xl+dxl*2,4),y=c(yl,yl-dyl,yl-dyl*2,yl-dyl*3),labels=c("OP","UA","Tree","Barn"),pos=4,cex=1.5)
lines(x=c(xl-dxl*2,xl-dxl*2,xl+200),y=c(yl-100,yl+dyl,yl+dyl))

## Sonic2 (SonicD)
plot(type="n",1, xlim=c(centre[1]-210,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
windrose(Sonic2[indMK_MFC_3], wd = "WD", ws = "Ustar", max_freq = 16, draw.grid = FALSE,
    ws_breaks = c(0,0.15,0.25,0.35,0.46),delta_wd = 12, add = TRUE, scale = 0.6,
    unit=expression("u"["*"]*" [m s"^-1*"]"), legend.x= NA, start=6, center = c(583707,210173))
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60',code=1)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
plot(Sonics_MK[Sonics_MK[1]=='Sonic2',],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
text(x=centre[1]-225,y=centre[2]+130,"UA-5.3h",cex=2,pos=4)
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
plot(type="n",1, xlim=c(centre[1]-210,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
windrose(SonicB[indMK_MFC_3], wd = "WD", ws = "Ustar", max_freq = 16, draw.grid = FALSE,
    ws_breaks = c(0,0.15,0.25,0.35,0.46),delta_wd = 12, add = TRUE, scale = 0.6, unit=expression("u"["*"]*" [m s"^-1*"]"),
    legend.x = "bottomright",legend.cex=1.5,legend.bty="n", start=6, center=c(583707,210173))
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60',code=1)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2,col='grey60')
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2,col='grey60')
plot(Sonics_MK[Sonics_MK[1]=='SonicB',],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
lines(x=c(bx+28,bx+28,bx+200),y=c(by-100,by+75,by+75))
text(x=centre[1]-225,y=centre[2]+130,"UA-8.6h",cex=2,pos=4)

dev.off()
}

##################
### Wind speed ###
##################

max(SonicA[indMK_MFC,"U_sonic"],na.rm=TRUE)
max(SonicB[indMK_MFC,"U_sonic"],na.rm=TRUE)
max(SonicC[indMK_MFC,"U_sonic"],na.rm=TRUE)
max(Sonic2[indMK_MFC,"U_sonic"],na.rm=TRUE)


png(file.path(PfadFigures,"Windrose_U_sonic.png"),width=30,height=30,unit="cm",res=600)
par(mfrow=c(2,2),mar=c(0,0,0,0))
## SonicC
plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,code=1)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
plot(Sonics_MK,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
windrose(SonicC[indMK_MFC_3], wd = "WD", ws = "U_sonic", max_freq = 16, draw.grid = FALSE,
    ws_breaks = c(0,1.5,2.5,3.5,6.5),delta_wd = 12, add = TRUE, scale = 0.6,legend.x=NA,
    center = c(mean(Sensors_MK[Sensors_MK[1]=="GF25",4]),mean(Sensors_MK[Sensors_MK[1]=="GF25",5])))
text(x=centre[1]-200,y=centre[2]+115,"3DUA-UW",cex=2,pos=4)

## SonicA
plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,code=1)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
plot(Sonics_MK,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
windrose(SonicA[indMK_MFC_3], wd = "WD", ws = "U_sonic", max_freq = 16, draw.grid = FALSE,
    ws_breaks = c(0,1.5,2.5,3.5,6.5),delta_wd = 12, add = TRUE, scale = 0.6,
    center = c(Sonics_MK[Sonics_MK[1] == "SonicA",4],Sonics_MK[Sonics_MK[1] == "SonicA",5]),legend.x = "NA")
text(x=centre[1]-200,y=centre[2]+115,"3DUA-50m",cex=2,pos=4)
## Legende
dyl <- 16
dxl <- 8
xl <- centre[1] + 62
yl <- centre[2] - 152 - dyl
Arrows(x0=xl-dxl,x1=xl+dxl,y0=yl,y1=yl,arr.type="triangle",arr.length=0.15,lwd=1.5,arr.adj=1)
Arrows(x0=xl+dxl,x1=xl-dxl,y0=yl,y1=yl,arr.type="T",arr.length=0.15,lwd=1.5,arr.adj=1)
points(x=xl,y=yl-dyl,pch=8,cex=1.5)
points(x=xl,y=yl-dyl*2,pch=21,cex=2.3,col='grey60',bg='grey95')
rect(xleft=xl-dxl,xright=xl+dxl,ybottom=yl-dyl*3-4,ytop=yl-dyl*3+4,col="grey95")
text(x=rep(xl+dxl*2,4),y=c(yl,yl-dyl,yl-dyl*2,yl-dyl*3),labels=c("OP","3DUA","Tree","Barn"),pos=4,cex=1.5)
lines(x=c(xl-dxl*2,xl-dxl*2,xl+200),y=c(yl-100,yl+dyl,yl+dyl))

## Sonic2 (SonicD)
plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,code=1)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
plot(Sonics_MK,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
windrose(Sonic2[indMK_MFC_3], wd = "WD", ws = "U_sonic", max_freq = 16, draw.grid = FALSE,
    ws_breaks = c(0,1.5,2.5,3.5,6.5),delta_wd = 12, add = TRUE, scale = 0.6,
    center = c(Sonics_MK[Sonics_MK[1] == "Sonic2",4],Sonics_MK[Sonics_MK[1] == "Sonic2",5]),unit=expression("u"["*"]*" [m s"^-1*"]"), 
    legend.x= NA)
text(x=centre[1]-200,y=centre[2]+115,"3DUA-100m",cex=2,pos=4)
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
plot(type="n",1, xlim=c(centre[1]-190,centre[1]+145),ylim=c(centre[2]-165,centre[2]+70),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2,code=1)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[9,][,4],y0=Sensors_MK[9,][,5],
	x1=Sensors_MK[10,][,4],y1=Sensors_MK[10,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
plot(Sonics_MK,sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
windrose(SonicB[indMK_MFC_3], wd = "WD", ws = "U_sonic", max_freq = 16, draw.grid = FALSE,
    ws_breaks = c(0,1.5,2.5,3.5,6.5),delta_wd = 12, add = TRUE, scale = 0.6,
    center = c(Sonics_MK[Sonics_MK[1] == "SonicB",4],Sonics_MK[Sonics_MK[1] == "SonicB",5]),unit=expression("u  [m s"^-1*"]"),legend.x = "bottomright",legend.cex=1.5,legend.bty="n")
lines(x=c(bx+28,bx+28,bx+200),y=c(by-100,by+75,by+75))
text(x=centre[1]-200,y=centre[2]+115,"3DUA-150m",cex=2,pos=4)

dev.off()



###############################
###############################
#####                     #####
#####    Contour plots    #####
#####                     #####
###############################
###############################

########################
### xy contour plots ###
########################

## SonicA

# NE
# png(filename = file.path(PfadFigures,"contour_XY_SonicA.png"),width=900, height=900, units="px")
# contourXY(XY_SonicA, fill = TRUE, levels = function(x)c(0.05, 0.1, 0.15, 0.2, 0.3, 0.5))
# contourXY(XY_SonicA, fill = TRUE, stats="CE")
# contourXY(XY_SonicA, fill = TRUE, levels = function(x) quantile(range(x, na.rm = TRUE), c(0.01, 0.05, 0.1, 0.5, 0.9)))
# contourXY(XY_SonicA, fill = TRUE, levels = function(x) (x*6*1E6*3600, c(0.03, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5)))

XY_SonicA_Conc <- XY_SonicA
XY_SonicA_Conc[[1]]$z <- XY_SonicA_Conc[[1]]$z * 6/mgs_to_kgh/681.43 # Fläche der Quelle (681.43 m2)
XY_SonicB_Conc <- XY_SonicB
XY_SonicB_Conc[[1]]$z <- XY_SonicB_Conc[[1]]$z * 6/mgs_to_kgh/681.43
XY_SonicC_Conc <- XY_SonicC
XY_SonicC_Conc[[1]]$z <- XY_SonicC_Conc[[1]]$z * 6/mgs_to_kgh/681.43
XY_Sonic2_Conc <- XY_Sonic2
XY_Sonic2_Conc[[1]]$z <- XY_Sonic2_Conc[[1]]$z * 6/mgs_to_kgh/681.43

### Schematische Darstellung
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


png(file.path(PfadFigures,"ContourXY_All_corrected.png"),width=30,height=30,unit="cm",res=600)
par(mfrow=c(2,2),mar=c(0,3.7,3.7,0))

## SonicC
plot(type="n",1, xlim=c(centre[1]-220,centre[1]+60),ylim=c(centre[2]-220,centre[2]+20),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
plot(Sonics_MK[Sonics_MK[,1] == "SonicC",],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
# contourXY(XY_SonicC_Conc, fill = TRUE,add=TRUE,showLegend=TRUE,lpos="bottomright")
contourXY(XY_SonicC_Conc, fill = TRUE,add=TRUE,showLegend=FALSE)
# text(x=centre[1]-180,y=centre[2]+35,"3DUA-UW",cex=2)
text(x=centre[1]-190,y=centre[2]+35,"UA-UW",cex=2)

## SonicA
par(mar=c(0,0,3.7,3.7))
plot(type="n",1, xlim=c(centre[1]-220,centre[1]+60),ylim=c(centre[2]-220,centre[2]+20),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
plot(Sonics_MK[Sonics_MK[,1] == "SonicA",],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
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
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
plot(Sonics_MK[Sonics_MK[,1] == "Sonic2",],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
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
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
plot(Sonics_MK[Sonics_MK[,1] == "SonicB",],sensors.text.args=list(labels=""),points.args = list(pch = 8, cex = 1.5, col ="black",lwd=3),add=TRUE)
contourXY(XY_SonicB_Conc, fill = TRUE,add=TRUE,showLegend=FALSE)
# text(x=centre[1]-180,y=centre[2]+35,"3DUA-150m",cex=2)
text(x=centre[1]-190,y=centre[2]+35,"UA-150m",cex=2)
# legend
legend("bottomright",legend=c("(0.06,0.3]","(0.3,0.6]","(0.6,3]","(3,5]",">5"),fill=paste0(ConcPalette(5),"4c"),title=expression("Conc. [mg m"^-3*"]"),cex=1.5,bty="n")
lines(x=c(bx+35,bx+35,bx+200),y=c(by-100,by+85,by+85))

dev.off()
######################


### second plume
## SonicC
XY_SonicC_Conc_shift <- XY_SonicC_Conc
XY_SonicC_Conc_shift[[1]]$x <- XY_SonicC_Conc[[1]]$x-65
XY_SonicC_Conc_shift[[1]]$y <- XY_SonicC_Conc[[1]]$y-5

# png(file.path(PfadFigures,"ContourXY_2plume_corrected.png"),width=30,height=30,unit="cm",res=600)
pdf(file.path(PfadFigures,"ContourXY_2plume_corrected.pdf"),width=8,height=8)

par(mar=c(0,0,0,0)+3.5)
plot(type="n",1, xlim=c(centre[1]-220,centre[1]+60),ylim=c(centre[2]-220,centre[2]+20),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Sources,polygon.args=list(col="grey95"),sources.text.args=list(labels=""),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="triangle",lwd=2,lty=1,arr.adj=1,arr.length=0.2)
Arrows(x0=Sensors_MK[1,][,4],y0=Sensors_MK[1,][,5],
	x1=Sensors_MK[2,][,4],y1=Sensors_MK[2,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[3,][,4],y0=Sensors_MK[3,][,5],
	x1=Sensors_MK[4,][,4],y1=Sensors_MK[4,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[5,][,4],y0=Sensors_MK[5,][,5],
	x1=Sensors_MK[6,][,4],y1=Sensors_MK[6,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
Arrows(x0=Sensors_MK[7,][,4],y0=Sensors_MK[7,][,5],
	x1=Sensors_MK[8,][,4],y1=Sensors_MK[8,][,5],arr.type="T",lwd=2,lty=1,arr.adj=1,code=1,arr.length=0.2)
contourXY(XY_SonicC, fill = TRUE,add=TRUE,showLegend=FALSE,col=brewer.pal(5,"Reds"))
# contourXY(XY_SonicC_Conc_shift, fill = TRUE,add=TRUE,col=brewer.pal(5,"Blues"),showLegend=FALSE)
text(x=centre[1]-198,y=centre[2]+35,"3DUA-1",cex=2)
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
dev.off()
######################

pdf(file.path(PfadFigures,"ContourXY_2plumeblue_corrected.pdf"),width=8,height=8)

par(mar=c(0,0,0,0)+3.5)
plot(type="n",1, xlim=c(centre[1]-220,centre[1]+60),ylim=c(centre[2]-220,centre[2]+20),xaxt="n",yaxt="n", xlab="",ylab="",asp=1)

contourXY(XY_SonicC_Conc, fill = TRUE,add=TRUE,col=brewer.pal(5,"Blues"),showLegend=FALSE)

dev.off()


######################################
######################################
#####                            #####
#####    Weatherstation plots    #####
#####                            #####
######################################
######################################

WS2_10min['19.03.2021 11:10 - 19.03.2021 11:30', c('relHum','Press','Temp','WS','WD_WS')] <- NA_real_ # set the values during power outage to NA
WS2_dt <- cbind(as.data.table(WS2_10min),st=st(WS2_10min))
# shading_release <- data.table(x1min = parse_date_time3("19.03.2021 10:30",tz="Etc/GMT-1"),x1max=parse_date_time3("19.03.2021 16:30",tz="Etc/GMT-1")
# 				,x2min = parse_date_time3("19.03.2021 21:30",tz="Etc/GMT-1"),x2max=parse_date_time3("20.03.2021 08:30",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf)

limits_WS  <- parse_date_time3(c("18.03.2021 16:00","22.03.2021 18:00"),tz="Etc/GMT-1")
limits_WSMK  <- parse_date_time3(c("18.03.2021 11:00","21.03.2021 14:00"),tz="Etc/GMT-1")
mabreaks <- parse_date_time3("18.03.2021 00:00",tz="Etc/GMT-1") + seq(0,5)*86400
mibreaks <- parse_date_time3("18.03.2021 12:00",tz="Etc/GMT-1") + seq(0,4)*86400
fill.col <- "grey70"

WS_Temp <- WS2_dt[st > "2021-03-18",{
 	ggplot(.SD,aes(x=st,y=Temp)) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill=fill.col,alpha=0.01) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill=fill.col,alpha=0.01) +
 	geom_rect(aes(xmin=parse_date_time3("22.03.2021 09:30",tz="Etc/GMT-1"),xmax=parse_date_time3("22.03.2021 14:30",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill=fill.col,alpha=0.01) +
 	geom_hline(yintercept = 0, colour="grey") +
 	geom_line() +
 	# geom_line(col="orange",lwd=1.2) +
 	scale_x_datetime(limits=limits_WS,expand=c(0.02,0.02)) +
 	ylim(-5,10) +
 	ylab("Temperature [°C]") +
 	xlab(NULL) +
 	theme_classic(base_size=18)
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_blank(),plot.margin=unit(c(0,0,0,0),units="lines"),legend.position = "none", text = element_text(size=30))
 }]


WS_TempMK <- WS2_dt[st >= parse_date_time3("18.03.2021 11:00",tz='Etc/GMT-1') & st <= parse_date_time3('21.03.2021 14:00',tz='Etc/GMT-1'),{
 	ggplot(.SD,aes(x=st,y=Temp)) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_hline(yintercept = 0, colour="grey") +
 	geom_line() +
 	# geom_line(col="orange",lwd=1.2) +
 	scale_x_datetime(limits=limits_WSMK,expand=c(0.02,0.02)) +
 	ylim(-5,10) +
 	ylab("Temperature [°C]") +
 	xlab(NULL) +
 	theme_bw(base_size=18) +
 	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3))
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_blank(),plot.margin=unit(c(0,0,0,0),units="lines"),legend.position = "none", text = element_text(size=30))
 }]

# WS_Press <- WS2_dt[st > "2021-03-18",{
#  	ggplot(.SD,aes(x=st,y=Press)) +
#  	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill=fill.col,alpha=0.01) +
#  	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill=fill.col,alpha=0.01) +
#  	geom_rect(aes(xmin=parse_date_time3("22.03.2021 09:30",tz="Etc/GMT-1"),xmax=parse_date_time3("22.03.2021 14:30",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill=fill.col,alpha=0.01) +
#  	geom_line() +
#  	# geom_line(col="purple",lwd=1.2) +
#  	scale_x_datetime(limits=limits_WS,expand=c(0.02,0.02)) +
#  	scale_y_continuous(breaks=seq(963,975,by=3),limits=c(962,975)) +
#  	ylab("Pressure [hPa]") +
#  	xlab(NULL) +
#  	theme_classic(base_size=18)
#   #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_blank(),plot.margin=unit(c(0,0,0,0),units="lines"),legend.position = "none", text = element_text(size=30))
#  }]

WS_WD <- WS2_dt[st > "2021-03-18",{
 	# ggplot(.SD,aes(x=st,y=(WD_WS + 180) %% 360 -180)) +
 	ggplot(.SD,aes(x=st,y=WD_WS)) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill=fill.col,alpha=0.01) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill=fill.col,alpha=0.01) +
 	geom_rect(aes(xmin=parse_date_time3("22.03.2021 09:30",tz="Etc/GMT-1"),xmax=parse_date_time3("22.03.2021 14:30",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill=fill.col,alpha=0.01) +
 	geom_line() +
 	scale_x_datetime(limits=limits_WS,expand=c(0.02,0.02)) +
	scale_y_continuous(breaks = seq(0, 360, by=90), limits=c(0,360)) +
 	ylab("Wind direction [°]") +
 	xlab(NULL) +
 	theme_classic(base_size=18)
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

WS_WS <- WS2_dt[st > "2021-03-18",{
 	ggplot(.SD,aes(x=st,y=WS)) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill=fill.col,alpha=0.01) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill=fill.col,alpha=0.01) +
 	geom_rect(aes(xmin=parse_date_time3("22.03.2021 09:30",tz="Etc/GMT-1"),xmax=parse_date_time3("22.03.2021 14:30",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill=fill.col,alpha=0.01) +
 	geom_line() +
 	# geom_line(col="darkgreen",lwd=1.2) +
 	scale_x_datetime(limits=limits_WS,expand=c(0.02,0.02),breaks="1 day",minor_breaks="6 hours",date_labels = "%b %d") +
 	ylab(expression("Wind speed [m s"^-1*"]")) +
 	xlab(NULL) +
 	theme_classic(base_size=18)
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.margin=unit(c(0,0,0,0),units="lines"),legend.position = "none"
  #   	, text = element_text(size=30))
    	# , text = element_text(size=30),axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
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


Meteo <- ggarrange(WS_Temp,WS_WD,WS_WS,ncol=1,align="hv")
Meteo_MK <- ggarrange(WS_TempMK,WS_WDMK,WS_WSMK,ncol=1,align="hv")
ggsave(file.path(PfadFigures, "Meteo_STO.png"), Meteo, width = 30, height = 30/1.57, units="cm")
ggsave(file.path(PfadFigures, "Meteo_STO_MK.png"), Meteo_MK, width = 30, height = 30/1.57, units="cm")


################################################
### Combination of Weather station and UA-UW ###
################################################

limits_WSUA  <- parse_date_time3(c("19.03.2021","21.03.2021"),tz="Etc/GMT-1")
mabreaks <- parse_date_time3("19.03.2021 00:00",tz="Etc/GMT-1") + seq(0,5)*86400
mibreaks <- parse_date_time3("19.03.2021 12:00",tz="Etc/GMT-1") + seq(0,4)*86400
fill.col <- "grey70"


WS_TempUA <- WS2_dt[st >= parse_date_time3("18.03.2021 20:00",tz='Etc/GMT-1') & st <= parse_date_time3('21.03.2021 04:00',tz='Etc/GMT-1'),{
 	ggplot(.SD,aes(x=st,y=Temp)) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 11:19",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 11:40",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 01:08",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("20.03.2021 01:21",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_hline(yintercept = 0, colour="grey") +
 	geom_line() +
 	# geom_line(col="orange",lwd=1.2) +
 	scale_x_datetime(limits=limits_WSUA,expand=c(0.02,0.02)) +
 	ylim(-5,10) +
 	ylab("Temperature [°C]") +
 	xlab(NULL) +
 	theme_bw(base_size=18) +
 	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3),
	 axis.text.x = element_blank(), axis.ticks.x = element_blank(),plot.margin = unit(c(0.2, 1, -2, 0), "cm"))
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
 	# geom_line(data=.SD[WD < 180], aes(y=WD + 360),col='red') +
 	# geom_line(data=.SD[WD > 355], aes(y=WD - 360)) +
 	# geom_line(data=.SD[WD < 180],aes(y=WD),size=0.5) +
 	scale_x_datetime(limits=limits_WSUA,expand=c(0.02,0.02)) +
	scale_y_continuous(breaks = seq(0, 360, by=90), limits=c(-90,450),expand=c(-0.15,-0.15)) +
 	ylab("Wind direction [°]") +
 	xlab(NULL) +
 	theme_bw(base_size=18) +
 	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3), 
 		axis.text.x = element_blank(), axis.ticks.x = element_blank(),plot.margin = unit(c(-3, 1, -3, 0), "cm"))
}]



Fig_WS_UAUW <- Result[Sonic == 'SonicC' & st >= parse_date_time3("19.03.2021",tz='Etc/GMT-1') & st <= parse_date_time3('21.03.2021',tz='Etc/GMT-1'),{
	ggplot(.SD,aes(x=st,y=U_sonic)) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 11:19",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 11:40",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 01:08",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("20.03.2021 01:21",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_line(size=0.5) +
 	scale_x_datetime(limits=limits_WSUA,expand=c(0.02,0.02)) +
 	# scale_x_datetime(limits=limits_WSUA,expand=c(0.02,0.02),date_labels='%d.%m.%y %H:%M') +
	scale_y_continuous(breaks = seq(0, 6, by=2), limits=c(0,6),expand=c(0.02,0.02)) +
 	ylab(expression("Wind speed [m s"^-1*"]")) +
 	xlab(NULL) +
 	theme_bw(base_size=18) +
 	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3),plot.margin = unit(c(-3, 1, 0, 0), "cm"))
}]

Meteo_MK_WSUA <- ggarrange(WS_TempUA,Fig_WD_UAUW,Fig_WS_UAUW,ncol=1,align='hv')

ggsave(file.path(PfadFigures, "Meteo_WSUA.png"), Meteo_MK_WSUA, width = 30, height = 30/1.8, units="cm")


########################################################
########################################################
#####                                              #####
#####    WD deviation to WS700 and Upwind Sonic    #####
#####                                              #####
########################################################
########################################################


#######################################################
### Comparison of planar fit with two axis rotation ###
#######################################################

SonicC_1min_combined <- cbind(vanDik=SonicC_planar_1min_vanDik[,'WD'],Wilczak=SonicC_planar_1min_Wilczak[,'WD'],twoaxis=SonicC_1min[st(SonicC_planar_1min_vanDik),'WD'])
SonicC_10min_combined <- cbind(vanDik=SonicC_planar_10min_vanDik[,'WD'],Wilczak=SonicC_planar_10min_Wilczak[,'WD'],twoaxis=SonicC[st(SonicC_planar_10min_vanDik),'WD'])


plot(SonicC_1min_combined[,'vanDik'] ~  SonicC_1min_combined[,'Wilczak'],col='blue')
plot(SonicC_10min_combined[,'vanDik'] ~  SonicC_10min_combined[,'Wilczak'],add=TRUE)
abline(0,1,col='red',lwd=2)
# both planar fits are the same and also no difference in the average

plot(SonicC_10min_combined[,'vanDik'] ~  SonicC_10min_combined[,'twoaxis'])
abline(0,1,col='red',lwd=2)
# two axis is not different from planar fit


round(summary((SonicC_1min_combined$'vanDik' - SonicC_1min_combined$'Wilczak')) %% 360,2)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # -0.33    0.00    0.00    0.00    0.01    0.35 
round(summary(SonicC_1min_combined$'vanDik' - SonicC_1min_combined$'twoaxis'),2)
# maybe 0.3 deg difference. That can be ignored

#######################################################

##########################################################
### Differences to Background Sonic or weather station ###
##########################################################

####### 1 min #########

Sonic2_1min[,"MK"] <-  NA_real_
Sonic2_1min[MK,"MK"] <-  "MK"
Sonic2_1min[QV3,"MK"] <-  "QV3"

SonicAll_1min <- merge(merge(
        merge(
            merge(cbind(as.data.table(SonicB_1min),start_interval=st(SonicB_1min),end_interval = et(SonicB_1min)), cbind(as.data.table(SonicA_1min),start_interval=st(SonicA_1min),end_interval = et(SonicA_1min)), by = c('start_interval', 'end_interval'), suffixes = c('', '.a'))
            , cbind(as.data.table(SonicC_1min),start_interval=st(SonicC_1min),end_interval = et(SonicC_1min)), by = c('start_interval', 'end_interval'), suffixes = c('', '.c'))
                , cbind(as.data.table(Sonic2_1min),start_interval=st(Sonic2_1min),end_interval = et(Sonic2_1min)), by = c('start_interval', 'end_interval'), suffixes = c('.b', '.2'))
                    , cbind(as.data.table(WS_1min),start_interval=st(WS_1min),end_interval=et(WS_1min)), by = c('start_interval','end_interval'))


SonicAll_1min[, ':=' ( # Turbulence parameters do not make sense for such a short time interval. At least 10min needed.
	dWD.a = (((WD.a - WD.c) %% 360) + 180) %% 360 -180,
	dWD.2 = (((WD.2 - WD.c) %% 360) + 180) %% 360 -180,
	dWD.b = (((WD.b - WD.c) %% 360) + 180) %% 360 -180,
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


SonicAll_1min_cast <- melt(SonicAll_1min,id.vars = c('start_interval','MK','WD.c','U_sonic.c','WD_corr','WS_0480'),
	measure.vars = c('dWD.a','dWD.2','dWD.b','dU_sonic.a','dU_sonic.2','dU_sonic.b',
		'dWDWS.a','dWDWS.2','dWDWS.b','dWDWS.c','dU_sonicWS.a','dU_sonicWS.2','dU_sonicWS.b','dU_sonicWS.c'))

SonicAll_1min_cast[grep('.a$', variable), Sonic := '3DUA-50m']
SonicAll_1min_cast[grep('.2$', variable), Sonic := '3DUA-100m']
SonicAll_1min_cast[grep('.b$', variable), Sonic := '3DUA-150m']
SonicAll_1min_cast[grep('.c$', variable), Sonic := '3DUA-UW']

SonicAll_1min_cast$Sonic <- factor(SonicAll_1min_cast$Sonic, levels=c("3DUA-50m","3DUA-100m","3DUA-150m","3DUA-UW"),
	labels=c("3DUA-50m","3DUA-100m","3DUA-150m","3DUA-UW"))

SonicAll_1min_cast[grep('WD', variable), top_variable := 'WD']
SonicAll_1min_cast[grep('U_sonic', variable), top_variable := 'U_sonic']

SonicAll_1min_cast[grep('WS',variable),comparison := 'WS']
SonicAll_1min_cast[!grep('WS',variable),comparison := 'UW']



####### 10min #########

SonicB <- SonicB[,names(SonicB)!='MK']
Sonic2[,"MK"] <-  NA_real_
Sonic2[MK,"MK"] <-  "MK"
Sonic2[QV3,"MK"] <-  "QV3"

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


SonicAll_10min_cast <- melt(SonicAll_10min,id.vars = c('start_interval','MK','WD.c','U_sonic.c','Ustar.c','sVu.c','sWu.c','WD_corr','WS_0480'),
	measure.vars = c('dWD.a','dWD.2','dWD.b','dU_sonic.a','dU_sonic.2','dU_sonic.b',
		'dUstar.a','dUstar.2','dUstar.b','dsVu.a','dsVu.2','dsVu.b','dsWu.a','dsWu.2','dsWu.b',
		'dWDWS.a','dWDWS.2','dWDWS.b','dWDWS.c','dU_sonicWS.a','dU_sonicWS.2','dU_sonicWS.b','dU_sonicWS.c'))

SonicAll_10min_cast[grep('.a$', variable), Sonic := '3DUA-50m']
SonicAll_10min_cast[grep('.2$', variable), Sonic := '3DUA-100m']
SonicAll_10min_cast[grep('.b$', variable), Sonic := '3DUA-150m']
SonicAll_10min_cast[grep('.c$', variable), Sonic := '3DUA-UW']

SonicAll_10min_cast$Sonic <- factor(SonicAll_10min_cast$Sonic, levels=c("3DUA-50m","3DUA-100m","3DUA-150m","3DUA-UW"),
	labels=c("3DUA-50m","3DUA-100m","3DUA-150m","3DUA-UW"))

SonicAll_10min_cast$Sonic_ord <- factor(SonicAll_10min_cast$Sonic, levels=c("3DUA-UW","3DUA-50m","3DUA-100m","3DUA-150m"),
	labels=c("UA-UW","UA-50m","UA-100m","UA-150m"))

SonicAll_10min_cast$Sonic_ordh <- factor(SonicAll_10min_cast$Sonic, levels=c("3DUA-UW","3DUA-50m","3DUA-100m","3DUA-150m"),
	labels=c("UA-UW","UA-2.0h","UA-5.3h","UA-8.6h"))

SonicAll_10min_cast[grep('WD', variable), top_variable := 'WD']
SonicAll_10min_cast[grep('U_sonic', variable), top_variable := 'U_sonic']
SonicAll_10min_cast[grep('Ustar', variable), top_variable := 'Ustar']
SonicAll_10min_cast[grep('sVu', variable), top_variable := 'sVu']
SonicAll_10min_cast[grep('sWu', variable), top_variable := 'sWu']

SonicAll_10min_cast[grep('WS',variable),comparison := 'WS']
SonicAll_10min_cast[!grep('WS',variable),comparison := 'UW']


#####################################
### Difference to weather station ###
#####################################

### 1min ###

Fig_dWDWS_1min <- SonicAll_1min_cast[MK == "MK" & top_variable == 'WD' & comparison =='WS', {
    ggplot(.SD, aes(x=(WD_corr + 180) %% 360 - 180, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    geom_smooth(method='gam') +
    scale_x_continuous(limits=c(-10,120),breaks=seq(0,120,by=20)) +
    scale_y_continuous(limits=c(-30,30), breaks=seq(-30,30,by=10)) +
    xlab(expression("WD"["WS"]*" [°]")) +
    facet_grid(Sonic ~.) +
    ylab(expression("WD"["3DUA-DW"]*"  -  WD"["WS"]*"  [°]")) + 
    theme_bw(base_size = 18) +
    theme(strip.background = element_rect(fill="white"), panel.grid = element_blank())
}]


Fig_dUWS_1min <- SonicAll_1min_cast[MK == "MK" & top_variable == 'U_sonic' & comparison == 'WS', {
    ggplot(.SD, aes(x=WS_0480, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    geom_smooth(method='gam') +
    xlab(expression("WS"["3DUA-UW"]*" [°]")) +
    facet_grid(Sonic ~.) +
    ylab(expression("WD"["3DUA-DW"]*"  -  WD"["3DUA-UW"]*"  [°]")) + 
    theme_bw(base_size = 18) +
    theme(strip.background = element_rect(fill="white"), panel.grid = element_blank())
}]

### 10min ###

Fig_dWDWS_10min <- SonicAll_10min_cast[MK == "MK" & top_variable == 'WD' & comparison =='WS', {
    ggplot(.SD, aes(x=(WD_corr + 180) %% 360 - 180, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    geom_smooth(method='gam') +
    scale_x_continuous(limits=c(-10,120),breaks=seq(0,120,by=20)) +
    scale_y_continuous(limits=c(-30,30), breaks=seq(-30,30,by=10)) +
    xlab(expression("WD"["WS"]*" [°]")) +
    facet_grid(Sonic ~.) +
    ylab(expression("WD"["3DUA-DW"]*"  -  WD"["WS"]*"  [°]")) + 
    theme_bw(base_size = 18) +
    theme(strip.background = element_rect(fill="white"), panel.grid = element_blank())
}]


Fig_dUWS_10min <- SonicAll_10min_cast[MK == "MK" & top_variable == 'U_sonic' & comparison == 'WS', {
    ggplot(.SD, aes(x=WS_0480, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    geom_smooth(method='gam') +
    xlab(expression("WS"["3DUA-UW"]*" [°]")) +
    facet_grid(Sonic ~.) +
    ylab(expression("WD"["3DUA-DW"]*"  -  WD"["3DUA-UW"]*"  [°]")) + 
    theme_bw(base_size = 18) +
    theme(strip.background = element_rect(fill="white"), panel.grid = element_blank())
}]



##################################
### difference to upwind sonic ###
##################################

### 1min ### 

Fig_dWD_1min <- SonicAll_1min_cast[MK == "MK" & top_variable == 'WD' & comparison =='UW', {
    ggplot(.SD, aes(x=(WD.c + 180) %% 360 - 180, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    geom_smooth(method='gam') +
    scale_x_continuous(limits=c(-10,120),breaks=seq(0,120,by=20)) +
    scale_y_continuous(limits=c(-30,30), breaks=seq(-30,30,by=10)) +
    xlab(expression("WD"["3DUA-UW"]*" [°]")) +
    facet_grid(Sonic ~.) +
    ylab(expression("WD"["3DUA-DW"]*"  -  WD"["3DUA-UW"]*"  [°]")) + 
    theme_bw(base_size = 18) +
    theme(strip.background = element_rect(fill="white"), panel.grid = element_blank())
}]

ggsave(file.path(PfadFigures, "Fig_dWD_1min.png"), Fig_dWD_1min, width = 20, height = 20,units="cm")

Fig_dUsonic_WD_1min <- SonicAll_1min_cast[MK == "MK" & top_variable == 'U_sonic' & comparison == 'UW', {
    ggplot(.SD, aes(x=(WD.c + 180) %% 360 - 180, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    geom_smooth(method='gam') +
    scale_x_continuous(limits=c(-10,120),breaks=seq(0,120,by=20)) +
    scale_y_continuous(limits=c(-2,2), breaks=round(seq(-2,2,by=1),0)) +
    xlab(expression("WD"["3DUA-UW"]*" [°]")) +
    facet_grid(Sonic ~.) +
    ylab(expression("u "["3DUA-DW"]*"  -  u "["3DUA-UW"]*"  [m s"^-1*"]")) + 
    theme_bw(base_size = 18) +
    theme(strip.background = element_rect(fill="white"), panel.grid = element_blank())
}]

ggsave(file.path(PfadFigures, "Fig_dUsonic_WD_1min.png"), Fig_dUsonic_WD_1min, width = 20, height = 20,units="cm")


### 10min ###

Fig_dWD_10min <- SonicAll_10min_cast[MK == "MK" & top_variable == 'WD' & comparison == 'UW', {
    ggplot(.SD, aes(x=(WD.c + 180) %% 360 - 180, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    geom_smooth(method='gam') +
    scale_x_continuous(limits=c(-10,120),breaks=seq(0,120,by=20)) +
    scale_y_continuous(limits=c(-30,30), breaks=seq(-30,30,by=10)) +
    xlab(expression("WD"["3DUA-UW"]*" [°]")) +
    facet_grid(Sonic ~.) +
    ylab(expression("WD"["3DUA-DW"]*"  -  WD"["3DUA-UW"]*"  [°]")) + 
    theme_bw(base_size = 18) +
    theme(strip.background = element_rect(fill="white"), panel.grid = element_blank())
}]

ggsave(file.path(PfadFigures, "Fig_dWD_10min.png"), Fig_dWD_10min, width = 20, height = 20,units="cm")

Fig_dWD_10min_woT <- SonicAll_10min_cast[MK == "MK" & top_variable == 'WD' & comparison == 'UW', {
    ggplot(.SD, aes(x=(WD.c + 180) %% 360 - 180, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    # geom_vline(xintercept=c(35,40,45,50,55,60,65,70)) +
    scale_x_continuous(limits=c(-10,120),breaks=seq(0,120,by=20)) +
    scale_y_continuous(limits=c(-30,30), breaks=seq(-30,30,by=10)) +
    xlab(expression("WD"["UA-UW"]*" [°]")) +
    facet_grid(Sonic_ordh ~.) +
    ylab(expression("WD"["UA-DW"]*"  -  WD"["UA-UW"]*"  [°]")) + 
    theme_bw(base_size = 18) +
    theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_blank(), strip.background = element_rect(fill="white"))
    # theme(strip.background = element_rect(fill="white"),panel.grid = element_blank())
}]

ggsave(file.path(PfadFigures, "Fig_dWD_10min_woT.png"), Fig_dWD_10min_woT, width = 20, height = 20,units="cm")


Fig_dUsonic_WD_10min <- SonicAll_10min_cast[MK == "MK" & U_sonic.c > 0.0 & top_variable == 'U_sonic' & comparison == 'UW', {
    ggplot(.SD, aes(x=(WD.c + 180) %% 360 - 180, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    # geom_smooth(method='gam') +
    scale_x_continuous(limits=c(-10,120),breaks=seq(0,120,by=20)) +
    scale_y_continuous(limits=c(-2,2), breaks=round(seq(-2,2,by=1),0)) +
    xlab(expression("WD"["UA-UW"]*" [°]")) +
    facet_grid(Sonic_ordh ~.) +
    # ylab(expression(paste(italic('u'),' '['UA-DW']*'  -  ',italic('u'),' '['UA-UW']*'  [m s'^-1*']'))) + 
    ylab(expression(paste(italic('u')['UA-DW']*'  -  ',italic('u')['UA-UW']*'  [m s'^-1*']'))) + 
    theme_bw(base_size = 18) +
    theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_blank(), strip.background = element_rect(fill="white"))
}]

ggsave(file.path(PfadFigures, "Fig_dUsonic_WD_10min.png"), Fig_dUsonic_WD_10min, width = 20, height = 25,units="cm")


Fig_dUstar_WD_10min <- SonicAll_10min_cast[MK == "MK" & U_sonic.c > 0.0 & top_variable == 'Ustar' & comparison == 'UW', {
    ggplot(.SD, aes(x=(WD.c + 180) %% 360 - 180, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    # geom_smooth(method='gam') +
    scale_x_continuous(limits=c(-10,120),breaks=seq(0,120,by=20)) +
    scale_y_continuous(limits=c(-0.3,0.3), breaks=round(seq(-0.3,0.3,by=0.1),1)) +
    xlab(expression("WD"["UA-UW"]*" [°]")) +
    facet_grid(Sonic_ordh ~.) +
    # ylab(expression("u*"["UA-DW"]*"  -  u*"["UA-UW"]*"  [m s"^-1*"]")) + 
    ylab(expression(paste(italic('u'['*'])['UA-DW']*'  -  ',italic('u'['*'])['UA-UW']*'  [m s'^-1*']'))) + 
    theme_bw(base_size = 18) +
    theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_blank(), strip.background = element_rect(fill="white"))
}]

ggsave(file.path(PfadFigures, "Fig_dUstar_WD_10min.png"), Fig_dUstar_WD_10min, width = 20, height = 25,units="cm")


Fig_dsVu_WD_10min <- SonicAll_10min_cast[MK == "MK" & U_sonic.c > 0.0 & top_variable == 'sVu' & comparison == 'UW', {
    ggplot(.SD, aes(x=(WD.c + 180) %% 360 - 180, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    # geom_smooth(method='gam') +
    scale_x_continuous(limits=c(-10,120),breaks=seq(0,120,by=20)) +
    scale_y_continuous(limits=c(-6,6), breaks=round(seq(-6,6,by=2),0)) +
    xlab(expression("WD"["UA-UW"]*" [°]")) +
    facet_grid(Sonic_ordh ~.) +
    # ylab(expression("sVu"["UA-DW"]*"  -  sVu"["UA-UW"]*"  [m s"^-1*"]")) + 
  	ylab(expression(italic(sigma[v])*italic(u['*'])^-1*phantom()['UA-DW']*'  -  '*italic(sigma[v])*italic(u['*'])^-1*phantom()['UA-UW'])) +
    theme_bw(base_size = 18) +
    theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_blank(), strip.background = element_rect(fill="white"))
}]

ggsave(file.path(PfadFigures, "Fig_dsVu_WD_10min.png"), Fig_dsVu_WD_10min, width = 20, height = 25,units="cm")


Fig_dsWu_WD_10min <- SonicAll_10min_cast[MK == "MK" & U_sonic.c > 0.0 & top_variable == 'sWu' & comparison == 'UW', {
    ggplot(.SD, aes(x=(WD.c + 180) %% 360 - 180, y= value)) +
    geom_hline(yintercept=0) +
    geom_point(size=0.5) +
    # geom_smooth(method='gam') +
    scale_x_continuous(limits=c(-10,120),breaks=seq(0,120,by=20)) +
    scale_y_continuous(limits=c(-4,4), breaks=round(seq(-6,6,by=2),0)) +
    xlab(expression("WD"["UA-UW"]*" [°]")) +
    facet_grid(Sonic_ordh ~.) +
    # ylab(expression("sWu"["UA-DW"]*"  -  sWu"["UA-UW"]*"  [m s"^-1*"]")) +
   	ylab(expression(italic(sigma[w])*italic(u['*'])^-1*phantom()['UA-DW']*'  -  '*italic(sigma[w])*italic(u['*'])^-1*phantom()['UA-UW'])) +
    theme_bw(base_size = 18) +
    theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_blank(), strip.background = element_rect(fill="white"))
}]

ggsave(file.path(PfadFigures, "Fig_dsWu_WD_10min.png"), Fig_dsWu_WD_10min, width = 20, height = 25,units="cm")


Fig_turbulence <- ggarrange(Fig_dUsonic_WD_10min,Fig_dUstar_WD_10min,Fig_dsVu_WD_10min,Fig_dsWu_WD_10min,ncol=2,nrow=2,align="hv")
ggsave(file.path(PfadFigures,'Fig_turbulence_paper.png'),Fig_turbulence,width=30, height=30, units='cm')

## calculate the mean offset for sVu and sWu for the 10min intervals
SonicAll_10min_cast[MK == "MK" & top_variable == 'sVu' & comparison == 'UW', mean(value,na.rm=TRUE), by = Sonic]





iWD <- c(15,100)
plot(Sources,Sensors_MK,sensors.text.arg=list(labels=''),lines.args=list(lwd=2,col="black",lty=3),points.args = list(pch = 20, cex = 2, col ="black"),polygon.args=list(col='grey60'),sources.text.arg=list(labels=''))
plot(Sonics_MK,sensors.text.arg=list(labels=''),points.args = list(pch = '*', cex = 3, col ="black"),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
lines_sec2xy(Sonics_MK,"SonicA",1,iWD[1],lty=2)
lines_sec2xy(Sonics_MK,"SonicA",1,iWD[2],lty=2)
lines_sec2xy(Sonics_MK,"Sonic2",1,iWD[1],lty=2)
lines_sec2xy(Sonics_MK,"Sonic2",1,iWD[2],lty=2)
lines_sec2xy(Sonics_MK,"SonicB",1,iWD[1],lty=2)
lines_sec2xy(Sonics_MK,"SonicB",1,iWD[2],lty=2)



PlotOnStaticMap(STO_Map)
plot(Sensors_MK_xy,sensors.text.arg=list(labels=''),lines.args=list(lwd=2,col="orange",lty=3),points.args = list(pch = 20, cex = 2, col ="orange"),add=TRUE)
plot(Sonics_MK_xy,sensors.text.arg=list(labels=''),points.args = list(pch = '*', cex = 3, col ="orange"),add=TRUE)
lines_sec2xy(Sonics_MK_xy,"SonicA",1,iWD[1],lty=2)
lines_sec2xy(Sonics_MK_xy,"SonicA",1,iWD[2],lty=2)
lines_sec2xy(Sonics_MK_xy,"Sonic2",1,iWD[1],lty=2)
lines_sec2xy(Sonics_MK_xy,"Sonic2",1,iWD[2],lty=2)
lines_sec2xy(Sonics_MK_xy,"SonicB",1,iWD[1],lty=2)
lines_sec2xy(Sonics_MK_xy,"SonicB",1,iWD[2],lty=2)
lines_sec2xy(Sonics_MK_xy,"SonicC",1,iWD[1],lty=2)
lines_sec2xy(Sonics_MK_xy,"SonicC",1,iWD[2],lty=2)



#########################################
#########################################
#####                               #####
#####    Scatterplots Sonic data    #####
#####                               #####
#########################################
#########################################


Sonic2[,"MK"] <-  NA_real_
Sonic2[,"Gas"] <-  NA_real_
Sonic2[MK,"MK"] <-  "MK"
Sonic2[QV3,"MK"] <-  "QV3"
Sonic2[indMK_MFC,"Gas"] <-  "Gas"



SonicAll <- merge(
        merge(
            merge(cbind(as.data.table(SonicB),start_interval=st(SonicB),end_interval = et(SonicB)), cbind(as.data.table(SonicA),start_interval=st(SonicA),end_interval = et(SonicA)), by = c('start_interval', 'end_interval'), suffixes = c('', '_50m'))
            , cbind(as.data.table(SonicC),start_interval=st(SonicC),end_interval = et(SonicC)), by = c('start_interval', 'end_interval'), suffixes = c('', '_UW')
            ), cbind(as.data.table(Sonic2),start_interval=st(Sonic2),end_interval = et(Sonic2)), by = c('start_interval', 'end_interval'), suffixes = c('_150m', '_100m')
        )


### Ustar
# entire MK
png(file.path(PfadFigures,"Scat_Ustar_MK.png"),width=30,height=30/2,unit="cm",res=600)
SonicAll[MK == "MK", {
    pairs(.SD, panel = function(x,y,...){
            m <- deming::deming(y~x)$coefficients
            points(x,y,...)
            # grid()
            abline(0,1)
            abline(m[1],m[2],col='red',lwd=2)
            },asp=1,xlim=c(0,0.51),ylim=c(0,0.51))
}, .SDcols = c('Ustar_50m', 'Ustar_100m', 'Ustar_150m', 'Ustar_UW')]
dev.off()

# Release
png(file.path(PfadFigures,"Scat_Ustar_Release.png"),width=30,height=30/2,unit="cm",res=600)
SonicAll[Gas == "Gas", {
    pairs(.SD, panel = function(x,y,...){
            m <- deming::deming(y~x)$coefficients
            points(x,y,...)
            # grid()
            abline(0,1)
            abline(m[1],m[2],col='red',lwd=2)
            },asp=1,xlim=c(0,0.51),ylim=c(0,0.51))
}, .SDcols = c('Ustar_50m', 'Ustar_100m', 'Ustar_150m', 'Ustar_UW')]
dev.off()



SonicAll_melt <- melt(SonicAll, id=c('WD_UW','Ustar_UW','U_sonic_UW','sUu_UW','sWu_UW','z0_UW','MK','Gas'),measure.vars=c('WD_50m','WD_100m','WD_150m'
	,'Ustar_50m','Ustar_100m','Ustar_150m'
	,'U_sonic_50m','U_sonic_100m','U_sonic_150m'
	,'sUu_50m','sUu_100m','sUu_150m'
	,'sWu_50m','sWu_100m','sWu_150m'
	,'z0_50m','z0_100m','z0_150m'))


SonicAll_melt[variable %in% c('WD_50m','WD_150m','WD_100m'),{
	ggplot(.SD, aes(x= (WD_UW +180) %%360 -180 ,y=(value+180) %% 360 - 180,colour=variable)) +
	geom_abline(slope=1,intercept=0,lwd=1.2) +
	geom_point() +
	geom_smooth(method='gam') +
	scale_x_continuous(limits=c(-20,120)) +
	scale_y_continuous(limits=c(-20,120)) +
	facet_grid(variable ~ .) +
	xlab('WD upwind sonic') +
	ylab('WD downwind sonic') +
	theme_bw() +
	theme(legend.position = 'none')
}]


Scat_Ustar_MK <- SonicAll_melt[variable %in% c('Ustar_50m','Ustar_150m','Ustar_100m') & MK == 'MK',{
	ggplot(.SD, aes(x= Ustar_UW,y=value,colour=variable)) +
	geom_abline(slope=1,intercept=0,lwd=1.2) +
	geom_point() +
	geom_smooth(method='lm') +
	facet_grid(variable ~ .) +
	ylab('Ustar (MK)') +
	theme_bw() +
	theme(legend.position = 'none')
}]

Scat_Ustar_Release <- SonicAll_melt[variable %in% c('Ustar_50m','Ustar_150m','Ustar_100m') & Gas == 'Gas',{
	ggplot(.SD, aes(x= Ustar_UW,y=value,colour=variable)) +
	geom_abline(slope=1,intercept=0,lwd=1.2) +
	geom_point() +
	geom_smooth(method='lm') +
	facet_grid(variable ~ .) +
	ylab('Ustar (Gas release)') +
	theme_bw() +
	theme(legend.position = 'none')
}]


Scat_u_MK <- SonicAll_melt[variable %in% c('U_sonic_50m','U_sonic_150m','U_sonic_100m') & MK == 'MK',{
	ggplot(.SD, aes(x= U_sonic_UW,y=value,colour=variable)) +
	geom_abline(slope=1,intercept=0,lwd=1.2) +
	geom_point() +
	geom_smooth(method='lm') +
	facet_grid(variable ~ .) +
	ylab('U_sonic (MK)') +
	theme_bw() +
	theme(legend.position = 'none')
}]


Scat_u_Release <- SonicAll_melt[variable %in% c('U_sonic_50m','U_sonic_150m','U_sonic_100m') & Gas == 'Gas',{
	ggplot(.SD, aes(x= U_sonic_UW,y=value,colour=variable)) +
	geom_abline(slope=1,intercept=0,lwd=1.2) +
	geom_point() +
	geom_smooth(method='lm') +
	facet_grid(variable ~ .) +
	ylab('U_sonic (Gas release)') +
	theme_bw() +
	theme(legend.position = 'none')
}]

Scat_Ustar <- ggarrange(Scat_Ustar_MK,Scat_Ustar_Release)
Scat_u <- ggarrange(Scat_u_MK,Scat_u_Release)
Scat_MK <- ggarrange(Scat_Ustar_MK,Scat_u_MK)
Scat_Release <- ggarrange(Scat_Ustar_Release,Scat_u_Release)

ggsave(file.path(PfadFigures,"Scat_Ustar.png"),Scat_Ustar,width=30, height=30/2,units="cm")
ggsave(file.path(PfadFigures,"Scat_u.png"),Scat_u,width=30, height=30/2,units="cm")
ggsave(file.path(PfadFigures,"Scat_MK.png"),Scat_MK,width=30, height=30/2,units="cm")
ggsave(file.path(PfadFigures,"Scat_Release.png"),Scat_Release,width=30, height=30/2,units="cm")



Scat_sUu_MK <- SonicAll_melt[variable %in% c('sUu_50m','sUu_150m','sUu_100m') & MK == 'MK',{
	ggplot(.SD, aes(x= sUu_UW,y=value,colour=variable)) +
	geom_abline(slope=1,intercept=0,lwd=1.2) +
	geom_point() +
	ylim(2,10) +
	geom_smooth(method='lm') +
	facet_grid(variable ~ .) +
	ylab('sUu (MK)') +
	theme_bw() +
	theme(legend.position = 'none')
}]

Scat_sUu_Release <- SonicAll_melt[variable %in% c('sUu_50m','sUu_150m','sUu_100m') & Gas == 'Gas',{
	ggplot(.SD, aes(x= sUu_UW,y=value,colour=variable)) +
	geom_abline(slope=1,intercept=0,lwd=1.2) +
	geom_point() +
	ylim(2,10) +
	geom_smooth(method='lm') +
	facet_grid(variable ~ .) +
	ylab('sUu (Gas release)') +
	theme_bw() +
	theme(legend.position = 'none')
}]

# sieht irgendwie konstant aus. 

Scat_sWu_Release <- SonicAll_melt[variable %in% c('sWu_50m','sWu_150m','sWu_100m') & Gas == 'Gas',{
	ggplot(.SD, aes(x= sWu_UW,y=value,colour=variable)) +
	geom_abline(slope=1,intercept=0,lwd=1.2) +
	geom_point() +
	ylim(1.1,1.7) +
	xlim(1.1,1.7) +
	geom_smooth(method='lm') +
	facet_grid(variable ~ .) +
	ylab('sWu (Gas release)') +
	theme_bw() +
	theme(legend.position = 'none')
}]


Scat_z0_Release <- SonicAll_melt[variable %in% c('z0_50m','z0_150m','z0_100m') & Gas == 'Gas',{
	ggplot(.SD, aes(x= z0_UW,y=value,colour=variable)) +
	geom_abline(slope=1,intercept=0,lwd=1.2) +
	geom_point() +
	ylim(0,0.07) +
	xlim(0,0.07) +
	geom_smooth(method='lm') +
	facet_grid(variable ~ .) +
	ylab('z0 (Gas release)') +
	theme_bw() +
	theme(legend.position = 'none')
}]



Scat_u_Release <- SonicAll_melt[variable %in% c('U_sonic_50m','U_sonic_150m','U_sonic_100m') & Gas == 'Gas',{
	ggplot(.SD, aes(x= U_sonic_UW,y=value,colour=variable)) +
	geom_abline(slope=1,intercept=0,lwd=1.2) +
	geom_point() +
	geom_smooth(method='lm') +
	facet_grid(variable ~ .) +
	ylab('U_sonic (Gas release)') +
	theme_bw() +
	theme(legend.position = 'none')
}]



#################################
#################################
#####                       #####
#####    Histogram plots    #####
#####                       #####
#################################
#################################

##############################
### Histogram of u* values ###
##############################
par(mfrow=c(2,2))
hist(na.omit(SonicC[indMK_MFC]$'Ustar'),main=NA,breaks=seq(0.05,0.55,by=0.025),col='#16D466')
hist(na.omit(SonicA[indMK_MFC]$'Ustar'),main=NA,breaks=seq(0.05,0.55,by=0.025),col='green')
hist(na.omit(Sonic2[indMK_MFC]$'Ustar'),main=NA,breaks=seq(0.05,0.55,by=0.025),col='#2C39EA22')
hist(na.omit(SonicB[indMK_MFC]$'Ustar'),main=NA,breaks=seq(0.05,0.55,by=0.025),col='#D2101022')

# on top of each other without the closest Sonic
par(mfrow=c(1,1))
hist(na.omit(SonicC[indMK_MFC]$'Ustar'),main=NA,breaks=seq(0.05,0.55,by=0.025),col='#D4CD16')
hist(na.omit(Sonic2[indMK_MFC]$'Ustar'),main=NA,breaks=seq(0.05,0.55,by=0.025),add=TRUE,col='#84F0E760')
hist(na.omit(SonicB[indMK_MFC]$'Ustar'),main=NA,breaks=seq(0.05,0.55,by=0.025),add=TRUE,col='#D2101070')

##############################

dt_Sonic_Gas <- rbind(as.data.table(cbind(SonicC[indMK_MFC,c('U_sonic','Ustar','z0','WD','sd_WD','L','sUu','sVu','sWu','d')],Sonic = 'UW'))
			,as.data.table(cbind(SonicA[indMK_MFC,c('U_sonic','Ustar','z0','WD','sd_WD','L','sUu','sVu','sWu','d')],Sonic = '50m'))
			,as.data.table(cbind(Sonic2[indMK_MFC,c('U_sonic','Ustar','z0','WD','sd_WD','L','sUu','sVu','sWu','d')],Sonic = '100m'))
			,as.data.table(cbind(SonicB[indMK_MFC,c('U_sonic','Ustar','z0','WD','sd_WD','L','sUu','sVu','sWu','d')],Sonic = '150m')))

dt_Sonic_MC <- rbind(as.data.table(cbind(SonicC[MK,c('U_sonic','Ustar','z0','WD','sd_WD','L','sUu','sVu','sWu','d')],Sonic = 'UW'))
			,as.data.table(cbind(SonicA[MK,c('U_sonic','Ustar','z0','WD','sd_WD','L','sUu','sVu','sWu','d')],Sonic = '50m'))
			,as.data.table(cbind(Sonic2[MK,c('U_sonic','Ustar','z0','WD','sd_WD','L','sUu','sVu','sWu','d')],Sonic = '100m'))
			,as.data.table(cbind(SonicB[MK,c('U_sonic','Ustar','z0','WD','sd_WD','L','sUu','sVu','sWu','d')],Sonic = '150m')))

dt_Sonic_Gas[,MC := 'Release']
dt_Sonic_MC[,MC := 'MC']

dt_Sonic <- rbind(dt_Sonic_Gas,dt_Sonic_MC)

### only during CH4 release in the MC ###

### Ustar

dt_Sonic[Sonic != '50m',{
	ggplot(na.omit(.SD), aes(x=Ustar,fill=Sonic,colour=Sonic,after_stat(density))) +
	geom_histogram(,position='dodge',colour='black') +
	geom_density(alpha=0.01,lwd=1) +
	scale_colour_manual(values = colSonics,aesthetics=c('colour','fill')) +
	facet_grid(. ~ MC) +
	theme_bw()
}]

Hist_Ustar <- dt_Sonic[,{
	ggplot(na.omit(.SD), aes(x=Ustar,fill=Sonic,colour=Sonic,after_stat(density))) +
	geom_histogram(,position='dodge',colour='black') +
	geom_density(alpha=0.01,lwd=1) +
	scale_colour_manual(values = colSonics,aesthetics=c('colour','fill')) +
	facet_grid(. ~ MC) +
	theme_bw()
}]

Hist_Ustar_vline <- dt_Sonic[,{
	ggplot(na.omit(.SD), aes(x=Ustar,fill=Sonic,colour=Sonic,after_stat(density))) +
	geom_histogram(,position='dodge',colour='black') +
	geom_density(alpha=0.01,lwd=1) +
	geom_vline(xintercept = c(0.15,0.25,0.35), lwd=1.3) +
	scale_colour_manual(values = colSonics,aesthetics=c('colour','fill')) +
	facet_grid(. ~ MC) +
	theme_bw()
}]

### U_sonic

dt_Sonic[Sonic != '50m',{
	ggplot(.SD, aes(x=U_sonic,fill=Sonic,colour=Sonic,after_stat(density))) +
	geom_histogram(position='dodge',colour='black') +
	geom_density(alpha=0.01,lwd=1) +
	facet_grid(. ~ MC) +
	theme_bw()
}]

Hist_Usonic <- dt_Sonic[,{
	ggplot(na.omit(.SD), aes(x=U_sonic,fill=Sonic,colour=Sonic,after_stat(density))) +
	geom_histogram(,position='dodge',colour='black') +
	geom_density(alpha=0.01,lwd=1) +
	scale_colour_manual(values = colSonics,aesthetics=c('colour','fill')) +
	facet_grid(. ~ MC) +
	theme_bw()
}]


### z0

Hist_z0 <- dt_Sonic[,{
	ggplot(na.omit(.SD), aes(x=z0,fill=Sonic,colour=Sonic,after_stat(density))) +
	geom_histogram(,position='dodge',colour='black') +
	geom_density(alpha=0.01,lwd=1) +
	scale_colour_manual(values = colSonics,aesthetics=c('colour','fill')) +
	xlim(0,0.15) +
	facet_grid(. ~ MC) +
	theme_bw()
}]


## WD
Hist_WD <- dt_Sonic[,{
	ggplot(na.omit(.SD), aes(x=(180+WD) %% 360 -180,fill=Sonic,colour=Sonic,after_stat(density))) +
	geom_histogram(position='dodge',colour='black') +
	geom_density(alpha=0.01,lwd=1) +
	scale_x_continuous(limits=c(-30,130),breaks=seq(-20,120,20)) +
	xlab('Wind direction') +
	scale_colour_manual(values = colSonics,aesthetics=c('colour','fill')) +
	facet_grid(. ~ MC) +
	theme_bw()
}]



## sd_WD
Hist_sdWD <- dt_Sonic[,{
	ggplot(na.omit(.SD), aes(x=sd_WD,fill=Sonic,colour=Sonic,after_stat(density))) +
	geom_histogram(position='dodge',colour='black') +
	geom_density(alpha=0.01,lwd=1) +
	scale_colour_manual(values = colSonics,aesthetics=c('colour','fill')) +
	xlim(0,35) +
	facet_grid(. ~ MC) +
	theme_bw()
}]


## L
Hist_L <- dt_Sonic[,{
	ggplot(na.omit(.SD), aes(x=(1/L),fill=Sonic,colour=Sonic,after_stat(density))) +
	geom_histogram(position='dodge',colour='black') +
	geom_density(alpha=0.01,lwd=1) +
	scale_colour_manual(values = colSonics,aesthetics=c('colour','fill')) +
	xlim(-0.6,0.6) +
	facet_grid(. ~ MC) +
	theme_bw()
}]

## sUu
Hist_sUu <- dt_Sonic[,{
	ggplot(na.omit(.SD), aes(x=sUu,fill=Sonic,colour=Sonic,after_stat(density))) +
	geom_histogram(position='dodge',colour='black') +
	geom_density(alpha=0.01,lwd=1) +
	scale_colour_manual(values = colSonics,aesthetics=c('colour','fill')) +
	xlim(1,7.5) +
	facet_grid(. ~ MC) +
	theme_bw()
}]

## sVu
Hist_sVu <- dt_Sonic[,{
	ggplot(na.omit(.SD), aes(x=sVu,fill=Sonic,colour=Sonic,after_stat(density))) +
	geom_histogram(position='dodge',colour='black') +
	geom_density(alpha=0.01,lwd=1) +
	scale_colour_manual(values = colSonics,aesthetics=c('colour','fill')) +
	xlim(1,7.5) +
	facet_grid(. ~ MC) +
	theme_bw()
}]

## sWu
Hist_sWu <- dt_Sonic[,{
	ggplot(na.omit(.SD), aes(x=sWu,fill=Sonic,colour=Sonic,after_stat(density))) +
	geom_histogram(position='dodge',colour='black') +
	geom_density(alpha=0.01,lwd=1) +
	scale_colour_manual(values = colSonics,aesthetics=c('colour','fill')) +
	xlim(0.5,4) +
	facet_grid(. ~ MC) +
	theme_bw()
}]


ggsave(file.path(PfadFigures,"Hist_Ustar.png"),Hist_Ustar,width=30, height=30/2,units="cm")
ggsave(file.path(PfadFigures,"Hist_Ustar_vline.png"),Hist_Ustar_vline,width=30, height=30/2,units="cm")
ggsave(file.path(PfadFigures,"Hist_Usonic.png"),Hist_Usonic,width=30, height=30/2,units="cm")
ggsave(file.path(PfadFigures,"Hist_z0.png"),Hist_z0,width=30, height=30/2,units="cm")
ggsave(file.path(PfadFigures,"Hist_WD.png"),Hist_WD,width=30, height=30/2,units="cm")
ggsave(file.path(PfadFigures,"Hist_sdWD.png"),Hist_sdWD,width=30, height=30/2,units="cm")
ggsave(file.path(PfadFigures,"Hist_L.png"),Hist_L,width=30, height=30/2,units="cm")
ggsave(file.path(PfadFigures,"Hist_sUu.png"),Hist_sUu,width=30, height=30/2,units="cm")
ggsave(file.path(PfadFigures,"Hist_sVu.png"),Hist_sVu,width=30, height=30/2,units="cm")
ggsave(file.path(PfadFigures,"Hist_sWu.png"),Hist_sWu,width=30, height=30/2,units="cm")


#################################
#################################
#####                       #####
#####    Wind field plot    #####
#####                       #####
#################################
#################################

SonicB[,"MK"] <-  NA_real_
SonicB[MK,"MK"] <-  "MK"
SonicB[QV3,"MK"] <-  "QV3"

SonicWF <- merge(merge(merge(cbind(as.data.table(SonicB[,c('WD','Ustar','U_sonic')]),st=st(SonicB)),
		cbind(as.data.table(SonicA[,c('WD','Ustar','U_sonic')]),st=st(SonicA)), by = c('st'),all=TRUE, suffixes = c('', '.a'))
            , cbind(as.data.table(SonicC[,c('WD','Ustar','U_sonic')]),st=st(SonicC)), by = c('st'),all=TRUE, suffixes = c('', '.c'))
                , cbind(as.data.table(Sonic2[,c('WD','Ustar','U_sonic')]),st=st(Sonic2)), by = c('st'),all=TRUE, suffixes = c('.b', '.2'))
                    # , cbind(as.data.table(WS_10min),st=st(WS_10min)), by = c('st','et'))

SonicWF[st > parse_date_time3('18.03.2021 11:00',tz='Etc/GMT-1') & st < parse_date_time3('21.03.2022 14:00', tz='Etc/GMT-1'),Campaign := 'MK']
SonicWF[st %in% Result[Campaign == 'MK' & !is.na(Q_MFC),unique(st)], Gas := 'YES']



wind_field <- function(GPS_sonic,data, r_circle=0, breaks = NA_real_, delta_wd=22.5,scale=10,start=0,wd='WD',ws='WS',a_lwd=3,cex_points=0.5
						,draw.circle=TRUE,draw.points=TRUE,draw.lines=FALSE,col_points='grey60',col_lines='grey60',col_circle='grey60',line_lwd=1
						,scale_line=10){
# browser()
if(is.na(delta_wd)){
	delta_wd <- 360/breaks
}
if((360 / delta_wd) %% 1 != 0) stop("360 degrees needs to be a multiplier of 'delta_wd'")

centre <- as.data.frame(GPS_sonic[,4:5])
# mean wind direction and mean wind speed for the different bins
mWD_dt <- lapply((start+seq(0, 360-delta_wd, by = delta_wd) %% 360), function(x) {
	if((x + delta_wd) <= 360){
		sub <- data[WD.c >= x & WD.c < x+delta_wd]} else {
		sub <- data[WD.c >= x | WD.c < x]
		}
	data.table('WD_UW' = (atan2(mean(sub[,U_sonic.c * sin(WD.c * pi/180)], na.rm=TRUE), mean(sub[,U_sonic.c * cos(WD.c * pi/180)],na.rm=TRUE)) * 180/pi) %% 360
	, 'WD_50m' = (atan2(mean(sub[,U_sonic.a * sin(WD.a * pi/180)], na.rm=TRUE), mean(sub[,U_sonic.a * cos(WD.a * pi/180)],na.rm=TRUE)) * 180/pi) %% 360
	, 'WD_100m' = (atan2(mean(sub[,U_sonic.2 * sin(WD.2 * pi/180)], na.rm=TRUE), mean(sub[,U_sonic.2 * cos(WD.2 * pi/180)],na.rm=TRUE)) * 180/pi) %% 360
	, 'WD_150m' = (atan2(mean(sub[,U_sonic.b * sin(WD.b * pi/180)], na.rm=TRUE), mean(sub[,U_sonic.b * cos(WD.b * pi/180)],na.rm=TRUE)) * 180/pi) %% 360
	, 'WS_UW' = sub[,mean(U_sonic.c,na.rm=TRUE)*scale]
	, 'WS_50m' = sub[,mean(U_sonic.a,na.rm=TRUE)*scale]
	, 'WS_100m' = sub[,mean(U_sonic.2,na.rm=TRUE)*scale]
	, 'WS_150m' = sub[,mean(U_sonic.b,na.rm=TRUE)*scale]
	,	bin = x+delta_wd/2)
})
mWD_dt <- rbindlist(mWD_dt)

# defining the arrowhead points
a_points <- lapply((seq(0,360-delta_wd,by=delta_wd)+delta_wd/2) %% 360,function(x)
	data.table(
	x0_UW=r_circle * sin(x*pi/180) + centre[4,1], y0_UW=r_circle*cos(x*pi/180) + centre[4,2]
	,x0_50m=r_circle * sin(x*pi/180) + centre[2,1], y0_50m=r_circle*cos(x*pi/180) + centre[2,2]
	,x0_100m=r_circle * sin(x*pi/180) + centre[1,1], y0_100m=r_circle*cos(x*pi/180) + centre[1,2]
	,x0_150m=r_circle * sin(x*pi/180) + centre[3,1], y0_150m=r_circle*cos(x*pi/180) + centre[3,2]
	,bin=x))
	# ,bin=(90-x) %% 360))
a_points <- rbindlist(a_points)

# merge both dt togther
dt_wf <- merge(a_points,mWD_dt)

dt_wf[, ':='(
  x1_UW = WS_UW * sin(WD_UW*pi/180) + x0_UW
 ,y1_UW = WS_UW * cos(WD_UW*pi/180) + y0_UW
 ,x1_50m = WS_50m * sin(WD_50m*pi/180) + x0_50m
 ,y1_50m = WS_50m * cos(WD_50m*pi/180) + y0_50m
 ,x1_100m = WS_100m * sin(WD_100m*pi/180) + x0_100m
 ,y1_100m = WS_100m * cos(WD_100m*pi/180) + y0_100m
 ,x1_150m = WS_150m * sin(WD_150m*pi/180) + x0_150m
 ,y1_150m = WS_150m * cos(WD_150m*pi/180) + y0_150m
 ,x2_UW = WD_UW/WD_UW * (10*scale_line * sin(bin*pi/180) + x0_UW)
 ,y2_UW = WD_UW/WD_UW * (10*scale_line * cos(bin*pi/180) + y0_UW)
 ,x2_50m = WD_UW/WD_UW * (10*scale_line * sin(bin*pi/180) + x0_50m)
 ,y2_50m = WD_UW/WD_UW * (10*scale_line * cos(bin*pi/180) + y0_50m)
 ,x2_100m = WD_UW/WD_UW * (10*scale_line * sin(bin*pi/180) + x0_100m)
 ,y2_100m = WD_UW/WD_UW * (10*scale_line * cos(bin*pi/180) + y0_100m)
 ,x2_150m = WD_UW/WD_UW * (10*scale_line * sin(bin*pi/180) + x0_150m)
 ,y2_150m = WD_UW/WD_UW * (10*scale_line * cos(bin*pi/180) + y0_150m)
)]


if(draw.circle == TRUE){
	plotrix::draw.circle(centre[2,1],centre[2,2],r_circle,border=col_circle)
	plotrix::draw.circle(centre[3,1],centre[3,2],r_circle,border=col_circle)
	plotrix::draw.circle(centre[4,1],centre[4,2],r_circle,border=col_circle)
	plotrix::draw.circle(centre[1,1],centre[1,2],r_circle,border=col_circle)
}
if(draw.points == TRUE){
	points(a_points[,.(x0_UW,y0_UW)],pch=20,cex=cex_points,col=col_points)
	points(a_points[,.(x0_50m,y0_50m)],pch=20,cex=cex_points,col=col_points)
	points(a_points[,.(x0_100m,y0_100m)],pch=20,cex=cex_points,col=col_points)
	points(a_points[,.(x0_150m,y0_150m)],pch=20,cex=cex_points,col=col_points)
}
if(draw.lines == TRUE){
arrows(x0=dt_wf[,x2_UW],y0=dt_wf[,y2_UW],x1=dt_wf[,x0_UW],y1=dt_wf[,y0_UW],length=0.1,lwd=line_lwd, col=col_lines)
arrows(x0=dt_wf[,x2_50m],y0=dt_wf[,y2_50m],x1=dt_wf[,x0_50m],y1=dt_wf[,y0_50m],length=0.1,lwd=line_lwd, col=col_lines)
arrows(x0=dt_wf[,x2_100m],y0=dt_wf[,y2_100m],x1=dt_wf[,x0_100m],y1=dt_wf[,y0_100m],length=0.1,lwd=line_lwd, col=col_lines)
arrows(x0=dt_wf[,x2_150m],y0=dt_wf[,y2_150m],x1=dt_wf[,x0_150m],y1=dt_wf[,y0_150m],length=0.1,lwd=line_lwd, col=col_lines)
}

colours <- rainbow(mWD_dt[!is.na(WD_UW),.N] - 1)
arrows(x0=dt_wf[,x1_UW],y0=dt_wf[,y1_UW],x1=dt_wf[,x0_UW],y1=dt_wf[,y0_UW],length=0.1,lwd=a_lwd, col=colours)
arrows(x0=dt_wf[,x1_50m],y0=dt_wf[,y1_50m],x1=dt_wf[,x0_50m],y1=dt_wf[,y0_50m],length=0.1,lwd=a_lwd, col=colours)
arrows(x0=dt_wf[,x1_100m],y0=dt_wf[,y1_100m],x1=dt_wf[,x0_100m],y1=dt_wf[,y0_100m],length=0.1,lwd=a_lwd, col=colours)
arrows(x0=dt_wf[,x1_150m],y0=dt_wf[,y1_150m],x1=dt_wf[,x0_150m],y1=dt_wf[,y0_150m],length=0.1,lwd=a_lwd, col=colours)
}

plot(Sources,Sensors_MK,xlim=c(583700,584100),ylim=c(210050,210400),sensors.text.arg=list(labels=''),lines.args=list(lwd=2,col="black",lty=3),points.args = list(pch = 20, cex = 2, col ="black"),polygon.args=list(col='grey60'),sources.text.arg=list(labels=''))
plot(Sonics_MK,sensors.text.arg=list(labels=''),points.args = list(pch = '*', cex = 3, col ="black"),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
# wind_field(GPS_sonic=Sonics_MK,data=SonicWF[Campaign=='MK'],delta_wd=6,r_circle=20)
wind_field(GPS_sonic=Sonics_MK,data=SonicWF[Gas=='YES'],delta_wd=10,r_circle=20,draw.points=FALSE,draw.lines=TRUE,scale_line=5,line_lwd=2)


plot(1,xlim=c(583700,584100),ylim=c(210050,210400),type='n',xaxt="n",yaxt="n", xlab="",ylab="",asp=1)
plot(Sources,polygon.args=list(col='grey60'),sources.text.arg=list(labels=''),add=TRUE)
plot(Baum,sensors.text.args=list(labels=""),points.args = list(pch = 21, cex = 3.8, col ="grey60",bg="grey95",lwd=2),add=TRUE)
plot(Sensors_MK,sensors.text.arg=list(labels=''),lines.args=list(lwd=2,col="black",lty=3),points.args = list(pch = 20, cex = 2, col ="black"),add=TRUE)
plot(Sonics_MK,sensors.text.arg=list(labels=''),points.args = list(pch = '*', cex = 3, col ="black"),add=TRUE)
# wind_field(GPS_sonic=Sonics_MK,data=SonicWF[Campaign=='MK'],delta_wd=6,r_circle=20)
wind_field(GPS_sonic=Sonics_MK,data=SonicWF[Gas=='YES'],delta_wd=2,r_circle=0,draw.points=TRUE,draw.lines=FALSE,scale_line=5,line_lwd=2)




rainbow(360/5 - 1)


#####################################
### Wind vector times series plot ###
#####################################

Result[Campaign == 'MK', stmin := (as.numeric(st) - as.numeric(st[1]))/60]

ggplot(Result[Sonic=='SonicC' & Campaign == 'MK'], aes(x = st, y = 0)) +
	  geom_segment(aes(xend = st * sin(WD * pi / 180) * U_sonic,
                   yend = cos(WD * pi / 180) * U_sonic),
               arrow = arrow(length = unit(0.15, "inches")), size = 1) +
  labs(x = "Time", y = "Wind Speed (m/s)", title = "Wind Vector Plot") +
  theme_minimal()

Result[Sonic=='SonicC' & Campaign == 'MK',{
	# browser()
	ggplot(.SD, aes(x = st, y = 0)) +
	geom_segment(aes(xend = st * U_sonic * cos((90 - WD) * pi / 180),
                   yend = U_sonic * sin((90 - WD) * pi / 180)),
               arrow = arrow(length = unit(0.15, "inches")), size = 1) +
	theme_minimal()
}]

plot(1,type='n',xlim=range(Result[Sonic=='SonicC' & Campaign == 'MK',st]),ylim=c(-100,100))
plot(1,type='n',xlim=c(0,360),ylim=c(-100,100))
plot(1)
addWD(WD=12,U=4)

##################################
##################################
#####                        #####
#####    Pressure and MFC    #####
#####                        #####
##################################
##################################

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

ggsave(Fig_MFC_Press_paper,file=file.path(PfadFigures,"Fig_MFC_Press_paper.png"),width=30,height=30/1.57, units="cm")



#################################
#################################
#####                       #####
#####    GasFinder plots    #####
#####                       #####
#################################
#################################


######################################
### uncorrected concentration data ###
######################################

# MFC_10min <- pool(MFC,granularity="10mins")
# MFC_1min <- pool(MFC,granularity="1mins")

# which(MFC_10min[QV2,"Q_MFC"] > 0)
# which(MFC_1min[QV2,"Q_MFC"] > 0)
# # start "03.09.2021 15:42" - 5min
# # end "03.09.2021 16:44" + 10min
# indQV2_MFC <- "03.09.2021 15:38 to 03.09.2021 16:54"
# which(MFC_10min[MK,"Q_MFC"] > 0)
# which(MFC_1min[MK,"Q_MFC"] > 0)
# # start_1 "18.09.2021 12:57" - 5min
# # end_1 "18.09.2021 13:50" + 10min
# indMK_MFC_1 <- "18.09.2021 12:52 to 18.09.2021 14:00"
# # start_2 "19.09.2021 10:30" - 5min
# # end_2 "19.09.2021 16:48" + 10min
# indMK_MFC_2 <- "19.09.2021 10:25 to 19.09.2021 16:58"
# # start_3 "19.09.2021 21:52" - 5min
# # end_3 "20.09.2021 06:52" + 10min
# indMK_MFC_3 <- "19.09.2021 21:47 to 20.09.2021 07:02"
# indMK_MFC <- cbind(indMK_MFC_1,indMK_MFC_2,indMK_MFC_3)
# which(MFC_10min[QV3,"Q_MFC"] > 0)
# which(MFC_1min[QV3,"Q_MFC"] > 0)
# # start "22.09.2021 09:45" - 5min
# # end "22.09.2021 14:20" + 10min
# indQV3_MFC <- "22.09.2021 09:40 to 22.09.2021 14:30"

# ind_MFC <- cbind(indQV2_MFC,indMK_MFC,indQV3_MFC) ## Zeiten, bei denen Quelle lief (5min vor release bis 10min nach release)


GFs10min_uncorr2 <- merge(GFs10min_uncorr,MFC[,"Q_MFC"])
# GFs10min_uncorr2[which(GFs10min_uncorr2[,"Q_MFC"] > 0) & GFs10min_uncorr2[,"MK"] == "MC") ]


dt_GF_uncorr <- cbind(as.data.table(GFs10min_uncorr2["05.03.2021 16:00 to "]),st=st(GFs10min_uncorr2["05.03.2021 16:00 to "]))
dt_GF_uncorr$MK <- factor(dt_GF_uncorr$MK, levels=c("QV1","QV2","MK","QV3"),labels=c("IC1","IC2","MC","IC3"))
# limits_Conc <- parse_date_time3("")
mbreaks_GF <- parse_date_time3("06.03.2021",tz="Etc/GMT-1") + c(0:20) * 86400

Conc_uncorr <-	dt_GF_uncorr[MK != "IC1",{
	ggplot(.SD, aes(x=st)) +
	geom_line(aes(y=GF16,col="GF16")) +
	geom_line(aes(y=GF17,col="GF17")) +
	geom_line(aes(y=GF18,col="GF18")) +
	geom_line(aes(y=GF25,col="GF25")) +
	geom_line(aes(y=GF26,col="GF26")) +
	ylab(expression("CH"[4]*" conc. [mg m"^-3*"]")) +
	# ylab(expression("CH"[4]*" Concentration [mg m"^-3*"]")) +
	xlab(NULL) +
	scale_x_datetime(breaks=mbreaks_GF[c(1,3,5,14,16,17,19,21)],date_labels = "%b %d",expand=c(0.02,0.02)) +
	# scale_x_datetime(breaks="2 days",date_labels = "%b %d",expand=c(0.02,0.02)) +
	scale_colour_manual(name=NULL,values=CH4Cols) + 
	facet_grid(~ MK, scales="free_x", space="free_x") +
	theme_classic(base_size=18)
	# theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.margin=unit(c(0,0,0,0),units="lines")
	# 	,legend.position = c(0.043,0.78), text = element_text(size=30))
}]


dC26 <- dt_GF_uncorr[(MK == "IC2" & is.na(Q_MFC) | Q_MFC == 0) |
		(MK == "MC" & is.na(Q_MFC) | Q_MFC == 0) | 
		(MK == "IC3"), {
		# (MK == "IC3" & is.na(Q_MFC) | Q_MFC == 0), {
# dC26 <- dt_GF_uncorr[MK != "IC1" & is.na(Q_MFC), {
	ggplot(.SD, aes(x=st)) +
	geom_hline(yintercept=0,col="black") +
	geom_point(aes(y=GF16-GF26),alpha=0.5,size=0.8, col=CH4Cols["GF16"]) +
	geom_smooth(aes(y=GF16-GF26,col="GF16"),se=FALSE,size=0.5) +
	geom_point(aes(y=GF17-GF26),alpha=0.5,size=0.8, col=CH4Cols["GF17"]) +
	geom_smooth(aes(y=GF17-GF26,col="GF17"),se=FALSE,size=0.5) +
	geom_point(aes(y=GF18-GF26),alpha=0.5,size=0.8, col=CH4Cols["GF18"]) +
	geom_smooth(aes(y=GF18-GF26,col="GF18"),se=FALSE,size=0.5) +
	geom_point(aes(y=GF25-GF26),alpha=0.5,size=0.8, col=CH4Cols["GF25"]) +
	geom_smooth(aes(y=GF25-GF26,col="GF25"),se=FALSE,size=0.5) +
	ylab(expression(Delta*"CH"[4]*" conc. to GF26 [mg m"^-3*"]")) +
	xlab(NULL) +
	scale_x_datetime(breaks=mbreaks_GF[c(1,3,5,14,16,17,19,21)],date_labels = "%b %d",expand=c(0.02,0.02)) +
	# scale_x_datetime(breaks="2 days",date_labels = "%b %d",expand=c(0.02,0.02)) +
	scale_colour_manual(name=NULL,values=CH4Cols[1:4]
		,labels=c(expression(Delta*"GF16"),expression(Delta*"GF17"),expression(Delta*"GF18"),expression(Delta*"GF25"))) + 
	facet_grid(~ MK, scales="free_x", space="free_x") +
	theme_classic(base_size = 18)
	# theme_bw(base_size = 20) +
 #    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none", text = element_text(size=30))
}]

GF_uncorr <- ggarrange(Conc_uncorr,dC26,ncol=1,align="hv")
ggsave(file.path(PfadFigures,"GF_uncorr.png"),GF_uncorr,width=30, height=30/1.7,units="cm")



########################
### Korrigierte Konz ###
########################

dt_GF_corr1 <- cbind(as.data.table(GFs10min_corr1["05.03.2021 16:00 to "]),st=st(GFs10min_corr1["05.03.2021 16:00 to "]))
dt_GF_corr1$MK <- factor(dt_GF_corr1$MK, levels=c("QV1","QV2","MK","QV3"),labels=c("IC1","IC2","MC","IC3"))

GFs10min$MK <- NA_real_
GFs10min[QV1,"MK"] <- "QV1"
GFs10min[QV2,"MK"] <- "QV2"
GFs10min[MK,"MK"] <- "MK"
GFs10min[QV3,"MK"] <- "QV3"
dt_GF <- cbind(as.data.table(GFs10min["05.03.2021 16:00 to "]),st=st(GFs10min["05.03.2021 16:00 to "]))
dt_GF$MK <- factor(dt_GF$MK, levels=c("QV1","QV2","MK","QV3"),labels=c("IC1","IC2","MC","IC3"))
mbreaks_c <- parse_date_time3("18.03.2021 21:00",tz="Etc/GMT-1") + c(0:5) * 3600*12
limits_c <- c(parse_date_time3("18.03.2021 15:00",tz="Etc/GMT-1"),parse_date_time3("21.03.2021 11:00",tz="Etc/GMT-1"))

Conc_corr1 <-	dt_GF_corr1[MK != "IC1",{
	ggplot(.SD, aes(x=st)) +
	geom_line(aes(y=GF16,col="GF16")) +
	geom_line(aes(y=GF17,col="GF17")) +
	geom_line(aes(y=GF18,col="GF18")) +
	geom_line(aes(y=GF25,col="GF25")) +
	geom_line(aes(y=GF26,col="GF26")) +
	ylab(expression("CH"[4]*" conc. [mg m"^-3*"]")) +
	xlab(NULL) +
	scale_y_continuous(breaks=seq(1.3,2.3,by=0.2),limits=c(1.23,2.3)) +
	scale_x_datetime(breaks="2 days",date_labels = "%b %d",expand=c(0.02,0.02)) +
	scale_colour_manual(name=NULL,values=CH4Cols) + 
	facet_grid(~ MK, scales="free_x", space="free_x") +
	theme_classic(base_size=18)
 #    theme(legend.position = "none", text = element_text(size=30))
}]

Conc <-	dt_GF[MK == "MC",{
# Conc <-	dt_GF[MK != "IC1",{
	ggplot(.SD, aes(x=st)) +
	geom_line(aes(y=GF16,col="GF16")) +
	geom_line(aes(y=GF17,col="GF17")) +
	geom_line(aes(y=GF18,col="GF18")) +
	geom_line(aes(y=GF25,col="GF25")) +
	geom_line(aes(y=GF26,col="GF26")) +
	ylab(expression("CH"[4]*" conc. [mg m"^-3*"]")) +
	xlab(NULL) +
	scale_y_continuous(breaks=seq(1.3,2.3,by=0.2),limits=c(1.23,2.3)) +
	scale_x_datetime(breaks=mbreaks_c,limits=limits_c,date_labels = "%H:%m",expand=c(0.02,0.02)) +
	scale_colour_manual(name=NULL,values=CH4Cols) + 
	facet_grid(~ MK, scales="free_x", space="free_x") +
	theme_classic(base_size = 18)
 #    theme(legend.position = c(0.043,0.78), text = element_text(size=30))
}]


GF_corr <- ggarrange(Conc_corr1,Conc,ncol=1,align="hv")
ggsave(file.path(PfadFigures,"GF_corr.png"),GF_corr,width=30, height=30/1.7, units="cm")


limits_REL <- c(parse_date_time3('2021-03-19 09:00',tz='Etc/GMT-1'), parse_date_time3('2021-03-20 09:00',tz='Etc/GMT-1'))
Conc_REL <-	dt_GF[,{
# Conc <-	dt_GF[MK != "IC1",{
	ggplot(.SD, aes(x=st)) +
	geom_line(aes(y=GF16,col="GF16")) +
	geom_line(aes(y=GF17,col="GF17")) +
	geom_line(aes(y=GF18,col="GF18")) +
	geom_line(aes(y=GF25,col="GF25")) +
	geom_line(aes(y=GF26,col="GF26")) +
	ylab(expression("CH"[4]*" conc. [mg m"^-3*"]")) +
	xlab(NULL) +
	scale_y_continuous(breaks=seq(1.3,2.3,by=0.1),limits=c(1.3,2.3)) +
	scale_x_datetime(limits=limits_REL,date_labels = "%H:%m",expand=c(0.02,0.02)) +
	# scale_x_datetime(breaks=mbreaks_c,limits=limits_REL,date_labels = "%H:%m",expand=c(0.02,0.02)) +
	scale_colour_manual(name=NULL,values=CH4Cols) + 
	# facet_grid(~ MK, scales="free_x", space="free_x") +
	theme_bw(base_size = 18)
 #    theme(legend.position = c(0.043,0.78), text = element_text(size=30))
}]


dConc_REL <- dt_GF[,{
	ggplot(.SD, aes(x=st)) +
	geom_line(aes(y=GF16-GF26,col="GF16")) +
	geom_line(aes(y=GF17-GF26,col="GF17")) +
	geom_line(aes(y=GF18-GF26,col="GF18")) +
	geom_line(aes(y=GF25-GF26,col="GF25")) +
	ylab(expression("delta CH"[4]*" conc. [mg m"^-3*"]")) +
	xlab(NULL) +
	scale_y_continuous(breaks=seq(-0.1,1,by=0.1),limits=c(-0.05,0.95)) +
	scale_x_datetime(breaks=mbreaks_c,limits=limits_REL,date_labels = "%H:%m",expand=c(0.02,0.02)) +
	scale_colour_manual(name=NULL,values=CH4Cols) + 
	theme_bw(base_size = 18)
}]


ggsave(file.path(PfadFigures,"GF_Release.png"),Conc_REL,width=30, height=30/1.7, units="cm")
ggsave(file.path(PfadFigures,"deltaGF_Release.png"),dConc_REL,width=30, height=30/1.7, units="cm")




################################
################################
#####                      #####
#####    Emissionsplots    #####
#####                      #####
################################
################################

# Result$Sonic_ord <- factor(Result$Sonic, levels=c("SonicA","Sonic2","SonicB","SonicC"),labels=c("SonicA","SonicD","SonicB","SonicC"))

# Result_Q <- data.table(rbind(
#  	cbind(melt(Result[Sonic_ord == "SonicA"], id.vars="st",measure.vars=c("Q_GF16","Q_GF17","Q_GF18","Q_GF25","Q_GF26"), variable.name="Q_GF", value.name="Emissions"),Sonic="SonicA")
#  	,cbind(melt(Result[Sonic_ord == "SonicB"], id.vars="st",measure.vars=c("Q_GF16","Q_GF17","Q_GF18","Q_GF25","Q_GF26"), variable.name="Q_GF", value.name="Emissions"),Sonic="SonicB")
#  	,cbind(melt(Result[Sonic_ord == "SonicC"], id.vars="st",measure.vars=c("Q_GF16","Q_GF17","Q_GF18","Q_GF25","Q_GF26"), variable.name="Q_GF", value.name="Emissions"),Sonic="SonicC")
#  	,cbind(melt(Result[Sonic_ord == "SonicD"], id.vars="st",measure.vars=c("Q_GF16","Q_GF17","Q_GF18","Q_GF25","Q_GF26"), variable.name="Q_GF", value.name="Emissions"),Sonic="SonicD")
# 	))
# Result_Q$GF <- factor(Result_Q$Q_GF, levels=c("Q_GF17","Q_GF18","Q_GF16","Q_GF25"),labels=c("GF17","GF18","GF16","GF25"))

# MFC_10min <- pool(MFC,granularity="10mins",st.to="19.03.2021 09:00",et.to="20.03.2021 09:00")
# MFC_dt <- data.table(cbind(MFC_10min,st=st(MFC_10min)))

# # Result[,Q_MFC2 := Q_MFC]
# # Result[Campaign == "MK" & is.na(Q_MFC), Q_MFC2 := 0]

# Fig_sameSonic <- Result[st > "2021-03-19 09:00" & st < "2021-03-20 09:00" & !is.na(Sonic_ord),{
# 	ggplot(.SD, aes(x=st)) +
# 	geom_hline(yintercept=0,linetype=3) +
# 	geom_line(aes(y=Q_GF16,colour="GF16"),size=0.6) +
# 	geom_line(aes(y=Q_GF17,colour="GF17"),size=0.6) +
# 	geom_line(aes(y=Q_GF18,colour="GF18"),size=0.6) +
# 	geom_line(aes(y=Q_GF25,colour="GF25"),size=0.6) +
# 	geom_point(aes(y=Q_GF16,colour="GF16"),size=1) +
# 	geom_point(aes(y=Q_GF17,colour="GF17"),size=1) +
# 	geom_point(aes(y=Q_GF18,colour="GF18"),size=1) +
# 	geom_point(aes(y=Q_GF25,colour="GF25"),size=1) +
# 	# geom_line(aes(y=Q_MFCex,colour="MFC"),size=1) +
# 	geom_line(data=MFC_dt[st > "2021-03-19 09:00" & st < "2021-03-20 09:00"],aes(x=st,y=Q_MFCex,colour="MFC"),size=1) +
# 	xlab(NULL) +
# 	ylab(expression("CH"[4]*" emissions [kg h"^-1*"]")) +
# 	facet_grid(Sonic_ord~.) +
# 	scale_colour_manual(name=NULL,values=c(CH4Cols[c(2,3,1,4)],MFC="black")) +
# 	theme_classic(base_size = 18)
#     # theme(legend.position = "right", text = element_text(size=30))
# }]

# ggsave(file.path(PfadFigures,"Emission_sameGF_diffSonic.png"),Fig_sameSonic,width = 30, height = 30/1.1, units="cm")


# Fig_sameGF <- Result_Q[st > "2021-03-19 09:00" & st < "2021-03-20 09:00" & GF %in% c("GF16","GF17","GF18","GF25"),{
# 	ggplot(.SD, aes(x=st,y=Emissions,colour=Sonic)) +
# 	geom_hline(yintercept=0,linetype=3) +
# 	geom_line(size=0.6) +
# 	geom_point(size=1) +
# 	geom_line(data=MFC_dt[st > "2021-03-19 09:00" & st < "2021-03-20 09:00"],aes(x=st,y=Q_MFCex,colour="MFC"),size=1) +
# 	xlab(NULL) +
# 	ylab(expression("CH"[4]*" emissions [kg h"^-1*"]")) +
# 	facet_grid(GF~.) +
# 	scale_colour_manual(name=NULL,values=c(SonicCols[1:4],MFC="black")) +
# 	theme_classic(base_size = 18)
#     # theme(legend.position = "none", text = element_text(size=30))
# }]

# ggsave(file.path(PfadFigures,"Emission_sameSonic_diffGF.png"),Fig_sameGF,width = 30, height = 30/1.1, units="cm")



#### QV3

# Fig_sameSonic_QV3 <- Result[st > "2021-03-22 07:00" & st < "2021-03-22 16:00" & Sonic_ord == "SonicD",{
# 	# browser()
# 	ggplot(.SD, aes(x=st)) +
# 	geom_hline(yintercept=0,linetype=3) +
# 	geom_line(aes(y=Q_GF16,colour="GF16"),size=0.6) +
# 	geom_line(aes(y=Q_GF17,colour="GF17"),size=0.6) +
# 	geom_line(aes(y=Q_GF18,colour="GF18"),size=0.6) +
# 	geom_line(aes(y=Q_GF25,colour="GF25"),size=0.6) +
# 	geom_line(aes(y=Q_GF26,colour="GF26"),size=0.6) +
# 	geom_point(aes(y=Q_GF16,colour="GF16"),size=1) +
# 	geom_point(aes(y=Q_GF17,colour="GF17"),size=1) +
# 	geom_point(aes(y=Q_GF18,colour="GF18"),size=1) +
# 	geom_point(aes(y=Q_GF25,colour="GF25"),size=1) +
# 	geom_point(aes(y=Q_GF26,colour="GF26"),size=1) +
# 	geom_line(aes(y=Q_MFCex,colour="MFC"),size=1) +
# 	xlab(NULL) +
# 	ylab(expression("CH"[4]*" emissions [kg h"^-1*"]")) +
# 	scale_x_datetime(breaks="2 hours",date_labels="%H:%M") +
# 	# facet_grid(Sonic_ord~.) +
# 	scale_colour_manual(name=NULL,values=c(CH4Cols,MFC="black")) +
# 	theme_classic(base_size = 18)
#     # theme(legend.position = "right", text = element_text(size=30))
# }]

# ggsave(file.path(PfadFigures,"Emission_sameGF_diffSonic_QV3.png"),Fig_sameSonic_QV3,width = 30, height = 30/3, units="cm")

# Result[Campaign == "MK" & (is.na(Q_MFC) | Q_MFC ==  0), summary(Q_GF16,Q_GF17,Q_GF18,Q_GF25)]
# Result[Campaign == "MK" & (is.na(Q_MFC) | Q_MFC ==  0), sd(c(Q_GF16,Q_GF17,Q_GF18,Q_GF25),na.rm=TRUE)]

# Result[Campaign=="MK",plot(Q_MFC,type="l")]
# Result[Campaign=="MK" & Sonic == "SonicA"][1:150,plot(Q_MFCex,type="l")]
# Result[Campaign=="MK" & Sonic == "SonicB"][80:230,plot(Q_MFCex,type="l")]
# Result[Campaign=="MK" & Sonic == "SonicC"][100:250,plot(Q_MFCex,type="l")]
# Result[Campaign=="MK" & Sonic == "Sonic2"][80:230,plot(Q_MFCex,type="l")]

# Result[st > "2021-03-19 10:00" & st < "2021-03-20 08:00" & Sonic == "Sonic2",.(st,Q_MFC,Q_MFCex,WD_WS2,Ustar)]
# Result[st > "2021-03-19 10:00" & st < "2021-03-19 12:00" & Sonic == "SonicB",.(st,Q_MFC,Q_MFCex,WD_WS2,Ustar)]


#############
### Paper ###
#############

#########################################################
### plot recovery rates instead of absoulte emissions ###
#########################################################

Result$Sonic_ord <- factor(Result$Sonic, levels=c("SonicC","SonicA","Sonic2","SonicB"),labels=c("3DUA-UW","3DUA-50m","3DUA-100m","3DUA-150m"))
Result$Sonic_ord2 <- factor(Result$Sonic, levels=c("SonicC","SonicA","Sonic2","SonicB"),labels=c("UA-UW","UA-2.0h","UA-5.3h","UA-8.6h"))

Result[Q_MFCex > 4,R_GF16 := Q_GF16 / Q_MFCex]
Result[Q_MFCex > 4,R_GF17 := Q_GF17 / Q_MFCex]
Result[Q_MFCex > 4,R_GF18 := Q_GF18 / Q_MFCex]
Result[Q_MFCex > 4,R_GF25 := Q_GF25 / Q_MFCex]
Result[Q_MFCex > 4,R_GF26 := Q_GF26 / Q_MFCex]


Result_R <- data.table(rbind(
 	cbind(melt(Result[Sonic_ord2 == "UA-UW"], id.vars=c('st','L','U_sonic','Ustar','WD'),measure.vars=c("R_GF16","R_GF17","R_GF18","R_GF25","R_GF26"), variable.name="R_GF", value.name="Recovery"),Sonic="UA-UW")
 	,cbind(melt(Result[Sonic_ord2 == "UA-50m"], id.vars=c('st','L','U_sonic','Ustar','WD'),measure.vars=c("R_GF16","R_GF17","R_GF18","R_GF25","R_GF26"), variable.name="R_GF", value.name="Recovery"),Sonic="UA-50m")
 	,cbind(melt(Result[Sonic_ord2 == "UA-100m"], id.vars=c('st','L','U_sonic','Ustar','WD'),measure.vars=c("R_GF16","R_GF17","R_GF18","R_GF25","R_GF26"), variable.name="R_GF", value.name="Recovery"),Sonic="UA-100m")
 	,cbind(melt(Result[Sonic_ord2 == "UA-150m"], id.vars=c('st','L','U_sonic','Ustar','WD'),measure.vars=c("R_GF16","R_GF17","R_GF18","R_GF25","R_GF26"), variable.name="R_GF", value.name="Recovery"),Sonic="UA-150m")
	))
Result_R$GF <- factor(Result_R$R_GF, levels=c("R_GF17","R_GF18","R_GF16","R_GF25"),labels=c("GF-50m","GF-100m","GF-150m","GF-200m"))
Result_R$OP <- factor(Result_R$R_GF, levels=c("R_GF17","R_GF18","R_GF16","R_GF25"),labels=c("OP-50m","OP-100m","OP-150m","OP-200m"))
Result_R$OPh <- factor(Result_R$R_GF, levels=c("R_GF17","R_GF18","R_GF16","R_GF25"),labels=c("OP-2.0h","OP-5.3h","OP-8.6h","OP-12h"))


Fig_recovery_sameSonic <- Result[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & !is.na(Sonic_ord),{
	ggplot(.SD, aes(x=st)) +
	geom_hline(yintercept=1,linetype=3) +
	geom_line(aes(y=R_GF16,colour="GF-150m"),size=0.6) +
	geom_line(aes(y=R_GF17,colour="GF-50m"),size=0.6) +
	geom_line(aes(y=R_GF18,colour="GF-100m"),size=0.6) +
	geom_line(aes(y=R_GF25,colour="GF-200m"),size=0.6) +
	geom_point(aes(y=R_GF16,colour="GF-150m"),size=1) +
	geom_point(aes(y=R_GF17,colour="GF-50m"),size=1) +
	geom_point(aes(y=R_GF18,colour="GF-100m"),size=1) +
	geom_point(aes(y=R_GF25,colour="GF-200m"),size=1) +
	xlab(NULL) +
	ylab('bLS recovery rate [%]') +
	# ylab(expression("CH"[4]*" recovery rate [%]")) +
	facet_grid(Sonic_ord~.) +
	scale_colour_manual(name=NULL,values=c(CH4Cols_GF[2:5])) +
	theme_bw(base_size = 18) + theme(strip.background =element_rect(fill="white"),panel.grid = element_blank())
    # theme(legend.position = "right", text = element_text(size=30))
}]

ggsave(file.path(PfadFigures,"Emission_sameGF_diffSonic_recovery_Paper.png"),Fig_recovery_sameSonic,width = 30, height = 30/1.1, units="cm")

# new names
Fig_recovery_disturbedSonic <- Result[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & !is.na(Sonic_ord2) & Sonic_ord2 != 'UA-UW',{
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
	ylab('bLS recovery rate [%]') +
	facet_grid(Sonic_ord2~.) +
	scale_colour_manual(name=NULL,values=c('OP-2.0h'='#F8766D','OP-5.3h'='#7CAE00','OP-8.6h'='#00BFC4','OP-12h'='#C77CFF'),breaks=c('OP-2.0h','OP-5.3h','OP-8.6h','OP-12h')) +
	theme_bw(base_size = 18) +
	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3), strip.background =element_rect(fill="white"),
	legend.title=element_blank(), legend.position='top')
    # theme(legend.position = "right", text = element_text(size=30))
}]

ggsave(file.path(PfadFigures,"Recovery_distubedSonics_Paper.png"),Fig_recovery_disturbedSonic,width = 30, height = 30/1.5, units="cm")

## same as above but different ylabel
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

ggsave(file.path(PfadFigures,"Recovery_distubedSonics_Paper_IDM.png"),Fig_recovery_disturbedSonic_IDM,width = 30, height = 30/1.5, units="cm")


# new names
limits_recovery <- c(parse_date_time3("2021-03-19 08:00",tz='Etc/GMT-1'), parse_date_time3("2021-03-20 09:00",tz='Etc/GMT-1'))

Fig_recovery_UWSonic <- Result[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & Sonic_ord2=='UA-UW',{
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
	scale_x_datetime(limits=limits_recovery) +
	ylab('bLS recovery rate [%]') +
	scale_colour_manual(name=NULL,values=c('OP-2.0h'='#F8766D','OP-5.3h'='#7CAE00','OP-8.6h'='#00BFC4','OP-12h'='#C77CFF'),breaks=c('OP-2.0h','OP-5.3h','OP-8.6h','OP-12h')) +
	theme_bw(base_size = 18) +
	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3),
	strip.background =element_rect(fill="white"), legend.title=element_blank(), legend.position='top',
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),plot.margin = unit(c(0, 0, 0, 0), "cm"))
}]

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
	scale_x_datetime(limits=limits_recovery) +
	ylab('IDM recovery rate [%]') +
	scale_colour_manual(name=NULL,values=c('OP-2.0h'='#F8766D','OP-5.3h'='#7CAE00','OP-8.6h'='#00BFC4','OP-12h'='#C77CFF'),breaks=c('OP-2.0h','OP-5.3h','OP-8.6h','OP-12h')) +
	theme_bw(base_size = 18) +
	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3),
	strip.background =element_rect(fill="white"), legend.title=element_blank(), legend.position='top',
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),plot.margin = unit(c(0, 0, 0, 0), "cm"))
}]

Result[!is.na(R_GF16) | !is.na(R_GF17) | !is.na(R_GF18) | !is.na(R_GF25),L_help := L]

Fig_stability_UWSonic <- Result[st > "2021-03-19 08:00" & st < "2021-03-20 09:00" & Sonic_ord2=='UA-UW',{
	ggplot(.SD, aes(x=st)) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 11:19",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 11:40",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 01:08",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("20.03.2021 01:21",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1"),ymin=-Inf,ymax=Inf),fill='grey95',alpha=0.15) +
	geom_hline(yintercept=0,linetype=3) +
	geom_line(aes(y=1/L),size=0.6) +
	geom_point(aes(y=1/L),size=1) +
	# geom_line(aes(y=1/L,colour='UA-UW'),size=0.6) +
	# geom_point(aes(y=1/L,colour='UA-UW'),size=1) +
	scale_y_continuous(breaks=seq(-0.2,0.4,by=0.05),limits=c(-0.4,0.4),expand=c(-0.108,-0.208)) +
	# scale_y_continuous(limits=c(-0.1,0.1)) +
	xlab(NULL) +
	scale_x_datetime(limits=limits_recovery) +
	ylab(expression('L'^-1*' [m'^-1*']')) +
	scale_colour_manual(name=NULL,values=c('UA-UW'='black')) +
	theme_bw(base_size = 18) +
	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3),
	legend.title=element_blank(),plot.margin = unit(c(0, 0, 0, 0), "null"))
}]


Fig_recovery_stability_UWSonic <- ggarrange(Fig_recovery_UWSonic,Fig_stability_UWSonic,align='v',ncol=1,heights=c(5,2))
Fig_recovery_stability_UWSonic_IDM <- ggarrange(Fig_recovery_UWSonic_IDM,Fig_stability_UWSonic,align='v',ncol=1,heights=c(5,2))

# ggarr_built <- ggplot_build(Fig_recovery_stability_UWSonic)

# # Access the coordinates and dimensions of the plot
# print(ggarr_built$layout)

# # Fig_rec_stab_UWSonic <- 
# ggplot_build(Fig_recovery_stability_UWSonic)
# Fig_recovery_stability_UWSonic + geom_segment(x=seq(-100,100,by=0.1),xend=seq(-100,100,by=0.1),y=seq(-100,100,by=0.1),inherit.aes = FALSE,linewidth=3)

# 	annotate("text", x = 0.287, y = 0.95, label = "|-------------------------------      Release 1      -------------------------------|", color = "black") +
# 	annotate("text", x = 0.727, y = 0.95, label = "|------------------------------------------------      Release 2      ----------------------------------------------------|", color = "black")


ggsave(file.path(PfadFigures,"recovery_stability_Paper.png"),Fig_recovery_stability_UWSonic,width = 30, height = 30/2, units="cm")
ggsave(file.path(PfadFigures,"recovery_stability_Paper_IDM.png"),Fig_recovery_stability_UWSonic_IDM,width = 30, height = 30/2, units="cm")



Result_R[L > 0,stability := 'stable']
Result_R[L < 0,stability := 'unstable']
Result_R[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & Sonic=='UA-UW' & !is.na(stability) & !is.na(OP),
.(mean=round(mean(Recovery,na.rm=TRUE),2),median=round(median(Recovery,na.rm=TRUE),2),sd=round(sd(Recovery,na.rm=TRUE),2)),by=.(stability,OP)][order(stability,OP)]

Result_R[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & Sonic=='UA-UW' & !is.na(stability) & !is.na(OP),
.(mean=round(mean(Recovery,na.rm=TRUE),2),median=round(median(Recovery,na.rm=TRUE),2),sd=round(sd(Recovery,na.rm=TRUE),2)),by=.(stability)][order(stability)]

Result_R[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & Sonic=='UA-UW' & !is.na(stability) & !is.na(OP),
.(mean=round(mean(Recovery,na.rm=TRUE),2),median=round(median(Recovery,na.rm=TRUE),2),sd=round(sd(Recovery,na.rm=TRUE),2)),by=.(OP)][order(OP)]

## recovery
Result_R[st > "2021-03-19 08:30" & st < "2021-03-20 08:00" & !is.na(stability) & !is.na(OP),
.(mean=round(mean(Recovery,na.rm=TRUE),2),median=round(median(Recovery,na.rm=TRUE),2),sd=round(sd(Recovery,na.rm=TRUE),2)),by=.(stability,OP)][order(stability,OP)]






names(DUACols) <- c('UA-UW','UA-50m','UA-100m','UA-150m')

Fig_recovery_sameGF <- Result_R[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & OP %in% c("OP-50m","OP-100m","OP-150m","OP-200m"),{
	ggplot(.SD, aes(x=st,y=Emissions,colour=Sonic_ord)) +
	geom_hline(yintercept=1,linetype=3) +
	geom_line(size=0.6) +
	geom_point(size=1) +
	xlab(NULL) +
	ylab('IDM recovery rate [%]') +
	# ylab(expression("CH"[4]*" recovery rate [%]")) +
	facet_grid(OP~.) +
	scale_colour_manual(name=NULL,values=DUACols) +
	theme_bw(base_size = 18) +
	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3), strip.background =element_rect(fill="white"), legend.title=element_blank())
}]

ggsave(file.path(PfadFigures,"Emission_sameSonic_diffGF_recovery_Paper.png"),Fig_recovery_sameGF,width = 30, height = 30/1.1, units="cm")

##### for IC2 rep. QV3:
CH4Cols_OP <- CH4Cols_GF
names(CH4Cols_OP) <- paste('OP',c('UW','50m','100m','150m','200m'),sep='-')

Fig_recovery_IC2 <- Result[st > "2021-03-22 09:30" & st < "2021-03-22 14:15" & Sonic == 'Sonic2',{
	ggplot(.SD, aes(x=st)) +
	geom_hline(yintercept=1,linetype=3) +
	geom_line(aes(y=R_GF16,colour="OP-150m"),size=0.6) +
	geom_line(aes(y=R_GF17,colour="OP-50m"),size=0.6) +
	geom_line(aes(y=R_GF18,colour="OP-100m"),size=0.6) +
	geom_line(aes(y=R_GF25,colour="OP-200m"),size=0.6) +
	geom_line(aes(y=R_GF26,colour="OP-UW"),size=0.6) +
	geom_point(aes(y=R_GF16,colour="OP-150m"),size=1) +
	geom_point(aes(y=R_GF17,colour="OP-50m"),size=1) +
	geom_point(aes(y=R_GF18,colour="OP-100m"),size=1) +
	geom_point(aes(y=R_GF25,colour="OP-200m"),size=1) +
	geom_point(aes(y=R_GF26,colour="OP-UW"),size=1) +
	xlab(NULL) +
	ylim(0,1) +
	ylab('IDM recovery rate [%]') +
	# ylab(expression("CH"[4]*" recovery rate [%]")) 
	scale_colour_manual(name=NULL,values=CH4Cols_OP,breaks=c('OP-UW','OP-50m','OP-100m','OP-150m','OP-200m')) +
	# scale_colour_manual(name=NULL,values=CH4Cols_OP) +
	theme_bw(base_size = 18) +
	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3), strip.background =element_rect(fill="white"), legend.title=element_blank())
    # theme(legend.position = "right", text = element_text(size=30))
}]

limits_recovery_IC2 <- c(parse_date_time3("2021-03-22 09:00",tz='Etc/GMT-1'), parse_date_time3("2021-03-22 15:00",tz='Etc/GMT-1'))

Result[Campaign == 'QV3' & Q_MFCex > 0,]

Fig_recovery_IC2 <- Result[st > "2021-03-22 09:30" & st < "2021-03-22 14:15" & Sonic_ord2=='UA-2.0h',{
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
	ylab('bLS recovery rate [%]') +
	scale_colour_manual(name=NULL,values=c('OP-2.0h'='#F8766D','OP-5.3h'='#7CAE00','OP-8.6h'='#00BFC4','OP-12h'='#C77CFF','OP-UW'='black'),breaks=c('OP-2.0h','OP-5.3h','OP-8.6h','OP-12h','OP-UW')) +
	theme_bw(base_size = 18) +
	theme(panel.grid.major = element_line(size = 0.6),panel.grid.minor = element_line(size = 0.3),
	strip.background =element_rect(fill="white"), legend.title=element_blank(), legend.position='top')
}]

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


ggsave(file.path(PfadFigures,"Recovery_IC2_Paper.png"),Fig_recovery_IC2,width = 30, height = 30/2, units="cm")
ggsave(file.path(PfadFigures,"Recovery_IC2_Paper_IDM.png"),Fig_recovery_IC2_IDM,width = 30, height = 30/2, units="cm")





##############################
##############################
#####                    #####
#####    WD variation    #####
#####                    #####
##############################
##############################


# sub_Result_var <- Result_var[st > "2021-03-19 09:00" & st < "2021-03-20 09:00" & Conc_corr == "P23" 
				# & (Sonic %in% unique(grep("SonicB_m",Sonic,value=TRUE)) | Sonic %in% unique(grep("SonicB_p",Sonic,value=TRUE))),]

sub_Result_var <- Result_var[st > "2021-03-19 09:00" & st < "2021-03-20 09:00" & Conc_corr == "P23" 
				& Sonic %in% unique(grep("SonicB_",Sonic,value=TRUE)),]

Result_Qvar <- data.table(
 	melt(sub_Result_var, id.vars=c("st","Sonic"),measure.vars=c("Q_GF17","Q_GF18","Q_GF16","Q_GF25"), variable.name="Q_GF", value.name="Emissions")
	)

Result_Qvar$GF <- factor(Result_Qvar$Q_GF, levels=c("Q_GF17","Q_GF18","Q_GF16","Q_GF25"),labels=c("GF17","GF18","GF16","GF25"))

Fig_sameGF <- Result_Qvar[,{
	ggplot(.SD, aes(x=st,y=Emissions,colour=Sonic)) +
	geom_hline(yintercept=0,linetype=3) +
	geom_line(size=1) +
	geom_point(size=0.6) +
	geom_line(data=MFC_dt[st > "2021-03-19 09:00" & st < "2021-03-20 09:00"],aes(x=st,y=Q_MFCex,colour="MFC"),size=1) +
	xlab(NULL) +
	scale_y_continuous(breaks=c(0,3,6,9),limits=c(-0.4662566,9)) +
	ylab(expression("CH"[4]*" emissions [kg h"^-1*"]")) +
	facet_grid(GF~.) +
	scale_colour_manual(name=NULL,values=c(WDvarCols,MFC="black")) +
	theme_classic(base_size = 18) +
    # labs(c)
    theme(legend.position="none")
}]

# ggsave(file.path(PfadFigures,"WDvar_GFs_SonicB.png"),Fig_sameGF,width = 30, height = 30/1.1, units="cm")
ggsave(file.path(PfadFigures,"WDvar_GFs_SonicB.png"),Fig_sameGF,width = 30, height = 30/0.97, units="cm")
# ggsave(file.path(PfadFigures,"WDvar_GFs_SonicB.pdf"),Fig_sameGF,width = 30, height = 30/1.1, units="cm")


#############
### Paper ###
#############

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

Result_Qvar$GF <- factor(Result_Qvar$Q_GF, levels=c("R_GF17","R_GF18","R_GF16","R_GF25"),labels=c("GF-50m","GF-100m","GF-150m","GF-200m"))
Result_Qvar$OP <- factor(Result_Qvar$Q_GF, levels=c("R_GF17","R_GF18","R_GF16","R_GF25"),labels=c("OP-50m","OP-100m","OP-150m","OP-200m"))


Fig_recovery_WDvar <- Result_Qvar[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & GF %in% c("GF-50m","GF-100m","GF-150m","GF-200m"),{
	ggplot(.SD, aes(x=st,y=Recovery,colour=Sonic)) +
	geom_hline(yintercept=1,linetype=3) +
	geom_line(linewidth=0.6) +
	geom_point(size=1) +
	xlab(NULL) +
	# scale_y_continuous(breaks=c(0,0.5,1,1.5),limits=c(0,1.75)) +
	ylab('bLS recovery rate [%]') +
	# ylab(expression("CH"[4]*" recovery rate [%]")) +
	facet_grid(GF~.) +
	scale_colour_manual(name=NULL,values=WDvarCols_10) +
	theme_bw(base_size = 18) + theme(legend.position='none',strip.background =element_rect(fill="white"),panel.grid = element_blank())
}]

ggsave(file.path(PfadFigures,"WDvar_recovery_Paper.png"),Fig_recovery_WDvar,width = 30, height = 30/0.97, units="cm")


# new names

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
	# scale_colour_viridis(name=NULL, discrete=TRUE) +
	theme_bw(base_size = 18) +
	theme(legend.position='none',strip.background =element_rect(fill="white"),panel.grid.major=element_line(size=0.6),panel.grid.minor=element_line(size=0.3))
}]

ggsave(file.path(PfadFigures,"WDvar_recovery_Paper.png"),Fig_recovery_WDvar10,width = 30, height = 30/0.97, units="cm")




# new colours

Viridis <- viridis_pal(alpha = 1, begin = 0, end = 1, option = "D")(20)
Viridis_WDvar <- c(Viridis[1:10],'black',Viridis[11:20])
names(Viridis_WDvar) <- paste0("SonicB",c(paste0("_m",1:10),"_0",paste0("_p",1:10)))
Result_Qvar[,Sonic_ord := factor(Sonic, levels = paste0("SonicB",c("_0",paste0("_m",1:10),paste0("_p",1:10))))]
Result_Qvar[, line_width := ifelse(Sonic == "SonicB_0", 1.2, 1)] # add line width column
Result_Qvar[, point_size := ifelse(Sonic == "SonicB_0", 1.5, 1.2)] # add line width column


Fig_recovery_WDvar10 <- Result_Qvar[st > "2021-03-19 10:30" & st < "2021-03-20 07:00" & GF %in% c("GF-50m","GF-100m","GF-150m","GF-200m")
								 & Sonic %in% names(Viridis_WDvar),{
	ggplot(.SD, aes(x=st,y=Recovery,colour=Sonic_ord)) +
	geom_hline(yintercept=1,linetype=3) +
	# geom_line(linewidth=0.6) +
	# geom_line(aes(linewidth=0.5)) +
	geom_line(aes(linewidth=line_width)) +
	scale_linewidth(range=c(0.6, 1.2)) +
	geom_point(aes(size=point_size)) +
	scale_size(range=c(0.6,1.5)) +
	xlab(NULL) +
	scale_y_continuous(breaks=c(0,0.5,1,1.5),limits=c(0,1.75)) +
	ylab('bLS recovery rate [%]') +
	facet_grid(OP ~ .) +
	scale_colour_manual(name=NULL,values=Viridis_WDvar) +
	# scale_colour_viridis(name=NULL, discrete=TRUE) +
	theme_bw(base_size = 18) + theme(legend.position='none',strip.background =element_rect(fill="white"),panel.grid = element_blank())
}]


ggsave(file.path(PfadFigures,"WDvar_recovery_Paper_Viridis.png"),Fig_recovery_WDvar10,width = 30, height = 30/0.97, units="cm")



##############################################
### Distances from the tree to the sensors ###
##############################################

##### calculate the middle of the GF path:
M_GF17 <- data.table(x=mean(Sensors_MK[Sensors_MK[,'Sensor Name'] == 'GF17',4]),y=mean(Sensors_MK[Sensors_MK[,'Sensor Name'] == 'GF17',5]))
M_GF18 <- data.table(x=mean(Sensors_MK[Sensors_MK[,'Sensor Name'] == 'GF18',4]),y=mean(Sensors_MK[Sensors_MK[,'Sensor Name'] == 'GF18',5]))
M_GF16 <- data.table(x=mean(Sensors_MK[Sensors_MK[,'Sensor Name'] == 'GF16',4]),y=mean(Sensors_MK[Sensors_MK[,'Sensor Name'] == 'GF16',5]))
M_GF25 <- data.table(x=mean(Sensors_MK[Sensors_MK[,'Sensor Name'] == 'GF25',4]),y=mean(Sensors_MK[Sensors_MK[,'Sensor Name'] == 'GF25',5]))

plot(Sonics_MK)
points(M_GF17,col='red')
points(M_GF18,col='red')
points(M_GF16,col='red')
points(M_GF25,col='red')

# difference between middle and sonics
Diff_GF17 <- sqrt((Sonics_MK[Sonics_MK[,'Sensor Name'] == 'SonicA',4]-M_GF17[,1])^2 +(Sonics_MK[Sonics_MK[,'Sensor Name'] == 'SonicA',5]-M_GF17[,2])^2)
Diff_GF18 <- sqrt((Sonics_MK[Sonics_MK[,'Sensor Name'] == 'Sonic2',4]-M_GF18[,1])^2 +(Sonics_MK[Sonics_MK[,'Sensor Name'] == 'Sonic2',5]-M_GF18[,2])^2)
Diff_GF16 <- sqrt((Sonics_MK[Sonics_MK[,'Sensor Name'] == 'SonicB',4]-M_GF16[,1])^2 +(Sonics_MK[Sonics_MK[,'Sensor Name'] == 'SonicB',5]-M_GF16[,2])^2)
# lets take the middle

# Distance to tree
D_GF17 <- round(sqrt((Baum[,4]-M_GF17[,1])^2 +(Baum[,5]-M_GF17[,2])^2),1)
D_GF18 <- round(sqrt((Baum[,4]-M_GF18[,1])^2 +(Baum[,5]-M_GF18[,2])^2),1)
D_GF16 <- round(sqrt((Baum[,4]-M_GF16[,1])^2 +(Baum[,5]-M_GF16[,2])^2),1)
D_GF25 <- round(sqrt((Baum[,4]-M_GF25[,1])^2 +(Baum[,5]-M_GF25[,2])^2),1)

# D/h
round(D_GF17/15,1)
round(D_GF18/15,1)
round(D_GF16/15,1)
round(D_GF25/15,1)


plot(Sources,Sonics_MK)


# converse to matrix
square_coords <- matrix(unlist(Sources[,2:3]), nrow=4, byrow=FALSE)

# Define the coordinates of the point
point_coords <-  c(583952.8, 210336.9) # I just copy&pate the coordinates

# Calculate the distances to the vertices
distances <- dist(rbind(point_coords, square_coords), diag=TRUE)

# Take the minimum distance
min_distance <- min(distances[1:4])
87.8 m




###################################
###################################
#####                         #####
#####    Stability classes    #####
#####                         #####
###################################
###################################


Sonic2[,"MK"] <-  NA_real_
Sonic2[MK,"MK"] <-  "MK"
Sonic2[QV3,"MK"] <-  "QV3"

SonicAll <- merge(merge(
        merge(
            merge(cbind(as.data.table(SonicB),start_interval=st(SonicB),end_interval = et(SonicB)), cbind(as.data.table(SonicA),start_interval=st(SonicA),end_interval = et(SonicA)), by = c('start_interval', 'end_interval'), suffixes = c('', '.a'))
            , cbind(as.data.table(SonicC),start_interval=st(SonicC),end_interval = et(SonicC)), by = c('start_interval', 'end_interval'), suffixes = c('', '.c'))
                , cbind(as.data.table(Sonic2),start_interval=st(Sonic2),end_interval = et(Sonic2)), by = c('start_interval', 'end_interval'), suffixes = c('.b', '.2'))
                    , cbind(as.data.table(WS_10min),start_interval=st(WS_10min),end_interval=et(WS_10min)), by = c('start_interval','end_interval'))


limits_REL <- c(parse_date_time3('2021-03-19 09:00',tz='Etc/GMT-1'), parse_date_time3('2021-03-20 09:00',tz='Etc/GMT-1'))


Result[!is.na(Sonic) & st >= limits_REL[1] & st <= limits_REL[2],{0
	ggplot(.SD, aes(x=st,y=1/L,group=Sonic,colour=Sonic)) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1")
		,ymin=-Inf,ymax=Inf),colour = 'grey75',fill='grey75',alpha=0.01) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1")
 		,ymin=-Inf,ymax=Inf),colour = 'grey75',fill='grey75',alpha=0.01) +
	geom_line() +
	geom_hline(yintercept=0) +
	theme_bw(base_size=16)
}]


Fig_stab <- Result[!is.na(Sonic) & st >= limits_REL[1] & st <= limits_REL[2],{
	ggplot(.SD, aes(x=st,y=1/L,group=Sonic,colour=Sonic)) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1")
		,ymin=-Inf,ymax=Inf),colour = 'grey75',fill='grey75',alpha=0.01) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1")
 		,ymin=-Inf,ymax=Inf),colour = 'grey75',fill='grey75',alpha=0.01) +
	geom_line() +
	scale_y_continuous(limits=c(-0.1,0.1),breaks=seq(-0.25,0.25,0.05)) +
	geom_hline(yintercept=0) +
	geom_hline(yintercept=c(-0.02,-0.005,-0.002,0.02,0.005,0.002),lty=3,col='grey60') +
	annotate('text',label=c('very unstable','unstable','near neutral unstable','neutral','near neutral stable','stable','very stable')
		,x=limits_REL[1],y=c(-0.03,-0.01,-0.003,0.001,0.004,0.015,0.03),hjust=0) +
	theme_bw(base_size=16)
}]

ggsave(file.path(PfadFigures, "Stability_classes_during_release.png"), Fig_stab, width = 30, height = 30/1.57, units="cm")


Result$Sonic_ord <- factor(Result$Sonic, levels=c("SonicC","SonicA","Sonic2","SonicB"),labels=c("UA-UW","UA-50m","UA-100m","UA-150m"))

Fig_stab <- Result[!is.na(Sonic_ord) & st > parse_date_time3("19.03.2021 09:00",tz="Etc/GMT-1") & st < parse_date_time3("20.03.2021 09:00",tz="Etc/GMT-1"),{
	ggplot(.SD, aes(x=st,y=1/L,colour=Sonic_ord)) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1")
		,ymin=-Inf,ymax=Inf),colour='grey95',fill='grey95',alpha=0.15) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1")
 		,ymin=-Inf,ymax=Inf),colour='grey95',fill='grey95',alpha=0.15) +
	geom_line() +
	geom_point() +
	ylab(expression('1/L [m'^-1*']')) +
	xlab('') +
	scale_y_continuous(limits=c(-0.1,0.1),breaks=seq(-0.25,0.25,0.05)) +
	geom_hline(yintercept=0) +
	theme_bw(base_size=16) +
	theme(legend.title = element_blank(),panel.grid.major=element_line(size=0.6),panel.grid.minor=element_line(size=0.3))
}]

ggsave(file.path(PfadFigures, "Stability_during_release.png"), Fig_stab, width = 30, height = 30/1.57, units="cm")

Result[Sonic_ord == 'UA-UW' & st > parse_date_time3("19.03.2021 09:00",tz="Etc/GMT-1") & st < parse_date_time3("19.03.2021 18:00",tz="Etc/GMT-1") & Q_MFCex > 0,summary(L)]
Result[Sonic_ord == 'UA-UW' & st > parse_date_time3("19.03.2021 19:00",tz="Etc/GMT-1") & st < parse_date_time3("20.03.2021 10:00",tz="Etc/GMT-1") & Q_MFCex > 0,summary(L)]
# -410.2 -9.5 -57.7
# 14.8 184.1 74.2











SonicAll <- merge(merge(
            merge(cbind(as.data.table(SonicB_30min),start_interval=st(SonicB_30min),end_interval = et(SonicB_30min)), cbind(as.data.table(SonicA_30min),start_interval=st(SonicA_30min),end_interval = et(SonicA_30min)), by = c('start_interval', 'end_interval'), suffixes = c('', '.a'))
            , cbind(as.data.table(SonicC_30min),start_interval=st(SonicC_30min),end_interval = et(SonicC_30min)), by = c('start_interval', 'end_interval'), suffixes = c('', '.c'))
                , cbind(as.data.table(Sonic2_30min),start_interval=st(Sonic2_30min),end_interval = et(Sonic2_30min)), by = c('start_interval', 'end_interval'), suffixes = c('.b', '.2'))


Sonic_30min <- rbind(
	cbind(as.data.table(SonicB_30min),st=st(SonicB_30min),Sonic = 'SonicB')
	, cbind(as.data.table(SonicA_30min),st=st(SonicA_30min),Sonic = 'SonicA')
    , cbind(as.data.table(SonicC_30min),st=st(SonicC_30min),Sonic = 'SonicC')
    , cbind(as.data.table(Sonic2_30min),st=st(Sonic2_30min),Sonic = 'Sonic2')
)

Fig_stab_30min <- Sonic_30min[!is.na(Sonic) & st >= limits_REL[1] & st <= limits_REL[2],{
	ggplot(.SD, aes(x=st,y=1/L,group=Sonic,colour=Sonic)) +
	geom_rect(aes(xmin=parse_date_time3("19.03.2021 10:20",tz="Etc/GMT-1"),xmax=parse_date_time3("19.03.2021 16:50",tz="Etc/GMT-1")
		,ymin=-Inf,ymax=Inf),colour = 'grey75',fill='grey75',alpha=0.01) +
 	geom_rect(aes(xmin=parse_date_time3("19.03.2021 21:40",tz="Etc/GMT-1"),xmax=parse_date_time3("20.03.2021 07:00",tz="Etc/GMT-1")
 		,ymin=-Inf,ymax=Inf),colour = 'grey75',fill='grey75',alpha=0.01) +
	geom_line() +
	scale_y_continuous(limits=c(-0.1,0.1),breaks=seq(-0.25,0.25,0.05)) +
	geom_hline(yintercept=0) +
	geom_hline(yintercept=c(-0.02,-0.005,-0.002,0.02,0.005,0.002),lty=3,col='grey60') +
	annotate('text',label=c('very unstable','unstable','near neutral unstable','neutral','near neutral stable','stable','very stable')
		,x=limits_REL[1],y=c(-0.03,-0.01,-0.003,0.001,0.004,0.015,0.03),hjust=0) +
	theme_bw(base_size=16)
}]

ggsave(file.path(PfadFigures, "Stability_classes_during_release_30min.png"), Fig_stab_30min, width = 30, height = 30/1.57, units="cm")




############################################
############################################
#####                                  #####
#####    WindTrax Forward modelling    #####
#####                                  #####
############################################
############################################

WT[,R_GF16 := V6/V6[1],by = .(V5,V2)]
WT[,R_GF17 := V11/V11[1],by = .(V5,V2)]
WT[,R_GF18 := V16/V16[1],by = .(V5,V2)]
WT[,R_GF25 := V21/V21[1],by = .(V5,V2)]

# Fig_WT_conc <- WT[,{
# 	ggplot(.SD, aes(x=V3)) +
# 	geom_line(aes(y=V6, colour="GF16")) +
# 	geom_line(aes(y=V11, colour="GF17")) +
# 	geom_line(aes(y=V16, colour="GF18")) +
# 	geom_line(aes(y=V21, colour="GF25")) +
# 	xlab("Release Height [m]") +
# 	ylab("Concentration") +
# 	facet_grid(V5 ~ V2) +
# 	scale_colour_manual(name=NULL,values=c(CH4Cols)) +
# 	theme_classic(base_size = 18)
# }]
WT$stab <- factor(WT$V2, levels=c("SU","N","SS"),labels=c("slightly unstable","neutral","slightly stable"))
WT$height <- factor(WT$V5, levels=c(1.6,3),labels=c("1.6 m","3 m"))
# WT$GFs <- factor(WT$V2, levels=c("GF17","GF18","GF16","GF25"))


Fig_WT_conc_rel <- WT[,{
	ggplot(.SD, aes(x=V3)) +
	geom_line(aes(y=R_GF17*100, col="GF17")) +
	geom_line(aes(y=R_GF18*100, col="GF18")) +
	geom_line(aes(y=R_GF16*100, col="GF16")) +
	geom_line(aes(y=R_GF25*100, col="GF25")) +
	xlab("Release Height [m]") +
	ylab("Concentration of groud release [%]") +
	facet_grid(height ~ stab) +
	scale_colour_manual(name=NULL,values=CH4Cols[c("GF17","GF18","GF16","GF25")]) +
	theme_classic(base_size = 18)
}]

# ggsave(Fig_WT_conc,file=file.path(PfadFigures,"WindTrax_Forward_conc_abs.png"),width=13*1.6,height=13)
ggsave(Fig_WT_conc_rel,file=file.path(PfadFigures,"WindTrax_Forward_conc_rel.png"),width=30,height=30/1.7,units="cm")


#############
### Paper ###
#############

Fig_WT_conc_rel_Paper <- WT[height == '1.6 m',{
	ggplot(.SD, aes(x=V3)) +
	geom_line(aes(y=R_GF17*100, col="OP-2")) +
	geom_line(aes(y=R_GF18*100, col="OP-3")) +
	geom_line(aes(y=R_GF16*100, col="OP-4")) +
	geom_line(aes(y=R_GF25*100, col="OP-5")) +
	xlab("Release Height [m]") +
	ylab("Concentration of groud release [%]") +
	facet_grid( ~ stab) +
	scale_colour_manual(name=NULL,values=CH4Cols_OP[2:5]) +
	theme_bw(base_size = 18) + theme(strip.background =element_rect(fill="white"),panel.grid = element_blank())
}]

# ggsave(Fig_WT_conc,file=file.path(PfadFigures,"WindTrax_Forward_conc_abs.png"),width=13*1.6,height=13)
ggsave(Fig_WT_conc_rel_Paper,file=file.path(PfadFigures,"WindTrax_Forward_conc_rel_Paper.png"),width=30,height=30/2,units="cm")



# ######################################
# ######################################
# #####                            #####
# #####    Scaling + Similarity    #####
# #####                            #####
# ######################################
# ######################################
# set.seed(77)

# dt_A <- data.table(z = seq(0,100,length.out=12),x = abs(c(seq(0,20,length.out=6),seq(20,0,length.out=6)) + rnorm(12,mean=0,sd=1)),obs = "V1")
# dt_B <- data.table(z = seq(0,40,by=5),x = abs(c(seq(0,4),seq(3,0)) + rnorm(9,mean=0,sd=0.5)),obs="V2")

# dt_A[,zs := z/(max(z)/max(dt_B[,z]))]
# dt_B[,xs := x/max(x)*dt_A[,max(x)]]
# dt_AB <- rbind(dt_A,dt_B,fill=TRUE)

# dt_AB[,{
# 	ggplot(.SD, aes(x=x,y=z,colour=obs,group=obs)) +
# 	geom_line()
# }]

# png(file.path(PfadFigures,"scaling.png"),width=24,height=12,unit="in",res=300)
# par(mfrow=c(1,2),mar=c(1.8,1.5,0.2,0.2),mgp=c(0.5,0.5,0),cex=2)
# plot(z ~ x, dt_A,type="o",lty=2,ylim=c(0,105),xlim=c(0,25),xaxt="n",yaxt="n",ann=FALSE)
# lines(z ~ x, dt_B,pch=19,type="o")
# mtext(expression(bar("a")),1,line=0.2,cex=2)
# mtext("z",2,line=0.2,cex=2)
# text("A",x=1,y=105)
# #
# points(z ~ x, dt_B,pch=19,type="o")
# plot(zs ~ x, dt_A,type="o",lty=2,ylim=c(0,105),xlim=c(0,25),xaxt="n",yaxt="n",ann=FALSE)
# points(z ~ xs, dt_B,pch=19,type="o")
# mtext(expression(bar("a")*" / a"["*"]),1,line=0.2,cex=2)
# mtext(expression("z / z"["i"]),2,line=0.2,cex=2)
# text("B",x=1,y=105)
# dev.off()


##################################
##################################
#####                        #####
#####    Turbulence plots    #####
#####                        #####
##################################
##################################


sC_dir <- dir(file.path(PfadDaten,"Sonic","SonicC"), full.names = TRUE)
sC_list <- lapply(sC_dir[10], readWindMaster_ascii)
sC_data <- rbindlist(sC_list)
sC_Reynold <- melt(sC_data, id.vars="Time",measure.vars=c("u","v","w"), variable.name="Axis", value.name="WS")
setkey(sC_Reynold,Time)

# sC_Reynold[1:100,{
Fig_Turb <- sC_Reynold[9E5:(9E5+36000),{
	ggplot(.SD,aes(x=Time,y=WS,group=Axis)) +
	geom_hline(yintercept=0) +
	geom_line() +
	ylab(expression("Wind speed [m s"^-1*"]")) +
	xlab(NULL) +
	# xlab("Time") +
	facet_grid(Axis~.) +
    theme_classic(base_size = 18)
}]


ggsave(file.path(PfadFigures, "Turbulence_uvw.jpg"), Fig_Turb, width = 30, height = 30/1.57, units="cm")



#########################################################################################################################################################
#########################################################################################################################################################
#-----------------------------------------                                                                ----------------------------------------------#
#-----------------------------------------             TABLES         TABLES        TABLES                ----------------------------------------------#
#-----------------------------------------             TABLES         TABLES        TABLES                ----------------------------------------------#
#-----------------------------------------             TABLES         TABLES        TABLES                ----------------------------------------------#
#-----------------------------------------                                                                ----------------------------------------------#
#########################################################################################################################################################
#########################################################################################################################################################


########################################
########################################
#####                              #####
#####    Tabelle Windrichtungen    #####
#####                              #####
########################################
########################################

Result$Sonic_ord <- factor(Result$Sonic, levels=c("SonicC","SonicA","Sonic2","SonicB"),labels=c("3DUA-UW","3DUA-50m","3DUA-100m","3DUA-150m"))
Result[L > 0, stability := 'stable']
Result[L < 0, stability := 'unstable']
## während MK
meanWD <- data.table(Sonic=Result[!is.na(Sonic_ord),mean(WD),by=Sonic_ord][,Sonic_ord], # this line is purly to have the names of the Sonic in the right order
	meanWD=round((360 + atan2(Result[!is.na(Sonic_ord) & Campaign == "MK",mean(U_sonic * sin(WD * pi/180),na.rm=TRUE),by=Sonic_ord][,V1]
	,Result[!is.na(Sonic_ord) & Campaign == "MK",mean(U_sonic * cos(WD * pi/180),na.rm=TRUE),by=Sonic_ord][,V1]) * 180/pi) %% 360,1))

meanWS <- Result[!is.na(Sonic_ord) & Campaign == "MK",round(mean(U_sonic,na.rm=TRUE),1), by=Sonic_ord]

meanWD
meanWS


## during CH4 release in MC by stability and sonic
Result[!is.na(Sonic_ord) & Campaign == "MK" & Q_MFCex > 4, .(meanWS=round(mean(U_sonic,na.rm=TRUE),1),
	meanWD=round((360 + atan2(mean(U_sonic * sin(WD * pi/180),na.rm=TRUE),mean(U_sonic * cos(WD * pi/180),na.rm=TRUE)) * 180/pi) %% 360,1),
	meanUstar=round(mean(Ustar,na.rm=TRUE),2),meanL=round(mean(1/L,na.rm=TRUE),2)),by=.(stability,Sonic_ord)][order(-stability,Sonic_ord)]
# only sonic
Result[!is.na(Sonic_ord) & Campaign == "MK" & Q_MFCex > 4, .(meanWS=round(mean(U_sonic,na.rm=TRUE),1),
	meanWD=round((360 + atan2(mean(U_sonic * sin(WD * pi/180),na.rm=TRUE),mean(U_sonic * cos(WD * pi/180),na.rm=TRUE)) * 180/pi) %% 360,1),
	meanUstar=round(mean(Ustar,na.rm=TRUE),2),meanL=round(mean(1/L,na.rm=TRUE),2)),by=Sonic_ord]


# Install and load the circular package
library(circular)

# Calculate sd of the wind direction by transorming it first to circular, calcualte the sd and then back to degrees
Result[Campaign == 'MK' & !is.na(Sonic_ord),round(sd.circular(circular(WD, units ='degrees'),na.rm=TRUE) * (180/pi),1),by=Sonic_ord]
Result[Campaign == 'MK' & Q_MFC > 0 & !is.na(Sonic_ord),round(sd.circular(circular(WD, units ='degrees'),na.rm=TRUE) * (180/pi),1),by=Sonic_ord]



## während MK und Gasfreisetzung

meanWD_rel <- data.table(Sonic=Result[!is.na(Sonic_ord),mean(WD),by=Sonic_ord][,Sonic_ord], # this line is purly to have the names of the Sonic in the right order
	meanWD=round((360 + atan2(Result[!is.na(Sonic_ord) & Campaign == "MK" & Q_MFC > 0,mean(U_sonic * sin(WD * pi/180),na.rm=TRUE),by=Sonic_ord][,V1]
		,Result[!is.na(Sonic_ord) & Campaign == "MK" & Q_MFC > 0,mean(U_sonic * cos(WD * pi/180),na.rm=TRUE),by=Sonic_ord][,V1]) * 180/pi) %% 360,1))

meanWS_rel <- Result[!is.na(Sonic_ord) & Campaign == "MK" & Q_MFC > 0,round(mean(U_sonic,na.rm=TRUE),1),by=Sonic_ord]

meanWD_rel
meanWS_rel


#################
### u* values ###
#################

##### During the entire MK:
meanUstar <- Result[!is.na(Sonic_ord) & Campaign == "MK",round(mean(Ustar,na.rm=TRUE),2), by=Sonic_ord]

##### Only during gas release:
meanUstar_rel <- Result[!is.na(Sonic_ord) & Campaign == "MK" & Q_MFC > 0,round(mean(Ustar,na.rm=TRUE),2),by=Sonic_ord]

meanUstar
meanUstar_rel


#################################
#################################
#####                       #####
#####    Table data loss    #####
#####                       #####
#################################
#################################



Result[Campaign == "MK" & Q_MFC > 0 & Sonic=="SonicA",]
sonic_names <- Result[!is.na(Sonic),unique(Sonic)]

ls_16 <- ls_17 <- ls_18 <- ls_25 <- vector(mode="list",4)
for(j in 1:4){
	ls_16[[j]]	<- as.data.table(Result[Campaign == "MK" & Sonic == sonic_names[j] & Q_MFC > 0 & is.na(Q_GF16),.N]/
		Result[Campaign == "MK" & Sonic == sonic_names[j] & Q_MFC > 0 ,.N])
	ls_17[[j]]	<- as.data.table(Result[Campaign == "MK" & Sonic == sonic_names[j] & Q_MFC > 0 & is.na(Q_GF17),.N]/
		Result[Campaign == "MK" & Sonic == sonic_names[j] & Q_MFC > 0 ,.N])
	ls_18[[j]]	<- as.data.table(Result[Campaign == "MK" & Sonic == sonic_names[j] & Q_MFC > 0 & is.na(Q_GF18),.N]/
		Result[Campaign == "MK" & Sonic == sonic_names[j] & Q_MFC > 0 ,.N])
	ls_25[[j]]	<- as.data.table(Result[Campaign == "MK" & Sonic == sonic_names[j] & Q_MFC > 0 & is.na(Q_GF25),.N]/
		Result[Campaign == "MK" & Sonic == sonic_names[j] & Q_MFC > 0 ,.N])
}

data.loss <- round(rbind(rbindlist(ls_17), rbindlist(ls_18),rbindlist(ls_16),rbindlist(ls_25)),2)
colnames(data.loss) <- c("loss")
data.loss[,NAME := c(paste0(rep("GF17_",4),sonic_names),paste0(rep("GF18_",4),
	sonic_names),paste0(rep("GF16_",4),sonic_names),paste0(rep("GF25_",4),sonic_names))]

data.loss[grep('SonicA',NAME), ':=' (Sonic = 'Sonic50m', GF = c(paste0('GF',c('50m','100m','150m','200m'))))]
data.loss[grep('SonicB',NAME), ':=' (Sonic = 'Sonic150m', GF = c(paste0('GF',c('50m','100m','150m','200m'))))]
data.loss[grep('SonicC',NAME), ':=' (Sonic = 'SonicUW', GF = c(paste0('GF',c('50m','100m','150m','200m'))))]
data.loss[grep('Sonic2',NAME), ':=' (Sonic = 'Sonic100m', GF = c(paste0('GF',c('50m','100m','150m','200m'))))]

data.loss

data.loss[Sonic=='Sonic50m']
data.loss[Sonic=='Sonic100m']
data.loss[Sonic=='Sonic150m']
data.loss[Sonic=='SonicUW']



round(mean(data.loss[grep('SonicC',NAME),loss]),2)
round(mean(data.loss[grep('SonicA',NAME),loss]),2)
round(mean(data.loss[grep('Sonic2',NAME),loss]),2)
round(mean(data.loss[grep('SonicB',NAME),loss]),2)



######################################
######################################
#####                            #####
#####    Table recovery rates    #####
#####                            #####
######################################
######################################


sonic_names <- Result[!is.na(Sonic),unique(Sonic)]

summary(Result[Campaign == "MK" & Q_MFC > 0 ,Q_GF16/Q_MFC])
summary(Result[Campaign == "MK" & Q_MFC > 0 ,Q_GF17/Q_MFC])
summary(Result[Campaign == "MK" & Q_MFC > 0 ,Q_GF18/Q_MFC])
summary(Result[Campaign == "MK" & Q_MFC > 0 ,Q_GF25/Q_MFC])

ls_16 <- ls_17 <- ls_18 <- ls_25 <- vector(mode="list",4)
for(j in 1:4){
	ls_16[[j]]	<- as.data.table(t(as.numeric(summary(Result[Campaign == "MK" & Sonic == sonic_names[j] & Q_MFC > 0 ,Q_GF16/Q_MFC]))))
	ls_17[[j]]	<- as.data.table(t(as.numeric(summary(Result[Campaign == "MK" & Sonic == sonic_names[j] & Q_MFC > 0 ,Q_GF17/Q_MFC]))))
	ls_18[[j]]	<- as.data.table(t(as.numeric(summary(Result[Campaign == "MK" & Sonic == sonic_names[j] & Q_MFC > 0 ,Q_GF18/Q_MFC]))))
	ls_25[[j]]	<- as.data.table(t(as.numeric(summary(Result[Campaign == "MK" & Sonic == sonic_names[j] & Q_MFC > 0 ,Q_GF25/Q_MFC]))))
}

recovery_rate <- round(rbind(rbindlist(ls_17), rbindlist(ls_18),rbindlist(ls_16),rbindlist(ls_25)),2)
colnames(recovery_rate) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max","NA")
recovery_rate[,NAME := c(paste0(rep("GF17_",4),sonic_names),paste0(rep("GF18_",4),
	sonic_names),paste0(rep("GF16_",4),sonic_names),paste0(rep("GF25_",4),sonic_names))]

recovery_rate
recovery_rate[,.(Mean,Median,NAME)]
recovery_rate[grep('SonicC',NAME),.(Mean,Median,NAME)]
recovery_rate[grep('SonicA',NAME),.(Mean,Median,NAME)]
recovery_rate[grep('Sonic2',NAME),.(Mean,Median,NAME)]
recovery_rate[grep('SonicB',NAME),.(Mean,Median,NAME)]


c(min(recovery_rate[,Mean]),max(recovery_rate[,Mean]))

ls_16 <- ls_17 <- ls_18 <- ls_25 <- ls_26 <- vector(mode="list",4)
for(j in 1:4){
	ls_16[[j]]	<- as.data.table(t(as.numeric(summary(Result[Campaign == "QV3" & Sonic == sonic_names[j] & Q_MFC > 0 ,Q_GF16/Q_MFC]))))[,1:6]
	ls_17[[j]]	<- as.data.table(t(as.numeric(summary(Result[Campaign == "QV3" & Sonic == sonic_names[j] & Q_MFC > 0 ,Q_GF17/Q_MFC]))))[,1:6]
	ls_18[[j]]	<- as.data.table(t(as.numeric(summary(Result[Campaign == "QV3" & Sonic == sonic_names[j] & Q_MFC > 0 ,Q_GF18/Q_MFC]))))[,1:6]
	ls_25[[j]]	<- as.data.table(t(as.numeric(summary(Result[Campaign == "QV3" & Sonic == sonic_names[j] & Q_MFC > 0 ,Q_GF25/Q_MFC]))))[,1:6]
	ls_26[[j]]	<- as.data.table(t(as.numeric(summary(Result[Campaign == "QV3" & Sonic == sonic_names[j] & Q_MFC > 0 ,Q_GF26/Q_MFC]))))[,1:6]
}

recovery_rate <- round(rbind(rbindlist(ls_17), rbindlist(ls_18),rbindlist(ls_16),rbindlist(ls_25),rbindlist(ls_26)),2)
colnames(recovery_rate) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max")
recovery_rate[,NAME := c(paste0(rep("GF17_",4),sonic_names),paste0(rep("GF18_",4),
	sonic_names),paste0(rep("GF16_",4),sonic_names),paste0(rep("GF25_",4),sonic_names),paste0(rep("GF26_",4),sonic_names))]

recovery_rate
recovery_rate[,.(Mean,Median,NAME)]


sd(recovery_rate[1:12,Mean])
sd(recovery_rate[1:12,Median])
sd(recovery_rate[,Mean])
sd(recovery_rate[,Median])

Result[st > "2021-03-19" & st < "2021-03-19 18:00" & Sonic=="SonicB", ]
Result[st > "2021-03-19" & st < "2021-03-19 18:00" & Q_MFC > 0 & Sonic=="SonicB",mean(Q_MFC,na.rm=TRUE)*.N/6] * 0.4

max_mgm3 <- 0.04/ (10 * (8.31445983 * (2 + 273.15)) / 16.04 / 965)

max_mgm3/Result[st > "2021-03-19" & st < "2021-03-19 18:00" & Sonic=="SonicB", mean(CQ_GF17,na.rm=TRUE)] * mgs_to_kgh




#################################################
### Distances to the tree instead of the barn ###
#################################################

Baum
Sensors_MK
Sonics_MK

d = √((x2 - x1)² + (y2 - y1)²)

Tree_distance <- sapply(1:4, function(i) sqrt((Baum[[4]] - Sonics_MK[[4]][i])^2 + (Baum[[5]] - Sonics_MK[[5]][i])^2))[1:3]
Tree_distance
Tree_distance/15

Senosors_centre <- rbindlist(lapply(c(1,3,5,7,9), function(i) data.table((Sensors_MK[[4]][i] + Sensors_MK[[4]][i+1])/2, (Sensors_MK[[5]][i] + Sensors_MK[[5]][i+1])/2)))
	
Tree_distance2 <- data.table(OP=unique(Sensors_MK[[1]]),
	rbindlist(lapply(1:5, function(i) sqrt((Baum[[4]] - Senosors_centre[,1][i])^2 + (Baum[[5]] - Senosors_centre[,2][i])^2)))
)

Tree_distance2[,hx := round(V1/15,1)]
Tree_distance2




################################################
### Precision calculation for the GasFinders ###
################################################

# Using the Equation 2 from Häni et al. (2020).
# Precision = 2.9 * MAD(dC) / (l_path * sqrt(2))

GFs10min
GFs10min[,'Campaign'] <- NA
GFs10min[QV2,'Campaign'] <- 'IC1'
GFs10min[MK,'Campaign'] <- 'MC'
GFs10min[QV3,'Campaign'] <- 'IC2'

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
OP_Precision
Conc_melt[!is.na(MFC) & !is.na(dCppm),.(N=.N),by=.(Campaign,variable)] # Number of measurements

OP_Precision_all <- sapply(Conc_melt[,unique(variable)], function(k){
		out <- Conc_melt[variable == k, round(2.9*mad(get('dCppm'),na.rm=TRUE) * 110 /sqrt(2),1)]
		names(out) <- k
		return(out)
	})	
OP_Precision_all
Conc_melt[!is.na(MFC) & !is.na(dCppm),.(N=.N),by=.(variable)] # Number of measurements



Conc[!is.na(MFC) & Campaign == 'MC' & !is.na(dCppm_OP2h),.N]


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

dt_Conc[Q_MFCex > 4 & !is.na(dCppm) & Campaign == 'MC',.(mean_ppmm=mean(dCppm*100),median_ppmm=median(dCppm*100),N=.N),by=variable]
dt_Conc[!is.na(Rate) & Campaign == 'MC',.(mean_ppmm=round(mean(dCppm*100),1),median_ppmm=round(median(dCppm*100),1),sd=round(sd(dCppm*100),1),N=.N),by=variable]
dt_Conc[!is.na(Rate) & Campaign == 'MC',.(mean_ppmm=round(mean(dCppm*100),1),median_ppmm=round(median(dCppm*100),1),sd=round(sd(dCppm*100),1),N=.N),by=.(variable,stability)]


Conc[Campaign == 'MC' & !is.na(dCppm_OP2h) & !is.na(MFC)]


---> noch Plots von anderen Scripts hier einfügen und allenfalls code dort löschen oder sogar ganzes Script


#########################################
### Plots für Diskussion mit Albrecht ###
#########################################

GFs_1min[,"dC16"] <- GFs_1min[,"GF16"] - GFs_1min[,"GF26"]
GFs_1min[,"dC18"] <- GFs_1min[,"GF18"] - GFs_1min[,"GF26"]

## Select Timerange
isub <- "20.03.2021 02:30 - 20.03.2021 05:30"
xabline <- rep(parse_date_time3("20.03.2021 02:00",tz="Etc/GMT-1"),16*5) + (seq(length.out = 16*5, by=60*5) - 1)

## SonicA
png(file.path(PfadFigures,"SonicA_Albrecht.png"),width=12,height=8,unit="in",res=300)
par(mfrow=c(3,1),mar=c(3,4,1,0))
plot(SonicA_1min[isub,"WD"],ylim=c(30,70),ylab="wind direction",grid.y=seq(0,90,by=5), col="orange",lwd=2)
abline(v=xabline,col="lightgrey",lty=3)
legend("topleft",legend="SonicA - 50 m")
#
plot(SonicA_1min[isub,"U_sonic"],ylim=c(1,6),ylab="wind speed [m/s]",col="darkgreen",lwd=2)
abline(v=xabline,col="lightgrey",lty=3)
legend("topleft",legend="SonicA - 50 m")
#
plot(GFs_1min[isub,"dC16"],col="red",ylim=c(0,0.7),ylab="dC [mg/m3]")
lines(GFs_1min[isub,"dC18"],col="blue")
abline(v=xabline,col="lightgrey",lty=3)
legend("topleft",legend=c("GF18 - 100 m","GF16 - 150 m"),fill=c("blue","red"))
dev.off()

## SonicB
png(file.path(PfadFigures,"SonicB_Albrecht.png"),width=12,height=8,unit="in",res=300)
par(mfrow=c(3,1),mar=c(3,4,1,0))
plot(SonicB_1min[isub,"WD"],ylim=c(30,70),ylab="wind direction",grid.y=seq(0,90,by=5), col="orange",lwd=2)
abline(v=xabline,col="lightgrey",lty=3)
legend("topleft",legend="SonicB - 150 m")
#
plot(SonicB_1min[isub,"U_sonic"],ylim=c(1,6),ylab="wind speed [m/s]",col="darkgreen",lwd=2)
abline(v=xabline,col="lightgrey",lty=3)
legend("topleft",legend="SonicB - 150 m")
#
plot(GFs_1min[isub,"dC16"],col="red",ylim=c(0,0.7),ylab="dC [mg/m3]")
lines(GFs_1min[isub,"dC18"],col="blue")
abline(v=xabline,col="lightgrey",lty=3)
legend("topleft",legend=c("GF18 - 100 m","GF16 - 150 m"),fill=c("blue","red"))
dev.off()

## SonicC
png(file.path(PfadFigures,"SonicC_Albrecht.png"),width=12,height=8,unit="in",res=300)
par(mfrow=c(3,1),mar=c(3,4,1,0))
plot(SonicC_1min[isub,"WD"],ylim=c(30,70),ylab="wind direction",grid.y=seq(0,90,by=5), col="orange",lwd=2)
abline(v=xabline,col="lightgrey",lty=3)
legend("topleft",legend="SonicC - background")
#
plot(SonicC_1min[isub,"U_sonic"],ylim=c(1,6),ylab="wind speed [m/s]",col="darkgreen",lwd=2)
abline(v=xabline,col="lightgrey",lty=3)
legend("topleft",legend="SonicC - background")
#
plot(GFs_1min[isub,"dC16"],col="red",ylim=c(0,0.7),ylab="dC [mg/m3]")
lines(GFs_1min[isub,"dC18"],col="blue")
abline(v=xabline,col="lightgrey",lty=3)
legend("topleft",legend=c("GF18 - 100 m","GF16 - 150 m"),fill=c("blue","red"))
dev.off()

## Sonic2
png(file.path(PfadFigures,"Sonic2_Albrecht.png"),width=12,height=8,unit="in",res=300)
par(mfrow=c(3,1),mar=c(3,4,1,0))
plot(Sonic2_1min[isub,"WD"],ylim=c(30,70),ylab="wind direction",grid.y=seq(0,90,by=5), col="orange",lwd=2)
abline(v=xabline,col="lightgrey",lty=3)
legend("topleft",legend="Sonic2 - 100 m")
#
plot(Sonic2_1min[isub,"U_sonic"],ylim=c(1,6),ylab="wind speed [m/s]",col="darkgreen",lwd=2)
abline(v=xabline,col="lightgrey",lty=3)
legend("topleft",legend="Sonic2 - 100 m")
#
plot(GFs_1min[isub,"dC16"],col="red",ylim=c(0,0.7),ylab="dC [mg/m3]")
lines(GFs_1min[isub,"dC18"],col="blue")
abline(v=xabline,col="lightgrey",lty=3)
legend("topleft",legend=c("GF18 - 100 m","GF16 - 150 m"),fill=c("blue","red"))
dev.off()

#########################################

## All Sonics
png(file.path(PfadFigures,"SonicAll_Albrecht.png"),width=12,height=8,unit="in",res=300)
par(mfrow=c(2,1),mar=c(3,4,1,0))
plot(SonicA_1min[isub,"WD"],ylim=c(30,70),ylab="wind direction",grid.y=seq(0,90,by=5), col="red",lwd=2)
lines(SonicB_1min[isub,"WD"], col="blue",lwd=2)
lines(SonicC_1min[isub,"WD"], col="black",lwd=2)
lines(Sonic2_1min[isub,"WD"], col="#FFEA00",lwd=2)
abline(v=xabline,col="lightgrey",lty=3)
# legend("topleft",legend=c("SonicA - 50 m","Sonic2 - 100 m","SonicB - 150 m","SonicC - background"),fill=c("red","#FFEA00","blue","black"))
#
plot(SonicA_1min[isub,"U_sonic"],ylim=c(1,6),ylab="wind speed [m/s]",col="red",lwd=2)
lines(SonicB_1min[isub,"U_sonic"], col="blue",lwd=2)
lines(SonicC_1min[isub,"U_sonic"], col="black",lwd=2)
lines(Sonic2_1min[isub,"U_sonic"], col="#FFEA00",lwd=2)
abline(v=xabline,col="lightgrey",lty=3)
legend("topright",legend=c("SonicA - 50 m","Sonic2 - 100 m","SonicB - 150 m","SonicC - background"),fill=c("red","#FFEA00","blue","black"))
dev.off()


## Select Timerange
isub2 <- "19.03.2021 13:00 - 19.03.2021 16:30"
xabline2 <- rep(parse_date_time3("19.03.2021 13:00",tz="Etc/GMT-1"),16*5) + (seq(length.out = 16*5, by=60*5) - 1)

## SonicA
png(file.path(PfadFigures,"SonicA_Albrecht_Nami.png"),width=12,height=8,unit="in",res=300)
par(mfrow=c(3,1),mar=c(3,4,1,0))
plot(SonicA_1min[isub2,"WD"],ylim=c(0,100),ylab="wind direction",grid.y=seq(0,110,by=5), col="orange",lwd=2)
abline(v=xabline2,col="lightgrey",lty=3)
legend("topleft",legend="SonicA - 50 m")
#
plot(SonicA_1min[isub2,"U_sonic"],ylim=c(1,7.5),ylab="wind speed [m/s]",col="darkgreen",lwd=2)
abline(v=xabline2,col="lightgrey",lty=3)
legend("topleft",legend="SonicA - 50 m")
#
plot(GFs_1min[isub2,"dC16"],col="red",ylim=c(-0.1,0.6),ylab="dC [mg/m3]")
lines(GFs_1min[isub2,"dC18"],col="blue")
abline(v=xabline2,col="lightgrey",lty=3)
legend("topleft",legend=c("GF18 - 100 m","GF16 - 150 m"),fill=c("blue","red"))
dev.off()

## SonicB
png(file.path(PfadFigures,"SonicB_Albrecht_Nami.png"),width=12,height=8,unit="in",res=300)
par(mfrow=c(3,1),mar=c(3,4,1,0))
plot(SonicB_1min[isub2,"WD"],ylim=c(0,100),ylab="wind direction",grid.y=seq(0,110,by=5), col="orange",lwd=2)
abline(v=xabline2,col="lightgrey",lty=3)
legend("topleft",legend="SonicB - 150 m")
#
plot(SonicB_1min[isub2,"U_sonic"],ylim=c(1,7.5),ylab="wind speed [m/s]",col="darkgreen",lwd=2)
abline(v=xabline2,col="lightgrey",lty=3)
legend("topleft",legend="SonicB - 150 m")
#
plot(GFs_1min[isub2,"dC16"],col="red",ylim=c(-0.1,0.6),ylab="dC [mg/m3]")
lines(GFs_1min[isub2,"dC18"],col="blue")
abline(v=xabline2,col="lightgrey",lty=3)
legend("topleft",legend=c("GF18 - 100 m","GF16 - 150 m"),fill=c("blue","red"))
dev.off()

## SonicC
png(file.path(PfadFigures,"SonicC_Albrecht_Nami.png"),width=12,height=8,unit="in",res=300)
par(mfrow=c(3,1),mar=c(3,4,1,0))
plot(SonicC_1min[isub2,"WD"],ylim=c(0,100),ylab="wind direction",grid.y=seq(0,110,by=5), col="orange",lwd=2)
abline(v=xabline2,col="lightgrey",lty=3)
legend("topleft",legend="SonicC - background")
#
plot(SonicC_1min[isub2,"U_sonic"],ylim=c(1,7.5),ylab="wind speed [m/s]",col="darkgreen",lwd=2)
abline(v=xabline2,col="lightgrey",lty=3)
legend("topleft",legend="SonicC - background")
#
plot(GFs_1min[isub2,"dC16"],col="red",ylim=c(-0.1,0.6),ylab="dC [mg/m3]")
lines(GFs_1min[isub2,"dC18"],col="blue")
abline(v=xabline2,col="lightgrey",lty=3)
legend("topleft",legend=c("GF18 - 100 m","GF16 - 150 m"),fill=c("blue","red"))
dev.off()

## Sonic2
png(file.path(PfadFigures,"Sonic2_Albrecht_Nami.png"),width=12,height=8,unit="in",res=300)
par(mfrow=c(3,1),mar=c(3,4,1,0))
plot(Sonic2_1min[isub2,"WD"],ylim=c(0,100),ylab="wind direction",grid.y=seq(0,110,by=5), col="orange",lwd=2)
abline(v=xabline2,col="lightgrey",lty=3)
legend("topleft",legend="Sonic2 - 100 m")
#
plot(Sonic2_1min[isub2,"U_sonic"],ylim=c(1,7.5),ylab="wind speed [m/s]",col="darkgreen",lwd=2)
abline(v=xabline2,col="lightgrey",lty=3)
legend("topleft",legend="Sonic2 - 100 m")
#
plot(GFs_1min[isub2,"dC16"],col="red",ylim=c(-0.1,0.6),ylab="dC [mg/m3]")
lines(GFs_1min[isub2,"dC18"],col="blue")
abline(v=xabline2,col="lightgrey",lty=3)
legend("topleft",legend=c("GF18 - 100 m","GF16 - 150 m"),fill=c("blue","red"))
dev.off()

#########################################

## All Sonics
png(file.path(PfadFigures,"SonicAll_Albrecht_Nami.png"),width=12,height=8,unit="in",res=300)
par(mfrow=c(2,1),mar=c(3,4,1,0))
plot(SonicA_1min[isub2,"WD"],ylim=c(0,100),ylab="wind direction",grid.y=seq(0,110,by=5), col="red",lwd=2)
lines(SonicB_1min[isub2,"WD"], col="blue",lwd=2)
lines(SonicC_1min[isub2,"WD"], col="black",lwd=2)
lines(Sonic2_1min[isub2,"WD"], col="#FFEA00",lwd=2)
abline(v=xabline2,col="lightgrey",lty=3)
# legend("topleft",legend=c("SonicA - 50 m","Sonic2 - 100 m","SonicB - 150 m","SonicC - background"),fill=c("red","#FFEA00","blue","black"))
#
plot(SonicA_1min[isub2,"U_sonic"],ylim=c(1,7.5),ylab="wind speed [m/s]",col="red",lwd=2)
lines(SonicB_1min[isub2,"U_sonic"], col="blue",lwd=2)
lines(SonicC_1min[isub2,"U_sonic"], col="black",lwd=2)
lines(Sonic2_1min[isub2,"U_sonic"], col="#FFEA00",lwd=2)
abline(v=xabline2,col="lightgrey",lty=3)
legend("topright",legend=c("SonicA - 50 m","Sonic2 - 100 m","SonicB - 150 m","SonicC - background"),fill=c("red","#FFEA00","blue","black"))
dev.off()
