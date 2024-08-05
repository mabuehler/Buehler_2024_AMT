
################################################
################################################
#####                                      #####
#####    Calculations for contour plots    #####
#####                                      #####
################################################
################################################

# Author: Marcel Bühler
# Date: August 5, 2024
# Contact: mb@bce.au.dk or Christoph Häni christoph.haeni@bfh.ch
# Description: This script calculates the plume contours in the XY and XZ plane. This script is not necessary to reproduce the findings of the publication.
#
# Note: This code was written by Marcel Bühler and is intended to follow the publication 'Applicability of the inverse dispersion method to measure emissions from animal housings' in AMT. 
# Please feel free to use and modify it, but attribution is appreciated.


#################
### Libraries ###
#################

library(ibts)
library(RgoogleMaps)
library(bLSmodelR)


#############
### Paths ###
#############

PathData <- "Path to /data"		
PathRSaves <- "Path to /RSaves"
PathFigures <- 'Path to /Figures'
Cat.Path <- 'Path to /Catalogs'


#################
### Functions ###
#################

source("https://raw.githubusercontent.com/hafl-gel/gel-scripts/main/wgs84-ch1903.r")
source(file.path(file.path(dirname(PathRSaves),"Other/contourXY.r")))
source(file.path(file.path(dirname(PathRSaves),"Other/contourXZ.r")))


#################
### load data ###
#################

	# Geometry
	load(file.path(PathRSaves,"Geometry.RData"))
	# Google map
	STO_Map <- ReadMapTile(file.path(PathFigures,"STO_GoogleMaps.png"))
	# Sonics
	load(file.path(PathRSaves, "Sonics_10min.RData"))


######################
### data treatment ###
######################

## convert Sensors/Sources to Map x/y
Source_xy <- ch_to_map(STO_Map,Source)
Sensors_xy <- ch_to_map(STO_Map,Sensors_MC)

## Define CH4 release in MC (+ 10 min)
	indSub <- c("19.03.2021 10:20 - 19.03.2021 16:40","19.03.2021 21:40 - 20.03.2021 07:00")
## Define variables for number chruncher	
	N_traj <- 1E6
	ncores <- 88
	MaxFetch <- 400


##################################
### Calculation of countour XY ###
##################################

dx <- 5 # 5m spacing
dy <- 5 # 5m spacing

centre <- colMeans(Source[,c(2,3)])

xlim <- c(centre[1]-250,centre[1]+145)
ylim <- c(centre[2]-200,centre[2]+70)

#### all (unfiltered) data during the release experiment

	genSonicC <- genInterval(cbind(
			setNames(as.data.frame(SonicC[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
			,st=st(SonicC[indSub])
			,et=et(SonicC[indSub])
			,Sonic="SonicC")
		,MaxFetch=MaxFetch,N=N_traj)

	genSonicA <- genInterval(na.omit(cbind( ## for some reason there is a NA value and thus the thing does not work
			setNames(as.data.frame(SonicA[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
			,st=st(SonicA[indSub])
			,et=et(SonicA[indSub])
			,Sonic="SonicA"))
		,MaxFetch=MaxFetch,N=N_traj)

	genSonicB <- genInterval(cbind(
			setNames(as.data.frame(SonicB[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
			,st=st(SonicB[indSub])
			,et=et(SonicB[indSub])
			,Sonic="SonicB")
		,MaxFetch=MaxFetch,N=N_traj)

	genSonic2 <- genInterval(cbind(
			setNames(as.data.frame(Sonic2[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
			,st=st(Sonic2[indSub])
			,et=et(Sonic2[indSub])
			,Sonic="Sonic2")
		,MaxFetch=MaxFetch,N=N_traj)


#################################
### Run contour plot function ###
#################################

png(file = tempfile(fileext = '.png')) # save temp. png file so that it doesn't open a plot window
 XY_SonicC <- contourXY(Source, z.meas = 1.6,xlim=xlim,ylim=ylim, #without map. either define an origin (then xlim and ylim are == 100, or define xlim and ylim)
    Interval = genSonicC[,], Cat.Path = Cat.Path, dx = dx, dy = dy,
    fill = TRUE,ncores=ncores)
dev.off() # close png
  saveRDS(XY_SonicC,paste0(PfadRSaves,"/XY_SonicC_MC.rds"))
#
png(file = tempfile(fileext = '.png'))
 XY_SonicB <- contourXY(Source, z.meas = 1.6,xlim=xlim,ylim=ylim,
    Interval = genSonicB[,], Cat.Path = Cat.Path, dx = dx, dy = dy,
    fill = TRUE,ncores=ncores)
dev.off()
  saveRDS(XY_SonicB,paste0(PfadRSaves,"/XY_SonicB_MC.rds"))
#
png(file = tempfile(fileext = '.png'))
 XY_SonicA <- contourXY(Source, z.meas = 1.6,xlim=xlim,ylim=ylim,
    Interval = genSonicA[,], Cat.Path = Cat.Path, dx = dx, dy = dy,
    fill = TRUE,ncores=ncores)
dev.off()
  saveRDS(XY_SonicA,paste0(PfadRSaves,"/XY_SonicA_MC.rds"))

png(file = tempfile(fileext = '.png'))
 XY_Sonic2 <- contourXY(Source, z.meas = 1.6,xlim=xlim,ylim=ylim,
    Interval = genSonic2[,], Cat.Path = Cat.Path, dx = dx, dy = dy,
    fill = TRUE,ncores=ncores)
dev.off()
  saveRDS(XY_Sonic2,paste0(PfadRSaves,"/XY_Sonic2_MC.rds"))


# # #### Plotting on Google Maps (This code might not work anymore)

#  XY_SonicC <- contourXY(Source, z.meas = 1.6, MyMap=STO_Map,
#     coord2wgs84 = function(x,y){
#       # browser()
#       list(lat=CH.to.WGS.lat(x,y),lon=CH.to.WGS.lng(x,y))},
#     Interval = genSonicC[,], Cat.Path = Cat.Path, dx = dx, dy = dy,
#     fill = TRUE,ncores=ncores)
#   saveRDS(XY_SonicC,paste0(PfadRSaves,"/XY_SonicC_MC_map.rds"))


#  XY_SonicB <- contourXY(Source, z.meas = 1.6, MyMap=STO_Map,
#     coord2wgs84 = function(x,y){
#       # browser()
#       list(lat=CH.to.WGS.lat(x,y),lon=CH.to.WGS.lng(x,y))},
#     Interval = genSonicB[,], Cat.Path = Cat.Path, dx = dx, dy = dy,
#     fill = TRUE,ncores=ncores)
#   saveRDS(XY_SonicB,paste0(PfadRSaves,"/XY_SonicB_MC_map.rds"))


#  XY_SonicA <- contourXY(Source, z.meas = 1.6, MyMap=STO_Map,
#     coord2wgs84 = function(x,y){
#       # browser()
#       list(lat=CH.to.WGS.lat(x,y),lon=CH.to.WGS.lng(x,y))},
#     Interval = genSonicA[,], Cat.Path = Cat.Path, dx = dx, dy = dy,
#     fill = TRUE,ncores=ncores)
#   saveRDS(XY_SonicA,paste0(PfadRSaves,"/XY_SonicA_MC_map.rds"))


#  XY_Sonic2 <- contourXY(Source, z.meas = 1.6, MyMap=STO_Map,
#     coord2wgs84 = function(x,y){
#       # browser()
#       list(lat=CH.to.WGS.lat(x,y),lon=CH.to.WGS.lng(x,y))},
#     Interval = genSonic2[,], Cat.Path = Cat.Path, dx = dx, dy = dy,
#     fill = TRUE,ncores=ncores)
#   saveRDS(XY_Sonic2,paste0(PfadRSaves,"/XY_Sonic2_MC_map.rds"))


##################################
### Calculation of countour XZ ###
##################################

## this code is nowhere used

pQ <-colMeans(Source[,2:3]) ## Punkt Quelle (Mitte)
pM <- colMeans(Sensors_MC[Sensors_MC[,1] == c("GF25"),4:5]) # Punkt Mitte GF Pfad

WD <- 53 # angle
length_g <- 171 # define a random length of the side opposite (the angle) to determine the adjacent
length_a <- length_g / tan(WD*pi/180)
pWD <- pQ - c(length_g,length_a) # point extension WD

# Function
lineM <- lm(c(pM[2],pQ[2]) ~ c(pM[1],pQ[1]))
lineWD <- lm(c(pWD[2],pQ[2]) ~ c(pWD[1],pQ[1]))

dm <- 30 # extension

# end points line
p1M <- c(pQ[1]+dm,(pQ[1]+dm) * coef(lineM)[2] + coef(lineM)[1])
p2M <- c(pM[1]-dm,(pM[1]-dm) * coef(lineM)[2] + coef(lineM)[1]) 

p1WD <- c(pQ[1]+dm,(pQ[1]+dm) * coef(lineWD)[2] + coef(lineWD)[1])
p2WD <- c(pWD[1]-dm,(pWD[1]-dm) * coef(lineWD)[2] + coef(lineWD)[1]) 

plot(Source,Sensors_MC)
lines(c(p1M[1],p2M[1]),c(p1M[2],p2M[2]),col="orange")
lines(c(p1WD[1],p2WD[1]),c(p1WD[2],p2WD[2]),col="lightblue")

lengthM <- round(((p1M[2]-p2M[2])^2+(p1M[1]-p2M[1])^2)^0.5,0)
lengthWD <- round(((p1WD[2]-p2WD[2])^2+(p1WD[1]-p2WD[1])^2)^0.5,0)
	

### run bLS
ncores <- 88

# Middle OP Path (45°)
XZ_out_M_SonicA <- contourXZ(Source,p1=p1M,p2=p2M,Interval=genSonicA[,],
  z.upper=15,nz=30,nx=lengthM,ncores=ncores,Cat.Path=Cat.Path)
saveRDS(XZ_out_M_SonicA,paste0(PfadRSaves,"/XZ_out_M_SonicA.rds"))

XZ_out_M_SonicB <- contourXZ(Source,p1=p1M,p2=p2M,Interval=genSonicB[,],
  z.upper=15,nz=30,nx=lengthM,ncores=ncores,Cat.Path=Cat.Path)
saveRDS(XZ_out_M_SonicB,paste0(PfadRSaves,"/XZ_out_M_SonicB.rds"))

XZ_out_M_SonicC <- contourXZ(Source,p1=p1M,p2=p2M,Interval=genSonicC[,],
  z.upper=15,nz=30,nx=lengthM,ncores=ncores,Cat.Path=Cat.Path)
saveRDS(XZ_out_M_SonicC,paste0(PfadRSaves,"/XZ_out_M_SonicC.rds"))

XZ_out_M_Sonic2 <- contourXZ(Source,p1=p1M,p2=p2M,Interval=genSonic2[,],
  z.upper=15,nz=30,nx=lengthM,ncores=ncores,Cat.Path=Cat.Path)
saveRDS(XZ_out_M_Sonic2,paste0(PfadRSaves,"/XZ_out_M_Sonic2.rds"))

# Mean wind direction (53°)
XZ_out_WD_SonicA <- contourXZ(Source,p1=p1WD,p2=p2WD,Interval=genSonicA[,],
  z.upper=15,nz=30,nx=lengthWD,ncores=ncores,Cat.Path=Cat.Path)
saveRDS(XZ_out_WD_SonicA,paste0(PfadRSaves,"/XZ_out_WD_SonicA.rds"))

XZ_out_WD_SonicB <- contourXZ(Source,p1=p1WD,p2=p2WD,Interval=genSonicB[,],
  z.upper=15,nz=30,nx=lengthWD,ncores=ncores,Cat.Path=Cat.Path)
saveRDS(XZ_out_WD_SonicB,paste0(PfadRSaves,"/XZ_out_WD_SonicB.rds"))

XZ_out_WD_SonicC <- contourXZ(Source,p1=p1WD,p2=p2WD,Interval=genSonicC[,],
  z.upper=15,nz=30,nx=lengthWD,ncores=ncores,Cat.Path=Cat.Path)
saveRDS(XZ_out_WD_SonicC,paste0(PfadRSaves,"/XZ_out_WD_SonicC.rds"))

XZ_out_WD_Sonic2 <- contourXZ(Source,p1=p1WD,p2=p2WD,Interval=genSonic2[,],
  z.upper=15,nz=30,nx=lengthWD,ncores=ncores,Cat.Path=CatPath)
saveRDS(XZ_out_WD_Sonic2,paste0(PfadRSaves,"/XZ_out_WD_Sonic2.rds"))

