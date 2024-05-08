
#############################################################
#############################################################
#####                                                   #####
#####    Emissionen berechnen und Filterung anwenden    #####
#####                                                   #####
#############################################################
#############################################################


#################################
### Header (Pfade, Libraries) ###
#################################

	library(ibts)
	library(RgoogleMaps)
	library(bLSmodelR)
	
PathData <- "Path to /data"		
PathRSaves <- "Path to /RSaves"

	PfadFigures <- file.path(dirname(PfadRSaves),"Figures")
	source("~/repos/3_Scripts/gel-scripts/wgs84-ch1903.r")
	Cat.Path <- paste0(PfadDaten,"/Catalogs")

##################
### Funktionen ###
##################

calcdC <- function(data, sens, bgdname, new_sens=sens, source=FALSE, subset = data[, seq_len(.N)]){
data[subset, paste0("dC_", new_sens) := get(sens) - get(bgdname) - if(source == FALSE) {0} else {
		get(paste("dC",source,sens,sep="_")) - get(paste("dC",source,bgdname,sep="_"))
	}]
}

calcQ <- function(data, sens, bgdname, new_sens = sens, source = FALSE, nthresh = 2, subset = data[, seq_len(.N)]){
		mult <- if(grepl("^GF",sens)) mgs_to_kgh else ugs_to_kgd
		data[subset, paste0("Q_", new_sens) := (get(sens) - get(bgdname) - if(source==FALSE) {0} else {
			get(paste("dC",source,sens,sep = "_")) - get(paste("dC",source,bgdname,sep = "_"))
			}) / ifelse(get(paste0("N_TD_", sens)) <= nthresh,
			NA_real_, get(paste0("CQ_",sens)))*mult]
	}

lines_sec2xy <- function(xyMK,sensor,node=1,wd,col="lightblue",lwd=2,...){
	GF <- xyMK[xyMK[,1] %in% sensor,]
	sens <- as.numeric(GF[GF[,3] == node,4:5])
	b <- tan((90 - wd)/180*pi)
	x <- if(wd <= 180) 600 else -600
	y <- sens[2] - (sens[1] - x)*b
	lines(c(sens[1],x),c(sens[2],y),col=col,lwd=lwd,...)
}

	Col_CH4 <- "#1e88e5"
	library(RColorBrewer)
	# display.brewer.all()
	CH4Cols <- brewer.pal(5,"Dark2")
	names(CH4Cols) <- paste0("GF",c(16:18,25:26))

	mgs_to_kgd <- (24*3600)/1E6
	mgs_to_kgh <- (3600)/1E6
	mgm2s_to_gm2h <- 3600/1E3


	### MKs
	QV1 <- " to 28.01.2021 12:00"
	QV2 <- "05.03.2021 to 10.03.2021 18:00"
	QV1u2 <- " to 10.03.2021 18:00"
	MK <- "18.03.2021 11:00 - 21.03.2021 14:00"
	QV3	<- "21.03.2021 14:00 to "

#######################
### Laden der Daten ###
#######################

	# Wetterstationen
	WS700_1 <- readRDS(file=file.path(PfadRSaves,"WS700_1_STO.rds"))
	WS700_2 <- readRDS(file=file.path(PfadRSaves,"WS700_2_STO.rds"))
	# Konzentrationen
	Conc_wo_10min <- readRDS(file=file.path(PfadRSaves,"STO_Conc_wo_10min.rds"))
	Conc_P23_10min <- readRDS(file=file.path(PfadRSaves,"STO_Conc_P23_10min.rds"))
	Conc_P4_10min <- readRDS(file=file.path(PfadRSaves,"STO_Conc_P4_10min.rds"))
	Conc_P6_10min <- readRDS(file=file.path(PfadRSaves,"STO_Conc_P6_10min.rds"))
	# MFC
	MFC <- readRDS(file=file.path(PfadRSaves,"STO_MFC.rds"))
	# Drucksensor
	CC_Sensor <- readRDS(file=file.path(PfadRSaves,"STO_CC_Sensor.rds"))
	### bLS Resultate
	# normal
	bLS_STO_MK_10min <- readRDS(file = file.path(PfadRSaves, "STO_bLS_MK_10min.rds"))
	bLS_STO_MK_30min <- readRDS(file = file.path(PfadRSaves, "STO_bLS_MK_30min.rds"))
	bLS_STO_QV3_10min <- readRDS(file = file.path(PfadRSaves, "STO_bLS_QV3_10min.rds"))
	bLS_STO_QV3_30min <- readRDS(file = file.path(PfadRSaves, "STO_bLS_QV3_30min.rds"))
	# WDvar
	bLS_STO_SonicA_WDvar_10min <- readRDS(file = file.path(PfadRSaves, "STO_bLS_SonicA_WDvar_10min.rds"))
	bLS_STO_SonicB_WDvar_10min <- readRDS(file = file.path(PfadRSaves, "STO_bLS_SonicB_WDvar_10min.rds"))
	bLS_STO_SonicC_WDvar_10min <- readRDS(file = file.path(PfadRSaves, "STO_bLS_SonicC_WDvar_10min.rds"))
	bLS_STO_Sonic2_WDvar_10min <- readRDS(file = file.path(PfadRSaves, "STO_bLS_Sonic2_WDvar_10min.rds"))
	# short building
	bLS_STO_short_MK_10min <- readRDS(file = file.path(PfadRSaves, "STO_bLS_short_MK_10min.rds"))
	bLS_STO_short_QV3_10min <- readRDS(file = file.path(PfadRSaves, "STO_bLS_short_QV3_10min.rds"))
	# Source shift
	bLS_STO_Source_shift_10min <- readRDS(file = file.path(PfadRSaves, 'STO_bLS_MK_Source_shift_10min.rds'))

#############################################
### merge the different bLS runs together ###
#############################################

	# updatePath(bLS_STO_MK_10min,Cat.Path)
	# updatePath(bLS_STO_QV3_10min,Cat.Path)
	# updatePath(bLS_STO_MK_30min,Cat.Path)
	# updatePath(bLS_STO_QV3_30min,Cat.Path)
	# updatePath(bLS_STO_SonicA_WDvar_10min,Cat.Path)
	# updatePath(bLS_STO_SonicB_WDvar_10min,Cat.Path)
	# updatePath(bLS_STO_SonicC_WDvar_10min,Cat.Path)
	# updatePath(bLS_STO_Sonic2_WDvar_10min,Cat.Path)
	# updatePath(bLS_STO_short_MK_10min,Cat.Path)
	# updatePath(bLS_STO_short_QV3_10min,Cat.Path)
	# updatePath(bLS_STO_Source_shift_10min,Cat.Path)

bLS_STO_MK_10min[,Campaign := "MK"]
bLS_STO_QV3_10min[,Campaign := "QV3"]
bLS_STO_MK_30min[,Campaign := "MK"]
bLS_STO_QV3_30min[,Campaign := "QV3"]
bLS_STO_short_MK_10min[,Campaign := "MK"]
bLS_STO_short_QV3_10min[,Campaign := "QV3"]


# join to a single object
bLS_STO_orig_10min <- join(bLS_STO_MK_10min,bLS_STO_QV3_10min) # normal run
bLS_STO_orig_30min <- join(bLS_STO_MK_30min,bLS_STO_QV3_30min) # normal run
bLS_STO_WDvar_10min <- join(bLS_STO_SonicA_WDvar_10min,bLS_STO_SonicB_WDvar_10min,
	bLS_STO_SonicC_WDvar_10min,bLS_STO_Sonic2_WDvar_10min) # WDvariation -15° to 15°
bLS_STO_short_10min <- join(bLS_STO_short_MK_10min,bLS_STO_short_QV3_10min) # Building variation


bLS_STO_WDvar_10min[,Campaign := "MK"]
bLS_STO_Source_shift_10min[,Campaign := "MK"]


Conc_wo_30min <- pool(Conc_wo_10min, granularity="30mins")
Conc_P23_30min <- pool(Conc_P23_10min, granularity="30mins")
Conc_P4_30min <- pool(Conc_P4_10min, granularity="30mins")
Conc_P6_30min <- pool(Conc_P6_10min, granularity="30mins")

>-------------------------------------------------------------------------------------<
>-------------- select the correct bLS run and concentration correction --------------<
>-------------------------------------------------------------------------------------<
## select time resolution
R <- "10min"
# R <- "30min"

## select bLS run
X <- "orig" # original without any WD correciton
# X <- "WDvar" # only MK release with -15 to 15° WD correction
# X <- "short" # different building lengths
# X <- 'Source_shift'

## select concentration correction
# Y <- "wo" # without second offset correction
Y <- "P23" # using period 2 and 3
# Y <- "P6" # using period 6
>-------------------------------------------------------------------------------------<
>-------------------------------------------------------------------------------------<


bLS_STO <- copy(get(paste("bLS_STO",X,R,sep="_")))
Conc_GFs <- copy(get(paste("Conc",Y,R,sep="_")))

####################
### Calculate CQ ###
####################

bLS_STO[,c("CQ","CQ_se") := list(CE/SourceArea,CE_se/SourceArea)]

########################
### dcast to Sensors ###
########################

bLS_cast_raw <- dcast(bLS_STO, st + et + rn + Source + SourceArea + Ustar + L + Zo + Sonic +
	sUu + sVu + sWu + bw + C0 + N0 + WD + sd_WD + U_sonic +	z_canopy + T_sonic +
	Campaign  ~ Sensor, value.var = list("CE","CE_se","CQ", "CQ_se","N_TD"))

###################################
### add all remaining variables ###
###################################

# Concentraiton
# bLS_GF <- add_data(bLS_cast_raw, cbind(Conc_GFs, st=st(Conc_GFs)), on ="st")
bLS_GF <- merge(bLS_cast_raw, cbind(Conc_GFs, st=st(Conc_GFs))["18.03.2021 to "], by="st", all=TRUE)

## Weather Stations
WS1_sub <- pool(merge(WS700_1[[2]]$data[,c("WS_0170","WS_0260","WS_0360","WS_0625","WS_0960")], WS700_1[[1]]$data[,c("WS_0160","WS_0480","WD_corr")]),granularity = paste0(R,"s"), st = "21.03.2021 13:00")
names(WS1_sub) <- c("Dew_Point_WS1","Rel_Hum_WS1","Press_WS1","Precip_WS1","Rad_WS1","Temp_WS1","U_WS1","WD_WS1")
WS2_sub <- pool(merge(WS700_2[[2]]$data["01.03.2021 to ",c("WS_0170","WS_0260","WS_0360","WS_0625","WS_0960")], WS700_2[[1]]$data["01.03.2021 to ",c("WS_0160","WS_0480","WD_corr")]),granularity = paste0(R,"s"), st = "18.03.2021 10:00")
names(WS2_sub) <- c("Dew_Point_WS2","Rel_Hum_WS2","Press_WS2","Precip_WS2","Rad_WS2","Temp_WS2","U_WS2","WD_WS2")
WS700 <- merge(WS2_sub,WS1_sub,all=TRUE)
WS700_dt <- data.table(WS700["18.03.2021 to "], st=st(WS700["18.03.2021 to "]))
# bLS_GF_WS <- merge(bLS_GF, cbind(WS700, st=st(WS700))["18.03.2021 to "], by="st", all=TRUE)
bLS_GF_WS <- merge(bLS_GF, WS700_dt, by="st", all=TRUE)

## MFC
MFC_pl <- pool(MFC,granularity=paste0(R,"s"),st="09.03.2021 15:30")
# bLS_GF_WS_MFC <- add_data(bLS_GF_WS, cbind(MFC_pl, st=st(MFC_pl)), on ="st")
# bLS_GF_WS_MFC <- merge(bLS_GF_WS, cbind(MFC, st=st(MFC))["18.03.2021 to "], by="st", all=TRUE)
bLS_GF_WS_MFC <- merge(bLS_GF_WS, cbind(MFC_pl, st=st(MFC_pl))["18.03.2021 to "], by ="st", all=TRUE)

# MFC["19.03.2021 - 20.03.2021"]$Q_MFC
# bLS_GF_WS_MFC[st >= "2021-03-19" & st < "2021-03-20"]

## Drucksensor
# CC_Sensor
CC_Sensor_pl <- pool(CC_Sensor,granularity=paste0(R,"s"),st="06.03.2021 08:30")
# bLS_cast <- add_data(bLS_GF_WS_MFC, cbind(CC_Sensor_pl, st=st(CC_Sensor_pl)), on ="st")
bLS_cast <- merge(bLS_GF_WS_MFC, cbind(CC_Sensor_pl, st=st(CC_Sensor_pl))["18.03.2021 to "], by ="st", all=TRUE)
# bLS_cast <- merge(bLS_GF_WS_MFC, cbind(CC_Sensor, st=st(CC_Sensor))["18.03.2021 to "], by="st", all=TRUE)

###################################
### add some additional columns ###
###################################

## separate by NE-SW line
bLS_cast[,Wind_dir := ifelse(WD > 135 & WD <= 315,"SW","NE")]
bLS_cast[,Wind_side := c("SW"="Luv","NE"="Lee")[Wind_dir]]

## ToD
bLS_cast[, ToD := hour(st) + minute(st)/60 + 0.25]

## Gasfreisetzung
# bLS_cast[Campaign == "MK", GasMK := ]


####################
### Calcualte dC ###
####################

calcdC(bLS_cast, "GF16", "GF26",subset = bLS_cast[,Campaign == "MK"])
calcdC(bLS_cast, "GF17", "GF26",subset = bLS_cast[,Campaign == "MK"])
calcdC(bLS_cast, "GF18", "GF26",subset = bLS_cast[,Campaign == "MK"])
calcdC(bLS_cast, "GF25", "GF26",subset = bLS_cast[,Campaign == "MK"])

calcdC(bLS_cast, "GF16", "GFall",subset = bLS_cast[,Campaign == "QV3"])
calcdC(bLS_cast, "GF17", "GFall",subset = bLS_cast[,Campaign == "QV3"])
calcdC(bLS_cast, "GF18", "GFall",subset = bLS_cast[,Campaign == "QV3"])
calcdC(bLS_cast, "GF25", "GFall",subset = bLS_cast[,Campaign == "QV3"])
calcdC(bLS_cast, "GF26", "GFall",subset = bLS_cast[,Campaign == "QV3"])

####################################
### Calculate emissions [kg / h] ###
####################################

calcQ(bLS_cast, "GF16", "GF26",subset = bLS_cast[,Campaign == "MK"])
calcQ(bLS_cast, "GF17", "GF26",subset = bLS_cast[,Campaign == "MK"])
calcQ(bLS_cast, "GF18", "GF26",subset = bLS_cast[,Campaign == "MK"])
calcQ(bLS_cast, "GF25", "GF26",subset = bLS_cast[,Campaign == "MK"])

calcQ(bLS_cast, "GF16", "GFall",subset = bLS_cast[,Campaign == "QV3"])
calcQ(bLS_cast, "GF17", "GFall",subset = bLS_cast[,Campaign == "QV3"])
calcQ(bLS_cast, "GF18", "GFall",subset = bLS_cast[,Campaign == "QV3"])
calcQ(bLS_cast, "GF25", "GFall",subset = bLS_cast[,Campaign == "QV3"])
calcQ(bLS_cast, "GF26", "GFall",subset = bLS_cast[,Campaign == "QV3"])

###################################################
### save object that is now ready for filtering ###
###################################################

saveRDS(bLS_cast,file=file.path(PfadRSaves,paste0("STO_Emiss_",X,"_",Y,"_",R,".rds")))

###################################################

rm(list = ls())