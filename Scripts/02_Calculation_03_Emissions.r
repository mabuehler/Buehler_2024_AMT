
######################################
######################################
#####                            #####
#####    Emission calculation    #####
#####                            #####
######################################
######################################

# Author: Marcel Bühler
# Date: August 4, 2024
# Contact: mb@bce.au.dk or Christoph Häni christoph.haeni@bfh.ch
# Description: This script calculates emissions and makes it ready for further use.
#
# Note: This code was written by Marcel Bühler and is intended to follow the publication 'Applicability of the inverse dispersion method to measure emissions from animal housing' in AMT. 
# Please feel free to use and modify it, but attribution is appreciated.


#################
### Libraries ###
#################

library(ibts)
library(bLSmodelR)
library(RColorBrewer)


#############
### Paths ###
#############

PathData <- "Path to /data"		
PathRSaves <- "Path to /RSaves"


#################
### Functions ###
#################

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


###################################
### Colours and unit conversion ###
###################################

Col_CH4 <- "#1e88e5"
# display.brewer.all()
CH4Cols <- brewer.pal(5,"Dark2")
names(CH4Cols) <- paste0("GF",c(16:18,25:26))

mgs_to_kgd <- (24*3600)/1E6
mgs_to_kgh <- (3600)/1E6
mgm2s_to_gm2h <- 3600/1E3


#################
### load data ###
#################

# Weather station
WS700_2 <- readRDS(file=file.path(PathRSaves,"WS700.rds"))
# Concentrations
Conc_wo_10min <- readRDS(file=file.path(PathRSaves,"Conc_wo_10min.rds"))
Conc_P23_10min <- readRDS(file=file.path(PathRSaves,"Conc_P23_10min.rds"))
Conc_P4_10min <- readRDS(file=file.path(PathRSaves,"Conc_P4_10min.rds"))
Conc_P6_10min <- readRDS(file=file.path(PathRSaves,"Conc_P6_10min.rds"))
# MFC
MFC <- readRDS(file=file.path(PathRSaves,"MFC.rds"))
# Pressure sensor
CC_Sensor <- readRDS(file=file.path(PathRSaves,"CC_Sensor_raw.rds"))
### bLS results
bLS_MC_10min <- readRDS(file = file.path(PathRSaves, "bLS_MC_10min.rds"))
bLS_IC2_10min <- readRDS(file = file.path(PathRSaves, "bLS_IC2_10min.rds"))
### bLS results of wind direction variation
# bLS_SonicA_WDvar_10min <- readRDS(file = file.path(PfadRSaves, "bLS_SonicA_WDvar_10min.rds"))
# bLS_SonicB_WDvar_10min <- readRDS(file = file.path(PfadRSaves, "bLS_SonicB_WDvar_10min.rds"))
# bLS_SonicC_WDvar_10min <- readRDS(file = file.path(PfadRSaves, "bLS_SonicC_WDvar_10min.rds"))
# bLS_Sonic2_WDvar_10min <- readRDS(file = file.path(PfadRSaves, "bLS_Sonic2_WDvar_10min.rds"))


#############################################
### merge the different bLS runs together ###
#############################################

bLS_MC_10min[,Campaign := "MC"]
bLS_IC2_10min[,Campaign := "IC2"]

# join to a single object
bLS_orig_10min <- join(bLS_MC_10min,bLS_IC2_10min) # normal run

## same for WD var
# bLS_WDvar_10min <- join(bLS_SonicA_WDvar_10min,bLS_SonicB_WDvar_10min,bLS_SonicC_WDvar_10min,bLS_Sonic2_WDvar_10min) # WD variation -10° to +10°
# bLS_WDvar_10min[,Campaign := "MC"]

>-------------------------------------------------------------------------------------<
>-------------- select the correct bLS run and concentration correction --------------<
>-------------------------------------------------------------------------------------<
## select time resolution
R <- "10min" # <-------------------- the one to go!
# R <- "30min" # data not provided

## select bLS run
X <- "orig" # original without any WD correciton <-------------------- the one to go!
	# X <- "WDvar" # only MC release with -10 to 10° WD correction # bLS provided but from this point on not further used
	# X <- "short" # different building lengths # data not provided
	# X <- 'Source_shift' # data not provided

## select concentration correction
# Y <- "wo" # without second offset correction
Y <- "P23" # using period 2 and 3  <-------------------- the one to go!
# Y <- "P6" # using period 6
# Y <- "P4" # using period 4
>-------------------------------------------------------------------------------------<
>-------------------------------------------------------------------------------------<

bLS <- copy(get(paste("bLS",X,R,sep="_")))
Conc_GFs <- copy(get(paste("Conc",Y,R,sep="_")))


####################
### Calculate CQ ###
####################

bLS[,c("CQ","CQ_se") := list(CE/SourceArea,CE_se/SourceArea)]


########################
### dcast to Sensors ###
########################

bLS_cast_raw <- dcast(bLS, st + et + rn + Source + SourceArea + Ustar + L + Zo + Sonic +
	sUu + sVu + sWu + bw + C0 + N0 + WD + sd_WD + U_sonic +	z_canopy + T_sonic +
	Campaign  ~ Sensor, value.var = list("CE","CE_se","CQ", "CQ_se","N_TD"))


###################################
### add all remaining variables ###
###################################

# Concentraiton
# bLS_GF <- add_data(bLS_cast_raw, cbind(Conc_GFs, st=st(Conc_GFs)), on ="st")
bLS_GF <- merge(bLS_cast_raw, cbind(Conc_GFs, st=st(Conc_GFs))["18.03.2021 to "], by="st", all=TRUE)

## Weather Stations
WS2_sub <- pool(merge(WS700_2[[2]]$data["01.03.2021 to ",c("WS_0170","WS_0260","WS_0360","WS_0625","WS_0960")], WS700_2[[1]]$data["01.03.2021 to ",c("WS_0160","WS_0480","WD_corr")]),granularity = paste0(R,"s"), st = "18.03.2021 10:00")
names(WS2_sub) <- c("Dew_Point_WS2","Rel_Hum_WS2","Press_WS2","Precip_WS2","Rad_WS2","Temp_WS2","U_WS2","WD_WS2")
WS700 <- WS2_sub
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

## pressure sensor
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

## ToD - Time of the day resp. hour of the day
bLS_cast[, ToD := hour(st) + minute(st)/60 + 0.25]


####################
### Calcualte dC ###
####################

calcdC(bLS_cast, "GF16", "GF26",subset = bLS_cast[,Campaign == "MC"])
calcdC(bLS_cast, "GF17", "GF26",subset = bLS_cast[,Campaign == "MC"])
calcdC(bLS_cast, "GF18", "GF26",subset = bLS_cast[,Campaign == "MC"])
calcdC(bLS_cast, "GF25", "GF26",subset = bLS_cast[,Campaign == "MC"])

calcdC(bLS_cast, "GF16", "GFall",subset = bLS_cast[,Campaign == "IC2"])
calcdC(bLS_cast, "GF17", "GFall",subset = bLS_cast[,Campaign == "IC2"])
calcdC(bLS_cast, "GF18", "GFall",subset = bLS_cast[,Campaign == "IC2"])
calcdC(bLS_cast, "GF25", "GFall",subset = bLS_cast[,Campaign == "IC2"])
calcdC(bLS_cast, "GF26", "GFall",subset = bLS_cast[,Campaign == "IC2"])


####################################
### Calculate emissions [kg / h] ###
####################################

calcQ(bLS_cast, "GF16", "GF26",subset = bLS_cast[,Campaign == "MC"])
calcQ(bLS_cast, "GF17", "GF26",subset = bLS_cast[,Campaign == "MC"])
calcQ(bLS_cast, "GF18", "GF26",subset = bLS_cast[,Campaign == "MC"])
calcQ(bLS_cast, "GF25", "GF26",subset = bLS_cast[,Campaign == "MC"])

calcQ(bLS_cast, "GF16", "GFall",subset = bLS_cast[,Campaign == "IC2"])
calcQ(bLS_cast, "GF17", "GFall",subset = bLS_cast[,Campaign == "IC2"])
calcQ(bLS_cast, "GF18", "GFall",subset = bLS_cast[,Campaign == "IC2"])
calcQ(bLS_cast, "GF25", "GFall",subset = bLS_cast[,Campaign == "IC2"])
calcQ(bLS_cast, "GF26", "GFall",subset = bLS_cast[,Campaign == "IC2"])


###################################################
### save object that is now ready for filtering ###
###################################################


saveRDS(bLS_cast,file=file.path(PathRSaves,paste0("Emiss_",X,"_",Y,"_",R,".rds")))

###################################################

rm(list = ls())

