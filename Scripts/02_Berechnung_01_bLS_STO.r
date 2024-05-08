

######################################################
######################################################
#####                                            #####
#####    bLS Runs Quellenversuch Stockerematt    #####
#####                                            #####
######################################################
######################################################


#################################
### Header (Pfade, Libraries) ###
#################################

	library(bLSmodelR)
	library(ibts)

PathData <- "Path to /data"		
PathRSaves <- "Path to /RSaves"
Cat.Path <- paste0(PfadDaten,"/Catalogs")

#################################


### load Sonic data and Geometry
load(file = file.path(PfadRSaves,"STO_Sonics_10min.RData"))
load(file = file.path(PfadRSaves,"STO_Sonics_30min.RData"))
load(file = file.path(PfadRSaves, "Geometry_STO.RData"))


	### MKs
	QV1 <- "20.01.2021 to 28.01.2021 12:00"
	QV2 <- "05.03.2021 to 10.03.2021 18:00"
	QV1u2 <- "20.01.2021 to 10.03.2021 18:00"
	MK <- "18.03.2021 11:00 - 21.03.2021 14:00"
	MK_rel <- "19.03.2021 10:00 to 20.03.2021 08:00"
	QV3	<- "21.03.2021 14:00 to 26.03.2021 12:00"
	MK_QV3 <- "18.03.2021 11:00 to 26.03.2021 12:00" 


############################
############################
#####                  #####
#####    10 min bLS    #####
#####                  #####
############################
############################

##############
### bLS MK ###
##############

SonicA_MK <- SonicA[MK]
SonicB_MK <- SonicB[MK]
SonicC_MK <- SonicC[MK]
Sonic2_MK <- Sonic2[MK]


N_traj <- 1E6
# N_traj <- 2.5E5
ncores <- 88
MaxFetch <- 400

indSub <- " to "


# MK
Sonic_A_MK <- genInterval(
	cbind(
		setNames(as.data.frame(SonicA_MK[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(SonicA_MK[indSub])
		,et=et(SonicA_MK[indSub])
		,Sonic="SonicA"
		,as.data.frame(SonicA_MK[indSub])[,!(names(SonicA_MK) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)

Sonic_B_MK <- genInterval(
	cbind(
		setNames(as.data.frame(SonicB_MK[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(SonicB_MK[indSub])
		,et=et(SonicB_MK[indSub])
		,Sonic="SonicB"
		,as.data.frame(SonicB_MK[indSub])[,!(names(SonicB_MK) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)

Sonic_C_MK <- genInterval(
	cbind(
		setNames(as.data.frame(SonicC_MK[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(SonicC_MK[indSub])
		,et=et(SonicC_MK[indSub])
		,Sonic="SonicC"
		,as.data.frame(SonicC_MK[indSub])[,!(names(SonicC_MK) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)

Sonic_2_MK <- genInterval(
	cbind(
		setNames(as.data.frame(Sonic2_MK[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(Sonic2_MK[indSub])
		,et=et(Sonic2_MK[indSub])
		,Sonic="Sonic2"
		,as.data.frame(Sonic2_MK[indSub])[,!(names(Sonic2_MK) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)



allSonic_MK <- rbind(Sonic_A_MK,Sonic_B_MK,Sonic_C_MK,Sonic_2_MK)

## check Zo + d
ind_MK <- (allSonic_MK$"Zo [m]" + allSonic_MK$"d [m]") < 1.22
InList_MK <- genInputList(Sensors_MK,Sources,allSonic_MK[which(ind_MK),]) # all Sonics

## run bLS model
Res_MK <- try(runbLS(InList_MK,Cat.Path,ncores = ncores))

# Number Crunsher things
if(inherits(Res_MK,"try-error")){
  system(paste0("echo \"Subject: bLS Run STO MK (error exit status!!!)\" | sendmail (enter your email address)"))
} else {
  system("echo \"Subject: bLS Run STO MK done\" | sendmail (enter your email address)")
	saveRDS(Res_MK,file=paste0(PfadRSaves,"/STO_bLS_MK_10min.rds"))
}
##


###############
### bLS QV3 ###
###############

SonicA_QV3 <- SonicA[QV3]
SonicB_QV3 <- SonicB[QV3]
SonicC_QV3 <- SonicC[QV3]
Sonic2_QV3 <- Sonic2[QV3]


N_traj <- 1E6
# N_traj <- 2.5E5
ncores <- 88
MaxFetch <- 400


indSub <- " to "


# QV3
Sonic_A_QV3 <- genInterval(
	cbind(
		setNames(as.data.frame(SonicA_QV3[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(SonicA_QV3[indSub])
		,et=et(SonicA_QV3[indSub])
		,Sonic="SonicA"
		,as.data.frame(SonicA_QV3[indSub])[,!(names(SonicA_QV3) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)

Sonic_B_QV3 <- genInterval(
	cbind(
		setNames(as.data.frame(SonicB_QV3[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(SonicB_QV3[indSub])
		,et=et(SonicB_QV3[indSub])
		,Sonic="SonicB"
		,as.data.frame(SonicB_QV3[indSub])[,!(names(SonicB_QV3) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)

Sonic_C_QV3 <- genInterval(
	cbind(
		setNames(as.data.frame(SonicC_QV3[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(SonicC_QV3[indSub])
		,et=et(SonicC_QV3[indSub])
		,Sonic="SonicC"
		,as.data.frame(SonicC_QV3[indSub])[,!(names(SonicC_QV3) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)

Sonic_2_QV3 <- genInterval(
	cbind(
		setNames(as.data.frame(Sonic2_QV3[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(Sonic2_QV3[indSub])
		,et=et(Sonic2_QV3[indSub])
		,Sonic="Sonic2"
		,as.data.frame(Sonic2_QV3[indSub])[,!(names(Sonic2_QV3) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)



allSonic_QV3 <- rbind(Sonic_A_QV3,Sonic_B_QV3,Sonic_C_QV3,Sonic_2_QV3)

### check Zo + d
ind_QV3 <- (allSonic_QV3$"Zo [m]" + allSonic_QV3$"d [m]") < 1.22
InList_QV3 <- genInputList(Sensors_QV3,Sources,allSonic_QV3[which(ind_QV3),]) # SonicC und Sonic2


## run bLS model
Res_QV3 <- try(runbLS(InList_QV3,Cat.Path,ncores = ncores))

# Number Crunsher things
if(inherits(Res_QV3,"try-error")){
  system(paste0("echo \"Subject: bLS Run STO QV3 (error exit status!!!)\" | sendmail (enter your email address)"))
} else {
  system("echo \"Subject: bLS Run STO QV3 done\" | sendmail (enter your email address)")
	saveRDS(Res_QV3,file=paste0(PfadRSaves,"/STO_bLS_QV3_10min.rds"))
}

#########
######



###############
### bLS QV3 ###
###############

SonicA_QV3 <- SonicA[QV3]
SonicB_QV3 <- SonicB[QV3]
SonicC_QV3 <- SonicC[QV3]
Sonic2_QV3 <- Sonic2[QV3]


N_traj <- 1E6
# N_traj <- 2.5E5
ncores <- 88
MaxFetch <- 400


indSub <- " to "


# QV3
Sonic_A_QV3 <- genInterval(
	cbind(
		setNames(as.data.frame(SonicA_QV3[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(SonicA_QV3[indSub])
		,et=et(SonicA_QV3[indSub])
		,Sonic="SonicA"
		,as.data.frame(SonicA_QV3[indSub])[,!(names(SonicA_QV3) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)

Sonic_B_QV3 <- genInterval(
	cbind(
		setNames(as.data.frame(SonicB_QV3[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(SonicB_QV3[indSub])
		,et=et(SonicB_QV3[indSub])
		,Sonic="SonicB"
		,as.data.frame(SonicB_QV3[indSub])[,!(names(SonicB_QV3) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)

Sonic_C_QV3 <- genInterval(
	cbind(
		setNames(as.data.frame(SonicC_QV3[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(SonicC_QV3[indSub])
		,et=et(SonicC_QV3[indSub])
		,Sonic="SonicC"
		,as.data.frame(SonicC_QV3[indSub])[,!(names(SonicC_QV3) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)

Sonic_2_QV3 <- genInterval(
	cbind(
		setNames(as.data.frame(Sonic2_QV3[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(Sonic2_QV3[indSub])
		,et=et(Sonic2_QV3[indSub])
		,Sonic="Sonic2"
		,as.data.frame(Sonic2_QV3[indSub])[,!(names(Sonic2_QV3) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)



allSonic_QV3 <- rbind(Sonic_A_QV3,Sonic_B_QV3,Sonic_C_QV3,Sonic_2_QV3)

### check Zo + d
ind_QV3 <- (allSonic_QV3$"Zo [m]" + allSonic_QV3$"d [m]") < 1.22
InList_QV3 <- genInputList(Sensors_QV3,Sources_short,allSonic_QV3[which(ind_QV3),]) # SonicC und Sonic2


## run bLS model
Res_QV3 <- try(runbLS(InList_QV3,Cat.Path,ncores = ncores))

# Number Crunsher things
if(inherits(Res_QV3,"try-error")){
  system(paste0("echo \"Subject: bLS Run STO QV3 (error exit status!!!)\" | sendmail (enter your email address)"))
} else {
  system("echo \"Subject: bLS Run STO QV3 done\" | sendmail (enter your email address)")
	saveRDS(Res_QV3,file=paste0(PfadRSaves,"/STO_bLS_short_QV3_10min.rds"))
}

#########
######