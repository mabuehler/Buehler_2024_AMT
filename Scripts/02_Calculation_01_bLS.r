
########################################
########################################
#####                              #####
#####    bLS run with bLSmodelR    #####
#####                              #####
########################################
########################################

# Author: Marcel Bühler
# Date: August 2, 2024
# Contact: mb@bce.au.dk or Christoph Häni christoph.haeni@bfh.ch
# Description: This script reads in the weather station data and makes it ready for further us.
#
# Note: This code was written by Marcel Bühler and is intended to follow the publication 'Applicability of the inverse dispersion method to measure emissions from animal housing' in AMT. 
# Please feel free to use and modify it, but attribution is appreciated.

#################
### Libraries ###
#################

library(bLSmodelR)
library(ibts)

#############
### Paths ###
#############

PathData <- "Path to /data"		
PathRSaves <- "Path to /RSaves"
Cat.Path <- paste0(PathData,"/Catalogs")

#################
### Campaigns ###
#################

MC <- "18.03.2021 11:00 - 21.03.2021 14:00"
IC2	<- "21.03.2021 14:00 to 26.03.2021 12:00"

#################
### load data ###
#################

### load Sonic data and Geometry
load(file = file.path(PathRSaves,"Sonics_10min.RData"))
load(file = file.path(PathRSaves, "Geometry.RData"))


############################
############################
#####                  #####
#####    10 min bLS    #####
#####                  #####
############################
############################

##############
### bLS MC ###
##############

SonicA_MC <- SonicA[MC]
SonicB_MC <- SonicB[MC]
SonicC_MC <- SonicC[MC]
Sonic2_MC <- Sonic2[MC]


N_traj <- 1E6
ncores <- 88
MaxFetch <- 400

indSub <- " to "


# MC
Sonic_A_MC <- genInterval(
	cbind(
		setNames(as.data.frame(SonicA_MC[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(SonicA_MC[indSub])
		,et=et(SonicA_MC[indSub])
		,Sonic="SonicA"
		,as.data.frame(SonicA_MC[indSub])[,!(names(SonicA_MC) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)

Sonic_B_MC <- genInterval(
	cbind(
		setNames(as.data.frame(SonicB_MC[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(SonicB_MC[indSub])
		,et=et(SonicB_MC[indSub])
		,Sonic="SonicB"
		,as.data.frame(SonicB_MC[indSub])[,!(names(SonicB_MC) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)

Sonic_C_MC <- genInterval(
	cbind(
		setNames(as.data.frame(SonicC_MC[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(SonicC_MC[indSub])
		,et=et(SonicC_MC[indSub])
		,Sonic="SonicC"
		,as.data.frame(SonicC_MC[indSub])[,!(names(SonicC_MC) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)

Sonic_2_MC <- genInterval(
	cbind(
		setNames(as.data.frame(Sonic2_MC[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(Sonic2_MC[indSub])
		,et=et(Sonic2_MC[indSub])
		,Sonic="Sonic2"
		,as.data.frame(Sonic2_MC[indSub])[,!(names(Sonic2_MC) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)



allSonic_MC <- rbind(Sonic_A_MC,Sonic_B_MC,Sonic_C_MC,Sonic_2_MC)

## check Zo + d
ind_MC <- (allSonic_MC$"Zo [m]" + allSonic_MC$"d [m]") < 1.22
InList_MC <- genInputList(Sensors_MC[Sensors_MC[,1] %in% c('GF16','GF17','GF18','GF25','GF26'),],Source,allSonic_MC[which(ind_MC),]) # all Sonics

## run bLS model
Res_MC <- try(runbLS(InList_MC,Cat.Path,ncores = ncores))

# Number Crunsher things
if(inherits(Res_MC,"try-error")){
  system(paste0("echo \"Subject: bLS Run MC (error exit status!!!)\" | sendmail (enter your email address)"))
} else {
  system("echo \"Subject: bLS Run MC done\" | sendmail (enter your email address)")
	saveRDS(Res_MC,file=paste0(PathRSaves,"/bLS_MC_10min.rds"))
}
##


###############
### bLS IC2 ###
###############

SonicA_IC2 <- SonicA[IC2]
SonicB_IC2 <- SonicB[IC2]
SonicC_IC2 <- SonicC[IC2]
Sonic2_IC2 <- Sonic2[IC2]


N_traj <- 1E6
ncores <- 88
MaxFetch <- 400


indSub <- " to "


# IC2
Sonic_A_IC2 <- genInterval(
	cbind(
		setNames(as.data.frame(SonicA_IC2[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(SonicA_IC2[indSub])
		,et=et(SonicA_IC2[indSub])
		,Sonic="SonicA"
		,as.data.frame(SonicA_IC2[indSub])[,!(names(SonicA_IC2) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)

Sonic_B_IC2 <- genInterval(
	cbind(
		setNames(as.data.frame(SonicB_IC2[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(SonicB_IC2[indSub])
		,et=et(SonicB_IC2[indSub])
		,Sonic="SonicB"
		,as.data.frame(SonicB_IC2[indSub])[,!(names(SonicB_IC2) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)

Sonic_C_IC2 <- genInterval(
	cbind(
		setNames(as.data.frame(SonicC_IC2[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(SonicC_IC2[indSub])
		,et=et(SonicC_IC2[indSub])
		,Sonic="SonicC"
		,as.data.frame(SonicC_IC2[indSub])[,!(names(SonicC_IC2) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)

Sonic_2_IC2 <- genInterval(
	cbind(
		setNames(as.data.frame(Sonic2_IC2[indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
		,st=st(Sonic2_IC2[indSub])
		,et=et(Sonic2_IC2[indSub])
		,Sonic="Sonic2"
		,as.data.frame(Sonic2_IC2[indSub])[,!(names(Sonic2_IC2) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
	,MaxFetch=MaxFetch,N=N_traj)



allSonic_IC2 <- rbind(Sonic_A_IC2,Sonic_B_IC2,Sonic_C_IC2,Sonic_2_IC2)

### check Zo + d
ind_IC2 <- (allSonic_IC2$"Zo [m]" + allSonic_IC2$"d [m]") < 1.22
InList_IC2 <- genInputList(Sensors_IC2[Sensors_IC2[,1] %in% c('GF16','GF17','GF18','GF25','GF26'),],Source,allSonic_IC2[which(ind_IC2),])


## run bLS model
Res_IC2 <- try(runbLS(InList_IC2,Cat.Path,ncores = ncores))

# Number Crunsher things
if(inherits(Res_IC2,"try-error")){
  system(paste0("echo \"Subject: bLS Run IC2 (error exit status!!!)\" | sendmail (enter your email address)"))
} else {
  system("echo \"Subject: bLS Run  IC2 done\" | sendmail (enter your email address)")
	saveRDS(Res_IC2,file=paste0(PathRSaves,"/bLS_IC2_10min.rds"))
}






################################################
################################################
#####                                      #####
#####    WD variation from -10° to +10°    #####
#####                                      #####
################################################
################################################

# This is only here as this was part of the initial submission. We recommend not to use it for much as it is quite arbitary.

N_traj <- 1E6
# N_traj <- 2.5E5
ncores <- 88
MaxFetch <- 400

## Define CH4 release time to reduce computational time
indSub <- "19.03.2021 10:00 to 20.03.2021 08:00"

Interval_list_A <-Interval_list_B <-Interval_list_C <-Interval_list_2 <- vector(mode="list",length=31)
Sonic_list_A <-Sonic_list_B <-Sonic_list_C <-Sonic_list_2 <- vector(mode="list",length=31)
deviation_WD <- c(-10:10)
Sonic_name <- c(paste0("m",c(10:1)),0,paste0("p",c(1:10)))

### big for loop 
for(i in 1:21){

	### SonicA ###
	Sonic_list_A[[i]] <- SonicA
	Sonic_list_A[[i]][,"WD"] <- (SonicA[,"WD"] + deviation_WD[i]) %% 360

	## Release
	Interval_list_A[[i]] <- genInterval(
		cbind(
			setNames(as.data.frame(Sonic_list_A[[i]][indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
			,st=st(Sonic_list_A[[i]][indSub])
			,et=et(Sonic_list_A[[i]][indSub])
			,Sonic=paste0("SonicA_",Sonic_name[i])
			,as.data.frame(Sonic_list_A[[i]][indSub])[,!(names(Sonic_list_A[[i]]) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
		,MaxFetch=MaxFetch,N=N_traj)

	### SonicB ###
	Sonic_list_B[[i]] <- SonicB
	Sonic_list_B[[i]][,"WD"] <- (SonicB[,"WD"] + deviation_WD[i]) %% 360

	## Release
	Interval_list_B[[i]] <- genInterval(
		cbind(
			setNames(as.data.frame(Sonic_list_B[[i]][indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
			,st=st(Sonic_list_B[[i]][indSub])
			,et=et(Sonic_list_B[[i]][indSub])
			,Sonic=paste0("SonicB_",Sonic_name[i])
			,as.data.frame(Sonic_list_B[[i]][indSub])[,!(names(Sonic_list_B[[i]]) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
		,MaxFetch=MaxFetch,N=N_traj)

	### SonicC ###
	Sonic_list_C[[i]] <- SonicC
	Sonic_list_C[[i]][,"WD"] <- (SonicC[,"WD"] + deviation_WD[i]) %% 360

	## Release
	Interval_list_C[[i]] <- genInterval(
		cbind(
			setNames(as.data.frame(Sonic_list_C[[i]][indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
			,st=st(Sonic_list_C[[i]][indSub])
			,et=et(Sonic_list_C[[i]][indSub])
			,Sonic=paste0("SonicC_",Sonic_name[i])
			,as.data.frame(Sonic_list_C[[i]][indSub])[,!(names(Sonic_list_C[[i]]) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
		,MaxFetch=MaxFetch,N=N_traj)

	### Sonic2 ###
	Sonic_list_2[[i]] <- Sonic2
	Sonic_list_2[[i]][,"WD"] <- (Sonic2[,"WD"] + deviation_WD[i]) %% 360

	## Release
	Interval_list_2[[i]] <- genInterval(
		cbind(
			setNames(as.data.frame(Sonic_list_2[[i]][indSub])[,c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic")],c("Ustar","L","Zo","WD","sd_WD","d","sUu","sVu","sWu","z_sWu"))
			,st=st(Sonic_list_2[[i]][indSub])
			,et=et(Sonic_list_2[[i]][indSub])
			,Sonic=paste0("Sonic2_",Sonic_name[i])
			,as.data.frame(Sonic_list_2[[i]][indSub])[,!(names(Sonic_list_2[[i]]) %in% c("Ustar","L","z0","WD","sd_WD","d","sUu","sVu","sWu","z_sonic"))])
		,MaxFetch=MaxFetch,N=N_traj)

}

## add together
SonicA_WDvar <- rbindlist(Interval_list_A)
SonicB_WDvar <- rbindlist(Interval_list_B)
SonicC_WDvar <- rbindlist(Interval_list_C)
Sonic2_WDvar <- rbindlist(Interval_list_2)
# assign again the correct 'classes' as it will otherwise not work ('Interval' was replaced by 'data.table')
class(SonicA_WDvar) <- c("Interval","data.frame") 
class(SonicB_WDvar) <- c("Interval","data.frame")
class(SonicC_WDvar) <- c("Interval","data.frame")
class(Sonic2_WDvar) <- c("Interval","data.frame")

## check Zo + d
ind_A_WDvar <- (SonicA_WDvar$"Zo [m]" + SonicA_WDvar$"d [m]") < 1.22
InList_A_WDvar <- genInputList(Sensors_MC[Sensors_MC[,1] %chin% c("GF16","GF17","GF18","GF25"),],Source,SonicA_WDvar[which(ind_A_WDvar),]) ## ohne GF26
#
ind_B_WDvar <- (SonicB_WDvar$"Zo [m]" + SonicB_WDvar$"d [m]") < 1.22
InList_B_WDvar <- genInputList(Sensors_MC[Sensors_MC[,1] %chin% c("GF16","GF17","GF18","GF25"),],Source,SonicB_WDvar[which(ind_B_WDvar),]) ## ohne GF26
#
ind_C_WDvar <- (SonicC_WDvar$"Zo [m]" + SonicC_WDvar$"d [m]") < 1.22
InList_C_WDvar <- genInputList(Sensors_MC[Sensors_MC[,1] %chin% c("GF16","GF17","GF18","GF25"),],Source,SonicC_WDvar[which(ind_C_WDvar),]) ## ohne GF26
#
ind_2_WDvar <- (Sonic2_WDvar$"Zo [m]" + Sonic2_WDvar$"d [m]") < 1.22
InList_2_WDvar <- genInputList(Sensors_MC[Sensors_MC[,1] %chin% c("GF16","GF17","GF18","GF25"),],Source,Sonic2_WDvar[which(ind_2_WDvar),]) ## ohne GF26
#


## run bLS model
Res_A_WDvar_10deg <- try(runbLS(InList_A_WDvar,Cat.Path,ncores = ncores))

# Number Crunsher things
if(inherits(Res_A_WDvar_10deg,"try-error")){
  system(paste0("echo \"Subject: bLS Run A_WDvar (error exit status!!!)\" | sendmail youremailaddress"))
} else {
  system("echo \"Subject: bLS Run A_WDvar done\" | sendmail youremailaddress")
	saveRDS(Res_A_WDvar_10deg,file=paste0(PathRSaves,"/bLS_SonicA_WDvar_10deg_10min.rds"))
}
#### run bLS model
Res_B_WDvar_10deg <- try(runbLS(InList_B_WDvar,Cat.Path,ncores = ncores))

# Number Crunsher things
if(inherits(Res_B_WDvar_10deg,"try-error")){
  system(paste0("echo \"Subject: bLS Run B_WDvar_10deg (error exit status!!!)\" | sendmail youremailaddress"))
} else {
  system("echo \"Subject: bLS Run B_WDvar_10deg done\" | sendmail youremailaddress")
	saveRDS(Res_B_WDvar_10deg,file=paste0(PathRSaves,"/bLS_SonicB_WDvar_10deg_10min.rds"))
}
##

#### run bLS model
Res_C_WDvar_10deg <- try(runbLS(InList_C_WDvar,Cat.Path,ncores = ncores))

# Number Crunsher things
if(inherits(Res_C_WDvar_10deg,"try-error")){
  system(paste0("echo \"Subject: bLS Run C_WDvar_10deg (error exit status!!!)\" | sendmail youremailaddress"))
} else {
  system("echo \"Subject: bLS Run C_WDvar_10deg done\" | sendmail youremailaddress")
	saveRDS(Res_C_WDvar_10deg,file=paste0(PathRSaves,"/bLS_SonicC_WDvar_10deg_10min.rds"))
}
##

#### run bLS model
Res_2_WDvar_10deg <- try(runbLS(InList_2_WDvar,Cat.Path,ncores = ncores))

# Number Crunsher things
if(inherits(Res_2_WDvar_10deg,"try-error")){
  system(paste0("echo \"Subject: bLS Run 2_WDvar_10deg (error exit status!!!)\" | sendmail youremailaddress"))
} else {
  system("echo \"Subject: bLS Run 2_WDvar_10deg done\" | sendmail marcel.buehler@bfh.ch")
	saveRDS(Res_2_WDvar_10deg,file=paste0(PathRSaves,"/bLS_Sonic2_WDvar_10deg_10min.rds"))
}
##



