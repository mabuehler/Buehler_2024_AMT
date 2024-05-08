

################################################
################################################
#####                                      #####
#####    Fraction covered by touchdowns    #####
#####                                      #####
################################################
################################################



#################################
### Header (Pfade, Libraries) ###
#################################

	library(bLSmodelR)
	library(ibts)

	PathData <- "Path to /data"		
	PathRSaves <- "Path to /RSaves"
	Cat.Path <- paste0(PfadDaten,"/Catalogs")

#################################

### load Data

	bLS_MK <- readRDS(file=paste0(PfadRSaves,"/STO_bLS_MK.rds"))
	bLS_QV3 <- readRDS(file=paste0(PfadRSaves,"/STO_bLS_QV3.rds"))
	bLS_REL <- readRDS(file=paste0(PfadRSaves,"/STO_bLS_REL.rds"))

ncores <- 88

#####################
### calculate FCT ###
#####################

	Res_MK <- try(calc_fct(bLS_MK, ncores = ncores))

	if(inherits(Res_MK,"try-error")){
	  system("echo \"Subject: FCT MK (error exit status!!!)\" | sendmail xxxx")
	} else {
	  system("echo \"Subject: FCT MK done\" | sendmail xxxx")
	  saveRDS(Res_MK,file=paste0(PfadRSaves,"/STO_bLS_FCT_MK.rds"))
	}

	Res_QV3 <- try(calc_fct(bLS_QV3, ncores = ncores))

	if(inherits(Res_QV3,"try-error")){
	  system("echo \"Subject: FCT QV3 (error exit status!!!)\" | sendmail xxxx")
	} else {
	  system("echo \"Subject: FCT QV3 done\" | sendmail xxxx")
	  saveRDS(Res_QV3,file=paste0(PfadRSaves,"/STO_bLS_FCT_QV3.rds"))
	}

####################