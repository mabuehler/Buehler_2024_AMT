


###########################################################
###########################################################
#####                                                 #####
#####    Data treatment ultrasonic anemometer data    #####
#####                                                 #####
###########################################################
###########################################################




###########################
### libraries and Paths ###
###########################

library(ibts)
library(openxlsx)
library(data.table)
library(readxl)
library(RgoogleMaps)
library(bLSmodelR)
library(ggplot2)

PathData <- "Path to /data"   
PathRSaves <- "Path to /RSaves"
PathFigures <- "~/repos/4_Projects/8_FerLAS/Figures"
source("https://raw.githubusercontent.com/hafl-gel/gel-scripts/main/sonic-turbulence.r")
source("https://raw.githubusercontent.com/hafl-gel/gel-scripts/main/wgs84-ch1903.r")



###################
### time offset ###
###################

time_offset_raw <- read_excel(file.path(dirname(PathRSaves),"Other/time_offsets.xlsx"))
time_offset <- as.data.table(time_offset_raw)
setnames(time_offset,c("Inst","RefCET","GF"))
time_offset[,RefTime := parse_date_time3(RefCET,tz="CET")]
time_offset[,GFTime := parse_date_time3(GF,tz="Etc/GMT-1")]
setkey(time_offset, Inst)
#
ZvSonicA <- time_offset["SonicA", mean(RefTime - GFTime)]
ZvSonicB <- time_offset["SonicB", mean(RefTime - GFTime)]
ZvSonicC <- time_offset["SonicC", mean(RefTime - GFTime)]
ZvSonic2 <- time_offset["Sonic2", mean(RefTime - GFTime)]


 ### Funktionen
shift_dt <- function(x,ibts,tz="Etc/GMT-1",Plot=FALSE){
	versatz <- as.numeric(x$RefTime - x$GFTime,units="secs")
	zeiten <- with_tz(x$RefTime,tz=tz)
	a <- versatz[-length(versatz)]
	b <- (versatz[-1] - a)/as.numeric(zeiten[-1] - zeiten[-length(zeiten)],units="secs")
	st_out <- st_in <- st(ibts)
	et_out <- et_in <- et(ibts)
	ind <- findInterval(st_in,zeiten,all.inside=TRUE)
	for(i in unique(ind)){
	  st_sub <- st_in[ind == i]
	  # browser()
	  st_out[ind == i] <- st_sub + round(a[i] + b[i]*as.numeric(st_sub - zeiten[i],units="secs"), 5)
	  et_sub <- et_in[ind == i]
	  et_out[ind == i] <- et_sub + round(a[i] + b[i]*as.numeric(et_sub - zeiten[i],units="secs"), 5)
	}
	if(Plot){
	  par(mfrow=c(2,1))
	  d_st <- st_out-st_in
	  d_et <- et_out-et_in
	  plot(ibts[,1],blank=TRUE,ylim=range(d_st),ylab=attr(d_st,"units"),main="st")
	  lines(st_in,d_st)
	  plot(ibts[,1],blank=TRUE,ylim=range(d_et),ylab=attr(d_et,"units"),main="et")
	  lines(et_in,d_et)
	}
	attr(ibts,"st") <- st_out
	attr(ibts,"et") <- et_out
	ibts
	 }

lines_sec2xy <- function(xyMK,sensor,node=1,wd,col="lightblue",lwd=2,...){
  GF <- xyMK[xyMK[,1] %in% sensor,]
  sens <- as.numeric(GF[GF[,3] == node,4:5])
  b <- tan((90 - wd)/180*pi)
  x <- if(wd <= 180) 600 else -600
  y <- sens[2] - (sens[1] - x)*b
  lines(c(sens[1],x),c(sens[2],y),col=col,lwd=lwd,...)
}

#####################
### canopy height ###
#####################

# uniform canopy height of 8 cm.  For Sonic C 20 cm
canopy <- 0.08
canopy_c <- 0.20


##############################
### Weather station offset ###
##############################

WS700_2 <- readRDS(file.path(PathRSaves, "WS700.rds"))
# 5min Offset Korrektur (noch nicht angewendet)
WS <- shift_time(WS700_2[[1]]$data,"-5mins") 
# remove no wind data
WS[which(WS$WS_0460 == 0),c("WS_0460","WD_corr")] <- NA_real_
# WS[which(WS$WS_0480 == 0),c("WS_0460","WD_corr")] <- NA_real_

  # Pfadlängen & Geometrie
  load(file.path(PathRSaves,"Coordinates_bLSmodelR.RData"))
  # Google map
  STO_Map <- ReadMapTile(file.path(PathFigures,"STO_GoogleMaps.png"))
  Sensors_MK_xy <- CH.to.map(STO_Map,Sensors_MK)
  Sensors_QV3_xy <- CH.to.map(STO_Map,Sensors_QV3)
  Sources_xy <- CH.to.map(STO_Map,Sources)
  Sonics_MK_xy <- CH.to.map(STO_Map,Sonics_MK)
  Sonics_QV3_xy <- CH.to.map(STO_Map,Sonics_QV3)
  WeatherStation_xy <- CH.to.map(STO_Map,WeatherStation)


#########################################################################
### Determine deviaiton to north with the help of the weather station ###
#########################################################################

  ###############
  ### Sonic A ###
  ###############

  sA_dir <- dir(file.path(PathData,"Sonic","SonicA"), full.names = TRUE)
  sA_list <- lapply(sA_dir[], readWindMaster_ascii)
  sA_data <- rbindlist(sA_list)
  d_t <- sA_data[, as.numeric(diff(Time)) / 2]
  md_t <- median(d_t)
  d_t[d_t > md_t*2] <- md_t 
  sA_ibts <- shift_dt(time_offset["SonicA"],as.ibts(sA_data[,-(1:2)], st = sA_data[, Time] - c(md_t, d_t), 
    et = sA_data[, Time] + c(d_t, md_t)))
  sA_ws <- sA_ibts[WS[, c("WS_0460", "WD_corr")]]
  sA_ws[,c("U_sA","WD_sA")] <- with(sA_ws,{
    out <- mapply(rotate_twoaxis,u=u, v=v, w=w, T=T, MoreArgs = list(c.system = "windmaster"))
    # browser()
    list(
      U = unlist(out["uprot",])
      ,WD = unlist(out["wd",])
      )
  })

  # Abweichung zu Nord
  sA_off <- 0
  sA_ws$WD_sAcorr <- (sA_ws$WD_sA + sA_off) %% 360
  colClasses(sA_ws)[c("WD_sA","WD_sAcorr")] <- "circ"
  sA_ws30 <- pool(sA_ws, "1mins")
  sA_ws30$WS_m_sA <- sA_ws30$WD_corr - sA_ws30$WD_sAcorr
  sA_ws30$sA_m_WS <- sA_ws30$WD_sAcorr - sA_ws30$WD_corr
  sA_ws30$WDdiff <- pmin(
    sA_ws30$WS_m_sA %% 360
    ,sA_ws30$sA_m_WS %% 360
    )
  sA_ws30$WDdiff[which(sA_ws30$WDdiff == sA_ws30$sA_m_WS)] <- -sA_ws30$WDdiff[which(sA_ws30$WDdiff == sA_ws30$sA_m_WS)]
  m_sA_1 <- MASS::rlm(WDdiff ~ 1, sA_ws30["21.01.2021 16:00 - 27.01.2021 13:00"], subset = U_sA > 0.8 & WDdiff > -69-30 & WDdiff < 69+30)
  m_sA_2 <- MASS::rlm(WDdiff ~ 1, sA_ws30["05.03.2021 15:00 - 10.03.2021 16:00"], subset = U_sA > 0.8 & WDdiff > 12-30 & WDdiff < 12+30)
  m_sA_3 <- MASS::rlm(WDdiff ~ 1, sA_ws30["19.03.2021 09:30 - 21.03.2021 12:00"], subset = U_sA > 0.8 & (WD_corr < 38 | WD_corr > 75) & WDdiff < 2+30 & WDdiff > 2-30) 
  m_sA_4 <- MASS::rlm(WDdiff ~ 1, sA_ws30["21.03.2021 15:10 - 26.03.2021 11:00"], subset = U_sA > 0.8 & (WD_corr < 38 | WD_corr > 75) & abs(WDdiff) < 22+30)
  wdoff_sA_1 <- (sA_off + coef(m_sA_1)) %% 360
  wdoff_sA_2 <- (sA_off + coef(m_sA_2)) %% 360
  wdoff_sA_3 <- (sA_off + coef(m_sA_3)) %% 360
  wdoff_sA_4 <- (sA_off + coef(m_sA_4)) %% 360
  # wdoff_sA_1 <- 294.1486
  # wdoff_sA_2 <- 13.5 
  # wdoff_sA_3 <- 359.4507 # habe das bLS nicht neu gerechnet, da nur 0.5 grad
  # wdoff_sA_4 <- 20.9562 # did not adapt this with the 1min data

  A1 <- as.data.table(sA_ws30["21.01.2021 16:00 - 27.01.2021 13:00"])
  A2 <- as.data.table(sA_ws30["05.03.2021 15:00 - 10.03.2021 16:00"])
  A3 <- as.data.table(sA_ws30["19.03.2021 09:30 - 21.03.2021 12:00"])
  A4 <- as.data.table(sA_ws30["21.03.2021 15:10 - 26.03.2021 11:00"])

  A1[,summary(WDdiff)]
  coef(MASS::rlm(WDdiff ~ 1, sA_ws30["21.01.2021 16:00 - 27.01.2021 13:00"]))
  plot(A1[,(WDdiff+180) %% 360 - 180])
  plot(A1[WS_0460 > 0.8,(WDdiff+180) %% 360 - 180])
  A1[WS_0460 > 0.8,summary(WDdiff)]
  plot(A1[WDdiff > -68-40 & WDdiff < -68+40,(WDdiff+180) %% 360 - 180])
  A1[WS_0460 > 0.8 & WDdiff > -68-40 & WDdiff < -68+40,summary(WDdiff)]

  A2[,summary((WDdiff+180) %% 360 - 180)]
  coef(MASS::rlm(WDdiff ~ 1, sA_ws30["05.03.2021 15:00 - 10.03.2021 16:00"]))
  plot(A2[,(WDdiff+180) %% 360 - 180])
  plot(A2[WS_0460 > 0.8,(WDdiff+180) %% 360 - 180])
  A2[WS_0460 > 0.8,summary(WDdiff)]
  plot(A2[WDdiff > 14-40 & WDdiff < 14+40,(WDdiff+180) %% 360 - 180])
  A2[WS_0460 > 0.8 & WDdiff > 14-40 & WDdiff < 14+40,summary(WDdiff)]

  A3[,summary((WDdiff+180) %% 360 - 180)]
  plot(WD_sAcorr ~ WD_corr,A3)
  points(WD_sAcorr ~ WD_corr,A3[U_sA < 0.8],col="red",pch=19)
  points(WD_sAcorr ~ WD_corr,A3[WD_corr > 38 & WD_corr < 75],col="red",pch=19)
  abline(v=c(38,75))
  abline(5,1)
  plot(A3[,(WDdiff+180) %% 360 - 180])
  plot(A3[WS_0460 > 0.8,(WDdiff+180) %% 360 - 180])
  A3[WS_0460 > 0.8,summary(WDdiff)]
  plot(A3[WDdiff > 2-40 & WDdiff < 2+40,(WDdiff+180) %% 360 - 180])
  plot(WDdiff~I((WD_corr+180)%% 360 - 180),A3)
  abline(v=c(38,75))
  A3[WS_0460 > 0.8 & WDdiff > -2-40 & WDdiff < -2+40,summary(WDdiff)]

  A4[,summary(WDdiff)]
  plot(WD_sAcorr ~ WD_corr,A4)
  points(WD_sAcorr ~ WD_corr,A4[U_sA < 0.5],col="red",pch=19)
  points(WD_sAcorr ~ WD_corr,A4[WD_corr > 38 & WD_corr < 62],col="red",pch=19)
  points(WD_sAcorr ~ WD_corr,A4[abs(WDdiff) > 22+30],col="red",pch=19)
  abline(v=c(38,62))


  iWDSaMK <- c(38,75)
  iWDSaQV3 <- c(38,75)
  PlotOnStaticMap(STO_Map)
  plot(Sonics_MK_xy[Sonics_MK_xy[,1] %in% c("SonicA"),],sensor.text.args=list(labels="",cex=2),points.args = list(pch = 20, cex = 2, col ="orange"),add=TRUE)
  plot(Sonics_QV3_xy[Sonics_QV3_xy[,1] %in% c("SonicA"),],sensor.text.args=list(labels="",cex=2),points.args = list(pch = 20, cex = 2, col ="lightblue"),add=TRUE)
  plot(WeatherStation_xy[WeatherStation_xy[,1] %in% c("WS2"),],sensor.text.args=list(labels="",cex=2),points.args = list(pch = 20, cex = 2, col ="red"),add=TRUE)
  lines_sec2xy(Sonics_MK_xy,"SonicA",1,iWDSaMK[1],lty=2, col="orange")
  lines_sec2xy(Sonics_MK_xy,"SonicA",1,iWDSaMK[2],lty=2, col="orange")
  lines_sec2xy(WeatherStation_xy,"WS2",1,iWDSaMK[1],lty=2, col="orange")
  lines_sec2xy(WeatherStation_xy,"WS2",1,iWDSaMK[2],lty=2, col="orange")
  # QV3
  lines_sec2xy(Sonics_QV3_xy,"SonicA",1,iWDSaQV3[1],lty=2)
  lines_sec2xy(Sonics_QV3_xy,"SonicA",1,iWDSaQV3[2],lty=2)
  lines_sec2xy(WeatherStation_xy,"WS2",1,iWDSaQV3[1],lty=2)
  lines_sec2xy(WeatherStation_xy,"WS2",1,iWDSaQV3[2],lty=2)


  ###############
  ### Sonic B ###
  ###############

  sB_dir <- dir(file.path(PathData,"Sonic","SonicB"), full.names = TRUE)
  sB_list <- lapply(sB_dir[], readWindMaster_ascii)
  sB_data <- rbindlist(sB_list)
  d_t <- sB_data[, as.numeric(diff(Time)) / 2]
  md_t <- median(d_t)
  d_t[d_t > md_t*2] <- md_t
### Es hat ein paar Intervalle, die nicht gehen. die müssen ausgeschlossen werden. Man kann aber das ibts Objekt icht printen (auch von SonicA und SonicC nicht)
# shift_dt(time_offset["SonicB"],as.ibts(sB_data[,-(1:2)], st = sB_data[, Time] - c(md_t, d_t), et = sB_data[, Time] + c(d_t, md_t)))
# print(as.numeric(sB_data[13796:13810,Time])- 1611159614.174289 ,digit = 5) 
# X <- sB_data[-(13797:13804), Time] - c(md_t, d_t[-(13797:13804)])
# Y <- sB_data[-(13797:13804), Time] + c(d_t[-(13797:13804)], md_t)
# min(as.numeric(diff(X)))
# min(as.numeric(diff(Y)))
  sB_ibts <- shift_dt(time_offset["SonicB"],as.ibts(sB_data[-(13797:13804),-(1:2)], st = sB_data[-(13797:13804), Time] - c(md_t, d_t[-(13797:13804)]),
      et = sB_data[-(13797:13804), Time] + c(d_t[-(13797:13804)], md_t)))
  sB_ws <- sB_ibts[WS[, c("WS_0460", "WD_corr")]]
  sB_ws[,c("U_sB","WD_sB")] <- with(sB_ws,{
    out <- mapply(rotate_twoaxis,u=u, v=v, w=w, T=T, MoreArgs = list(c.system = "windmaster"))
    # browser()
    list(
      U = unlist(out["uprot",])
      ,WD = unlist(out["wd",])
      )
  })

  # Abweichung zu Nord
  sB_off <- 0
  sB_ws$WD_sBcorr <- (sB_ws$WD_sB + sB_off) %% 360
  colClasses(sB_ws)[c("WD_sB","WD_sBcorr")] <- "circ"
  sB_ws30 <- pool(sB_ws, "1mins")
  sB_ws30$WS_m_sB <- sB_ws30$WD_corr - sB_ws30$WD_sBcorr
  sB_ws30$sB_m_WS <- sB_ws30$WD_sBcorr - sB_ws30$WD_corr
  sB_ws30$WDdiff <- pmin(
    sB_ws30$WS_m_sB %% 360
    ,sB_ws30$sB_m_WS %% 360
    )
  sB_ws30$WDdiff[which(sB_ws30$WDdiff == sB_ws30$sB_m_WS)] <- -sB_ws30$WDdiff[which(sB_ws30$WDdiff == sB_ws30$sB_m_WS)]
  m_sB_1 <- MASS::rlm(WDdiff ~ 1, sB_ws30["21.01.2021 15:00 - 27.01.2021 14:00"], subset = U_sB > 0.8 & abs(WDdiff) > 122-40 & abs(WDdiff) < 122+40)
  m_sB_2 <- MASS::rlm(WDdiff ~ 1, sB_ws30["18.03.2021 15:30 - 22.03.2021 00:00"], subset = U_sB > 0.8 & (WD_corr < 42 | WD_corr > 75) & abs(WDdiff) < 4+30)
  m_sB_3 <- MASS::rlm(WDdiff ~ 1, sB_ws30["22.03.2021 - 26.03.2021 16:00"], subset = U_sB > 0.5 & (WD_corr < 40 | WD_corr > 75) & abs(WDdiff) < 53+40)
  wdoff_sB_1 <- (sB_off + coef(m_sB_1)) %% 360
  wdoff_sB_2 <- (sB_off + coef(m_sB_2)) %% 360
  wdoff_sB_3 <- (sB_off + coef(m_sB_3)) %% 360
  # wdoff_sB_1 <- 123.3213 
  # wdoff_sB_2 <- 0.0974 
  # wdoff_sB_3 <- 312.44

  B1 <- as.data.table(sB_ws30["21.01.2021 15:00 - 27.01.2021 14:00"])
  B2 <- as.data.table(sB_ws30["18.03.2021 15:30 - 22.03.2021 00:00"])
  B3 <- as.data.table(sB_ws30["22.03.2021 - 26.03.2021 16:00"])

  B1[,summary(WDdiff)]
  plot(WD_sBcorr ~ WD_corr,B1)
  points(WD_sBcorr ~ WD_corr,B1[U_sB < 0.8],col="red",pch=19)
  points(WD_sBcorr ~ WD_corr,B1[WDdiff > 123+40 | WDdiff < 123-40],col="red",pch=19)

  B2[,summary(abs(WDdiff))]
  plot(WD_sBcorr ~ WD_corr,B2)
  points(WD_sBcorr ~ WD_corr,B2[U_sB < 0.8],col="red",pch=19)
  points(WD_sBcorr ~ WD_corr,B2[WD_corr > 42 & WD_corr < 75],col="red",pch=19)
  points(WD_sBcorr ~ WD_corr,B2[abs(WDdiff) > 4+30],col="red",pch=19)
  abline(v=c(42,75))
  abline(3,1)
  plot(WDdiff~I((WD_corr+180)%% 360 - 180),B2)
  plot(WDdiff~I((WD_corr+180)%% 360 - 180),B2,xlim=c(-40,120))
  abline(v=c(38,75))

  B3[,summary(abs(WDdiff))]
  plot(WD_sBcorr ~ WD_corr,B3)
  points(WD_sBcorr ~ WD_corr,B3[U_sB < 0.5],col="red",pch=19)
  points(WD_sBcorr ~ WD_corr,B3[WD_corr > 40 & WD_corr < 69],col="red",pch=19)
  points(WD_sBcorr ~ WD_corr,B3[abs(WDdiff) > 53+40 | abs(WDdiff) < 53-40],col="red",pch=19)


  iWDSbMK <- c(42,75)
  iWDSbQV3 <- c(40,75)
  PlotOnStaticMap(STO_Map)
  plot(Sonics_MK_xy[Sonics_MK_xy[,1] %in% c("SonicB"),],sensor.text.args=list(labels="",cex=2),points.args = list(pch = 20, cex = 2, col ="orange"),add=TRUE)
  plot(Sonics_QV3_xy[Sonics_QV3_xy[,1] %in% c("SonicB"),],sensor.text.args=list(labels="",cex=2),points.args = list(pch = 20, cex = 2, col ="lightblue"),add=TRUE)
  plot(WeatherStation_xy[WeatherStation_xy[,1] %in% c("WS2"),],sensor.text.args=list(labels="",cex=2),points.args = list(pch = 20, cex = 2, col ="red"),add=TRUE)
  lines_sec2xy(Sonics_MK_xy,"SonicB",1,iWDSbMK[1],lty=2, col="orange")
  lines_sec2xy(Sonics_MK_xy,"SonicB",1,iWDSbMK[2],lty=2, col="orange")
  lines_sec2xy(WeatherStation_xy,"WS2",1,iWDSbMK[1],lty=2, col="orange")
  lines_sec2xy(WeatherStation_xy,"WS2",1,iWDSbMK[2],lty=2, col="orange")
  # QV3
  lines_sec2xy(Sonics_QV3_xy,"SonicB",1,iWDSbQV3[1],lty=2)
  lines_sec2xy(Sonics_QV3_xy,"SonicB",1,iWDSbQV3[2],lty=2)
  lines_sec2xy(WeatherStation_xy,"WS2",1,iWDSbQV3[1],lty=2)
  lines_sec2xy(WeatherStation_xy,"WS2",1,iWDSbQV3[2],lty=2)



  ###############
  ### Sonic C ###
  ###############

  sC_dir <- dir(file.path(PathData,"Sonic","SonicC"), full.names = TRUE)
  sC_list <- lapply(sC_dir, readWindMaster_ascii)
  sC_data <- rbindlist(sC_list)
  d_t <- sC_data[, as.numeric(diff(Time)) / 2]
  md_t <- median(d_t)
  d_t[d_t > md_t*2] <- md_t 
  sC_ibts <- shift_dt(time_offset["SonicC"],as.ibts(sC_data[,-(1:2)], st = sC_data[, Time] - c(md_t, d_t), 
    et = sC_data[, Time] + c(d_t, md_t)))
  sC_ws <- sC_ibts[WS[, c("WS_0460", "WD_corr")]]
  sC_ws[,c("U_sC","WD_sC")] <- with(sC_ws,{
    out <- mapply(rotate_twoaxis,u=u, v=v, w=w, T=T, MoreArgs = list(c.system = "windmaster"))
    # browser()
    list(
      U = unlist(out["uprot",])
      ,WD = unlist(out["wd",])
      )
  })

  # Abweichung zu Nord
  sC_off <- 0
  sC_ws$WD_sCcorr <- (sC_ws$WD_sC + sC_off) %% 360
  colClasses(sC_ws)[c("WD_sC","WD_sCcorr")] <- "circ"
  sC_ws30 <- pool(sC_ws, "1mins") ## 10min statt 30min da wenig Daten
  sC_ws30$WS_m_sC <- sC_ws30$WD_corr - sC_ws30$WD_sCcorr
  sC_ws30$sC_m_WS <- sC_ws30$WD_sCcorr - sC_ws30$WD_corr
  sC_ws30$WDdiff <- pmin(
    sC_ws30$WS_m_sC %% 360
    ,sC_ws30$sC_m_WS %% 360
    )
  sC_ws30$WDdiff[which(sC_ws30$WDdiff == sC_ws30$sC_m_WS)] <- -sC_ws30$WDdiff[which(sC_ws30$WDdiff == sC_ws30$sC_m_WS)]
  # m_sC <- MASS::rlm(WDdiff ~ 1, sC_ws30, subset = U_sC > 2)
  m_sC_1 <- MASS::rlm(WDdiff ~ 1, sC_ws30["18.03.2021 15:10 - 21.03.2021 12:00"], subset = WS_0460 > 0.8 & (WD_corr < 39 | WD_corr < 62) & WDdiff < 6+30 & WDdiff > 6-30)
  m_sC_2 <- MASS::rlm(WDdiff ~ 1, sC_ws30["21.03.2021 15:10 - 26.03.2021 12:00"], subset = WS_0460 > 0.5 & (WD_corr < 40 | WD_corr > 78) & abs(WDdiff) < 40+40)
  wdoff_sC_1 <- (sC_off + coef(m_sC_1)) %% 360
  wdoff_sC_2 <- (sC_off + coef(m_sC_2)) %% 360
  # wdoff_sC_1 <- 5.326178
  # wdoff_sC_2 <- 39.65245


  C1 <- as.data.table(sC_ws30["18.03.2021 15:10 - 21.03.2021 12:00"])
  C2 <- as.data.table(sC_ws30["21.03.2021 15:10 - 26.03.2021 12:00"])

  C1[,summary(abs(WDdiff))]
  plot(WD_sCcorr ~ WD_corr,C1)
  points(WD_sCcorr ~ WD_corr,C1[U_sC < 0.8],col="red",pch=19)
  points(WD_sCcorr ~ WD_corr,C1[abs(WDdiff) > 6+30],col="red",pch=19)
  points(WD_sCcorr ~ WD_corr,C1[WD_corr > 39 & WD_corr < 62],col="red",pch=19)
  abline(v=c(39,62))
  plot(WDdiff~I((WD_corr+180)%% 360 - 180),C1)
  plot(WDdiff~I((WD_corr+180)%% 360 - 180),C1,xlim=c(-40,120))
  abline(v=c(38,75))


  C2[,summary(abs(WDdiff))]
  plot(WD_sCcorr ~ WD_corr,C2)
  points(WD_sCcorr ~ WD_corr,C2[U_sC < 0.5],col="red",pch=19)
  points(WD_sCcorr ~ WD_corr,C2[WD_corr > 40 & WD_corr < 78],col="red",pch=19)
  points(WD_sCcorr ~ WD_corr,C2[abs(WDdiff) > 40+40],col="red",pch=19)
  abline(v=c(40,78))
  abline(-40,1)


  iWDScMK <- c(39,58)
  iWDScQV3 <- c(40,78)
  PlotOnStaticMap(STO_Map)
  plot(Sonics_MK_xy[Sonics_MK_xy[,1] %in% c("SonicC"),],sensor.text.args=list(labels="",cex=2),points.args = list(pch = 20, cex = 2, col ="orange"),add=TRUE)
  plot(Sonics_QV3_xy[Sonics_QV3_xy[,1] %in% c("SonicC"),],sensor.text.args=list(labels="",cex=2),points.args = list(pch = 20, cex = 2, col ="lightblue"),add=TRUE)
  plot(WeatherStation_xy[WeatherStation_xy[,1] %in% c("WS2"),],sensor.text.args=list(labels="",cex=2),points.args = list(pch = 20, cex = 2, col ="red"),add=TRUE)
  lines_sec2xy(Sonics_MK_xy,"SonicC",1,iWDScMK[1],lty=2, col="orange")
  lines_sec2xy(Sonics_MK_xy,"SonicC",1,iWDScMK[2],lty=2, col="orange")
  lines_sec2xy(WeatherStation_xy,"WS2",1,iWDScMK[1],lty=2, col="orange")
  lines_sec2xy(WeatherStation_xy,"WS2",1,iWDScMK[2],lty=2, col="orange")
  # QV3
  lines_sec2xy(Sonics_QV3_xy,"SonicC",1,iWDScQV3[1],lty=2)
  lines_sec2xy(Sonics_QV3_xy,"SonicC",1,iWDScQV3[2],lty=2)
  lines_sec2xy(WeatherStation_xy,"WS2",1,iWDScQV3[1],lty=2)
  lines_sec2xy(WeatherStation_xy,"WS2",1,iWDScQV3[2],lty=2)






  ###############
  ### Sonic 2 ###
  ###############

  s2_dir <- dir(file.path(PathData,"Sonic","Sonic2"), full.names = TRUE)
  s2_list <- lapply(s2_dir, readWindMaster_old_ascii)
  s2_data <- rbindlist(s2_list)
  d_t <- s2_data[, as.numeric(diff(Time)) / 2]
  md_t <- median(d_t)
  d_t[d_t > md_t*2] <- md_t
  i_dt <- which(d_t == 0) ## hat d_t Werte, welche 0 sind. Diese muss ich ausschliessen, damit es funktioniert.
  s2_ibts <- shift_dt(time_offset["Sonic2"],as.ibts(s2_data[-i_dt,-(1:2)], st = s2_data[-i_dt, Time] - c(md_t, d_t[-i_dt]), 
    et = s2_data[-i_dt, Time] + c(d_t[-i_dt], md_t)))
  s2_ws <- s2_ibts[WS[, c("WS_0460", "WD_corr")]]
  s2_ws[,c("U_s2","WD_s2")] <- with(s2_ws,{
    out <- mapply(rotate_twoaxis,u=u, v=v, w=w, T=T, MoreArgs = list(c.system = "windmaster"))
    # browser()
    list(
      U = unlist(out["uprot",])
      ,WD = unlist(out["wd",])
      )
  })

  # Abweichung zu Nord
  s2_off <- 0
  s2_ws$WD_s2corr <- (s2_ws$WD_s2 + s2_off) %% 360
  colClasses(s2_ws)[c("WD_s2","WD_s2corr")] <- "circ"
  s2_ws30 <- pool(s2_ws, "1mins") ## 10min statt 30min da wenig Daten
  s2_ws30$WS_m_s2 <- s2_ws30$WD_corr - s2_ws30$WD_s2corr
  s2_ws30$s2_m_WS <- s2_ws30$WD_s2corr - s2_ws30$WD_corr
  s2_ws30$WDdiff <- pmin(
    s2_ws30$WS_m_s2 %% 360
    ,s2_ws30$s2_m_WS %% 360
    )
  s2_ws30$WDdiff[which(s2_ws30$WDdiff == s2_ws30$s2_m_WS)] <- -s2_ws30$WDdiff[which(s2_ws30$WDdiff == s2_ws30$s2_m_WS)]
  # m_s2 <- MASS::rlm(WDdiff ~ 1, s2_ws30, subset = U_s2 > 2)
  m_s2_1 <- MASS::rlm(WDdiff ~ 1, s2_ws30["07.12.2020 10:10 - 15.12.2020 01:00"], subset = U_s2 > 0.5 & abs(WDdiff) > 102-30 & abs(WDdiff) < 102+30)
  m_s2_2 <- MASS::rlm(WDdiff ~ 1, s2_ws30["18.01.2021 15:00 - 27.01.2021 12:00"], subset = U_s2 > 0.8 & abs(WDdiff) > 84-30 & abs(WDdiff) < 84+30)
  m_s2_3 <- MASS::rlm(WDdiff ~ 1, s2_ws30["18.03.2021 16:00 - 21.03.2021 13:00"], subset = U_s2 > 0.8 & (WD_corr < 38 | WD_corr > 75) & WDdiff > -112-30 & WDdiff < -112+30)
  m_s2_4 <- MASS::rlm(WDdiff ~ 1, s2_ws30["21.03.2021 15:00 - 24.03.2021 00:00"], subset = U_s2 > 0.5 & (WD_corr < 38 | WD_corr > 75) & WDdiff > -6-30 & WDdiff < -6+30)
  wdoff_s2_1 <- (s2_off + coef(m_s2_1)) %% 360
  wdoff_s2_2 <- (s2_off + coef(m_s2_2)) %% 360
  wdoff_s2_3 <- (s2_off + coef(m_s2_3)) %% 360
  wdoff_s2_4 <- (s2_off + coef(m_s2_4)) %% 360
  # wdoff_s2_1 <- 259.7256 
  # wdoff_s2_2 <- 83.94708
  # wdoff_s2_3 <- 247.1232
  # wdoff_s2_4 <- 353.3588


  D1 <- as.data.table(s2_ws30["07.12.2020 10:10 - 15.12.2020 01:00"])
  D2 <- as.data.table(s2_ws30["18.01.2021 15:00 - 27.01.2021 12:00"])
  D3 <- as.data.table(s2_ws30["18.03.2021 16:00 - 21.03.2021 13:00"])
  D4 <- as.data.table(s2_ws30["21.03.2021 15:00 - 24.03.2021 00:00"])

  D1[,summary(abs(WDdiff))]
  plot(WD_s2corr ~ WD_corr,D1)
  points(WD_s2corr ~ WD_corr,D1[U_s2 < 0.5],col="red",pch=19)
  points(WD_s2corr ~ WD_corr,D1[abs(WDdiff) > 102+30 | abs(WDdiff) < 102-30],col="red",pch=19)
  abline(102,1)

  D2[,summary(abs(WDdiff))]
  plot(WD_s2corr ~ WD_corr,D2)
  points(WD_s2corr ~ WD_corr,D2[U_s2 < 0.8],col="red",pch=19)
  points(WD_s2corr ~ WD_corr,D2[abs(WDdiff) > 84+30 | abs(WDdiff) < 84-30],col="red",pch=19)
  abline(-84,1)

  D3[,summary(abs(WDdiff))]
  plot(WD_s2corr ~ WD_corr,D3)
  points(WD_s2corr ~ WD_corr,D3[U_s2 < 0.8],col="red",pch=19)
  points(WD_s2corr ~ WD_corr,D3[WD_corr > 38 & WD_corr < 75],col="red",pch=19)
  abline(v=c(38,75))
  abline(112,1)

  D4[,summary(abs(WDdiff))]
  plot(WD_s2corr ~ WD_corr,D4)
  points(WD_s2corr ~ WD_corr,D4[U_s2 < 0.5],col="red",pch=19)
  points(WD_s2corr ~ WD_corr,D4[WD_corr > 38 & WD_corr < 75],col="red",pch=19)
  points(WD_s2corr ~ WD_corr,D4[abs(WDdiff) > 6+30],col="red",pch=19)
  abline(v=c(38,62))
  abline(6,1)

  iWDS2MK <- c(38,75)
  iWDS2QV3 <- c(38,75)
  PlotOnStaticMap(STO_Map)
  plot(Sonics_MK_xy[Sonics_MK_xy[,1] %in% c("Sonic2"),],sensor.text.args=list(labels="",cex=2),points.args = list(pch = 20, cex = 2, col ="orange"),add=TRUE)
  plot(Sonics_QV3_xy[Sonics_QV3_xy[,1] %in% c("Sonic2"),],sensor.text.args=list(labels="",cex=2),points.args = list(pch = 20, cex = 2, col ="lightblue"),add=TRUE)
  plot(WeatherStation_xy[WeatherStation_xy[,1] %in% c("WS2"),],sensor.text.args=list(labels="",cex=2),points.args = list(pch = 20, cex = 2, col ="red"),add=TRUE)
  lines_sec2xy(Sonics_MK_xy,"Sonic2",1,iWDS2MK[1],lty=2, col="orange")
  lines_sec2xy(Sonics_MK_xy,"Sonic2",1,iWDS2MK[2],lty=2, col="orange")
  lines_sec2xy(Sonics_MK_xy,"Sonic2",1,75,lty=2, col="pink")
  lines_sec2xy(WeatherStation_xy,"WS2",1,iWDS2MK[1],lty=2, col="orange")
  lines_sec2xy(WeatherStation_xy,"WS2",1,iWDS2MK[2],lty=2, col="orange")
  lines_sec2xy(WeatherStation_xy,"WS2",1,75,lty=2, col="pink")
  # QV3
  lines_sec2xy(Sonics_QV3_xy,"Sonic2",1,iWDS2QV3[1],lty=2)
  lines_sec2xy(Sonics_QV3_xy,"Sonic2",1,iWDS2QV3[2],lty=2)
  lines_sec2xy(WeatherStation_xy,"WS2",1,iWDS2QV3[1],lty=2)
  lines_sec2xy(WeatherStation_xy,"WS2",1,iWDS2QV3[2],lty=2)


##############################
##############################
#####                    #####
#####    Sonic 10 min    #####
#####                    #####
##############################
##############################



##############################
### Start Sonic evaluation ###
##############################


###############
### Sonic A ###
###############
# "21.01.2021 16:00 - 27.01.2021 13:00" 
# "05.03.2021 15:00 - 10.03.2021 16:00"
# "19.03.2021 09:30 - 21.03.2021 12:00"
# "21.03.2021 15:10 - 26.03.2021 11:00"

  wdoff_sA_1 <- 294.1486
  wdoff_sA_2 <- 13.5 
  wdoff_sA_3 <- 359.4507
  wdoff_sA_4 <- 20.9562


  SonicA.1_raw <- evalSonic(
    file_directory = file.path(PathData,"Sonic","SonicA")
    , start_time = "21.01.2021 16:00"
    , end_time = "27.01.2021 13:00"
    , add_time = ZvSonicA
    , tz_sonic = "Etc/GMT-1"
    , avg_period = "10mins"
    , z_canopy = canopy
    , z_sonic = 2.34
    , dev_north = wdoff_sA_1
    )
  ### as ibts
  SonicA.1 <- as.ibts(SonicA.1_raw,st="start_interval",et="end_interval",tz="Etc/GMT-1")
  colClasses(SonicA.1)["WD"] <- "circ"

 
  SonicA.2_raw <- evalSonic(
    file_directory = file.path(PathData,"Sonic","SonicA")
    , start_time = "05.03.2021 15:00"
    , end_time = "10.03.2021 16:00"
    , add_time = ZvSonicA
    , tz_sonic = "Etc/GMT-1"
    , avg_period = "10mins"
    , z_canopy = canopy
    , z_sonic = 2.34
    , dev_north = wdoff_sA_2
    )
  ### as ibts
  SonicA.2 <- as.ibts(SonicA.2_raw,st="start_interval",et="end_interval",tz="Etc/GMT-1")
  colClasses(SonicA.2)["WD"] <- "circ"

  SonicA.3_raw <- evalSonic(
    file_directory = file.path(PathData,"Sonic","SonicA")
    , start_time = "19.03.2021 09:30"
    , end_time = "21.03.2021 12:00"
    , add_time = ZvSonicA
    , tz_sonic = "Etc/GMT-1"
    , avg_period = "10mins"
    , z_canopy = canopy
    , z_sonic = 2.16
    , dev_north = wdoff_sA_3
    # , create_graphs = TRUE
    # , save_directory = paste0(dirname(file.path(PathData,"Sonic","SonicA")),"/TS_covar_plots_SonicA")
    # , add_name = "SonicA"
    )
  ### as ibts
  SonicA.3 <- as.ibts(SonicA.3_raw,st="start_interval",et="end_interval",tz="Etc/GMT-1")
  colClasses(SonicA.3)["WD"] <- "circ"

  SonicA.4_raw <- evalSonic(
    file_directory = file.path(PathData,"Sonic","SonicA")
    , start_time = "21.03.2021 15:10"
    , end_time = "26.03.2021 11:00"
    , add_time = ZvSonicA
    , tz_sonic = "Etc/GMT-1"
    , avg_period = "10mins"
    , z_canopy = canopy
    , z_sonic = 0.83
    , dev_north = wdoff_sA_4
    )
  ### as ibts
  SonicA.4 <- as.ibts(SonicA.4_raw,st="start_interval",et="end_interval",tz="Etc/GMT-1")
  colClasses(SonicA.4)["WD"] <- "circ"

  SonicA <- rbind(SonicA.1,SonicA.2,SonicA.3,SonicA.4)
# saveRDS(SonicA,file = file.path(PathRSaves,"STO_SonicA_10min.rds"))


# plot(WS700_2[[1]]$data["01.01.2021 to ","WD_corr"])
# lines(SonicA[,"WD"],col="red")

###############
### Sonic B ###
###############
# "21.01.2021 15:00 - 27.01.2021 14:00"
# "18.03.2021 15:30 - 22.03.2021 00:00"
# "22.03.2021 - 26.03.2021 16:00"

  wdoff_sB_1 <- 123.3213
  wdoff_sB_2 <- 0.0974
  wdoff_sB_3 <- 312.44

  SonicB.1_raw <- evalSonic(
    file_directory = file.path(PathData,"Sonic","SonicB")
    , start_time = "21.01.2021 00:00"
    , end_time = "27.01.2021 13:00"
    , add_time = ZvSonicB
    , tz_sonic = "Etc/GMT-1"
    , avg_period = "10mins"
    , z_canopy = canopy
    , z_sonic = 2.34
    , dev_north = wdoff_sB_1
    )
  ### as ibts
  SonicB.1 <- as.ibts(SonicB.1_raw,st="start_interval",et="end_interval",tz="Etc/GMT-1")
  colClasses(SonicB.1)["WD"] <- "circ"

  SonicB.2_raw <- evalSonic(
    file_directory = file.path(PathData,"Sonic","SonicB")
    , start_time = "18.03.2021 15:30"
    , end_time = "22.03.2021"
    , add_time = ZvSonicB
    , tz_sonic = "Etc/GMT-1"
    , avg_period = "10mins"
    , z_canopy = canopy
    , z_sonic = 2.16
    , dev_north = wdoff_sB_2
    # , create_graphs = TRUE
    # , save_directory = paste0(dirname(file.path(PathData,"Sonic","SonicB")),"/TS_covar_plots_SonicB")
    # , add_name = "SonicB"
    )
  ### as ibts
  SonicB.2 <- as.ibts(SonicB.2_raw,st="start_interval",et="end_interval",tz="Etc/GMT-1")
  colClasses(SonicB.2)["WD"] <- "circ"

  SonicB.3_raw <- evalSonic(
    file_directory = file.path(PathData,"Sonic","SonicB")
    , start_time = "22.03.2021 00:00"
    , end_time = "26.03.2021 16:00"
    , add_time = ZvSonicB
    , tz_sonic = "Etc/GMT-1"
    , avg_period = "10mins"
    , z_canopy = canopy
    , z_sonic = 1.35
    , dev_north = wdoff_sB_3
    )

  ### as ibts
  SonicB.3 <- as.ibts(SonicB.3_raw,st="start_interval",et="end_interval",tz="Etc/GMT-1")
  colClasses(SonicB.3)["WD"] <- "circ"

  SonicB <- rbind(SonicB.1,SonicB.2,SonicB.3)
# saveRDS(SonicB,file = file.path(PathRSaves,"STO_SonicB_10min.rds"))


# plot(WS700_2[[1]]$data["01.01.2021 to ","WD_corr"])
# lines(SonicB[,"WD"],col="red")

###############
### Sonic C ###
###############
# "18.03.2021 15:10 - 21.03.2021 12:00"
# "21.03.2021 15:10 - 26.03.2021 12:00"

  wdoff_sC_1 <- 5.326178
  wdoff_sC_2 <- 39.65245

  SonicC.1_raw <- evalSonic(
    file_directory = file.path(PathData,"Sonic","SonicC")
    , start_time = "18.03.2021 15:10"
    , end_time = "21.03.2021 12:00"
    , add_time = ZvSonicC
    , tz_sonic = "Etc/GMT-1"
    , avg_period = "10mins"
    , z_canopy = canopy_c
    , z_sonic = 2.16
    , dev_north = wdoff_sC_1
    # , create_graphs = TRUE
    # , save_directory = paste0(dirname(file.path(PathData,"Sonic","SonicC")),"/TS_covar_plots_SonicC")
    # , add_name = "SonicC"
    )
  ### as ibts
  SonicC.1 <- as.ibts(SonicC.1_raw,st="start_interval",et="end_interval",tz="Etc/GMT-1")
  colClasses(SonicC.1)["WD"] <- "circ"

  SonicC.2_raw <- evalSonic(
    file_directory = file.path(PathData,"Sonic","SonicC")
    , start_time = "21.03.2021 15:10"
    , end_time = "26.03.2021 12:00"
    , add_time = ZvSonicC
    , tz_sonic = "Etc/GMT-1"
    , avg_period = "10mins"
    , z_canopy = canopy
    , z_sonic = 3.00
    , dev_north = wdoff_sC_2
    )
  ### as ibts
  SonicC.2 <- as.ibts(SonicC.2_raw,st="start_interval",et="end_interval",tz="Etc/GMT-1")
  colClasses(SonicC.2)["WD"] <- "circ"

  SonicC <- rbind(SonicC.1,SonicC.2)
# saveRDS(SonicC,file = file.path(PathRSaves,"STO_SonicC_10min.rds"))


# plot(WS700_2[[1]]$data["01.01.2021 to ","WD_corr"])
# lines(SonicC[,"WD"],col="red")

###############
### Sonic 2 ###
###############
# "07.12.2020 10:10 - 15.12.2020 01:00"
# "18.01.2021 15:00 - 27.01.2021 12:00"
# "18.03.2021 16:00 - 21.03.2021 13:00"
# "21.03.2021 15:00 - 24.03.2021 00:00"

  wdoff_s2_1 <- 259.7256 
  wdoff_s2_2 <- 83.94708
  wdoff_s2_3 <- 247.1232  
  wdoff_s2_4 <- 353.3588

  Sonic2.1_raw <- evalSonic(
    file_directory = file.path(PathData,"Sonic","Sonic2")
    , start_time = "07.12.2020 10:10"
    , end_time = "15.12.2020 01:00"
    , add_time = ZvSonic2
    , tz_sonic = "Etc/GMT-1"
    , avg_period = "10mins"
    , z_canopy = canopy
    , z_sonic = 2.34 # nicht sicher ob korrekt
    , dev_north = wdoff_s2_1
    )
  ### as ibts
  Sonic2.1 <- as.ibts(Sonic2.1_raw,st="start_interval",et="end_interval",tz="Etc/GMT-1")
  colClasses(Sonic2.1)["WD"] <- "circ"

  Sonic2.2_raw <- evalSonic(
    file_directory = file.path(PathData,"Sonic","Sonic2")
    , start_time = "18.01.2021 15:00"
    , end_time = "27.01.2021 12:00"
    , add_time = ZvSonic2
    , tz_sonic = "Etc/GMT-1"
    , avg_period = "10mins"
    , z_canopy = canopy
    , z_sonic = 2.34
    , dev_north = wdoff_s2_2
    )
  ### as ibts
  Sonic2.2 <- as.ibts(Sonic2.2_raw,st="start_interval",et="end_interval",tz="Etc/GMT-1")
  colClasses(Sonic2.2)["WD"] <- "circ"

  Sonic2.3_raw <- evalSonic(
    file_directory = file.path(PathData,"Sonic","Sonic2")
    , start_time = "18.03.2021 16:00"
    , end_time = "21.03.2021 13:00"
    , add_time = ZvSonic2
    , tz_sonic = "Etc/GMT-1"
    , avg_period = "10mins"
    , z_canopy = canopy
    , z_sonic = 2.16
    , dev_north = wdoff_s2_3
    # , create_graphs = TRUE
    # , save_directory = paste0(dirname(file.path(PathData,"Sonic","Sonic2")),"/TS_covar_plots_Sonic2")
    # , add_name = "Sonic2"
    )
  ### as ibts
  Sonic2.3 <- as.ibts(Sonic2.3_raw,st="start_interval",et="end_interval",tz="Etc/GMT-1")
  colClasses(Sonic2.3)["WD"] <- "circ"

  Sonic2.4_raw <- evalSonic(
    file_directory = file.path(PathData,"Sonic","Sonic2")
    , start_time = "21.03.2021 15:00"
    , end_time = "24.03.2021 00:00"
    , add_time = ZvSonic2
    , tz_sonic = "Etc/GMT-1"
    , avg_period = "10mins"
    , z_canopy = canopy
    , z_sonic = 2.00
    , dev_north = wdoff_s2_4
    )
  ### as ibts
  Sonic2.4 <- as.ibts(Sonic2.4_raw,st="start_interval",et="end_interval",tz="Etc/GMT-1")
  colClasses(Sonic2.4)["WD"] <- "circ"

  Sonic2 <- rbind(Sonic2.1,Sonic2.2,Sonic2.3,Sonic2.4)
# saveRDS(Sonic2,file = file.path(PathRSaves,"STO_Sonic2_10min.rds"))

# plot(WS700_2[[1]]$data[,"WD_corr"])
# lines(Sonic2[,"WD"],col="red")

#################
### Save Data ###
#################

# save all in one 
save(SonicA,SonicB,SonicC,Sonic2, file = file.path(PathRSaves,"Sonics_10min.RData"))
