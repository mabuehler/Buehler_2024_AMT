
contourXZ <- function (
  Source,stats="CE",p1=NULL,p2=NULL,Ustar=NULL,L=NULL,Zo=NULL,Suu=NULL,Svu=NULL,
  Swu=NULL,z_sWu=NULL,WD=NULL,d=NULL, pDist = 100,
  z.lower=max(Interval[,"Zo [m]"]*2+Interval[,"d [m]"]),z.upper=10,nz=10,
  logz=FALSE,nx=20,levels=function(x)quantile(range(x,na.rm=TRUE),c(0.01,0.05,0.1,0.5,0.9)),
  zlim=c(z.lower,z.upper),Rug=TRUE,locatorWindow=c(x0=-200,x1=200,y0=-200,y1=200),
  ncores=NULL,Cat.Path=NULL,Overview=TRUE,add=FALSE,Tolerances=NULL,Sensors=NULL,
  Interval=NULL,Model=NULL,MaxFetch=NULL,N0=NULL,addWR=TRUE,WRpos=ceiling(WDmean/90),
  addSource=TRUE,SourceCol="darkgreen",addSB=TRUE,SBpos=1,labels=sprintf("%1.2f",levels(xz)),
  fill=FALSE,v_d=NULL,quiet=FALSE,...) {

  if (!inherits(Source, "contourXZ")) {

    if (inherits(Source, "InputList")) {
      if(is.null(Source$Sources))stop("At least one source must be provided!")
      # Interval:
      if(is.null(Interval)){
        Interval <- Source$Interval
      }
      # Model
      if(is.null(Model))Model <- Source$Model
      # Tolerances
      if(is.null(Tolerances))Tolerances <- Source$Tolerances
      # Sources
      Source <- Source$Sources
    } else if(inherits(Source,"Sources")){   
      # Interval:
      if(is.null(Interval)){
        if(length(wn <- which(sapply(list(Ustar,L,Zo,Suu,Svu,Swu,z_sWu,WD,d,N0,MaxFetch),is.null))))stop(paste("Insufficient windstatistics and windprofiles!\nPlease supply:\n",paste(c("- Ustar","- L","- Zo","- Suu","- Svu","- Swu","- z_sWu","- WD","- d")[wn],collapse="\n"),"\n",sep=""))
        Interval <- genInterval(Data=data.frame(Ustar,L,Zo,Suu,Svu,Swu,z_sWu,WD,d,N0,MaxFetch,"","",stringsAsFactors=FALSE))
      }
      # Model
      if(is.null(Model))Model <- genModel()
      # Tolerances
      if(is.null(Tolerances))Tolerances <- genTolerances()
    } else {
      stop("Argument Source must be one of the following classes: contourXZ, InputList or Sources!")
    }

    if(!is.function(levels)){
      sub_levels <- eval(substitute(levels))
      levels <- function(x)sub_levels
    }
    
    if(length(wn <- which(!sapply(list(Ustar,L,Zo,Suu,Svu,Swu,z_sWu,WD,d,N0,MaxFetch),is.null)))){
      IntervalIn <- do.call(genInterval, list(Ustar = Ustar, L = L, Zo = Zo,
        sUu = Suu, sVu = Svu, sWu = Swu, z_sWu = z_sWu, WD = WD, d = d, N0 = N0,
        MaxFetch = MaxFetch)[wn])
      Interval[, wn] <- IntervalIn[, wn]
    }
    WDmean <- (360 + atan2(sum(sin(Interval[,"WD [deg N]"]/180*pi)*Interval[,"Ustar [m/s]"])/sum(Interval[,"Ustar [m/s]"]),sum(cos(Interval["WD [deg N]"]/180*pi)*Interval[,"Ustar [m/s]"])/sum(Interval[,"Ustar [m/s]"]))/pi*180) %% 360
    Heights <- seq(z.lower,z.upper,length.out=nz)
    # Interval[1,10] <- paste0("Sensor",seq(nz),collapse=",")
    if(!(p2n <- is.null(p2)) && length(p2) == 1){
      if(is.numeric(p2)){
        wdrad <- (90 - p2) / 180 * pi
        p2 <- p1 + c(cos(wdrad), sin(wdrad)) * pDist[1]
      } else {
        wdrad <- (90 - WDmean) / 180 * pi
        p2 <- p1 + c(cos(wdrad), sin(wdrad)) * pDist[1]
      }
      if(length(pDist) == 2){
        p1 <- p1 + c(cos(wdrad), sin(wdrad)) * pDist[2]        
      }
    }
    if((p1n<-is.null(p1))|p2n){
      sm <- colMeans(Source[,2:3])
      siteMap(Source,main=paste0("Please locate ",paste(c("starting point p1","ending point p2")[which(c(p1n,p2n))],collapse=" and ")),xlim=locatorWindow[1:2]+sm[1],ylim=locatorWindow[3:4]+sm[2])
      addWindrose(WDmean)
      addScaleBar()
      p12 <- locator(p1n+p2n)
      graphics.off()
      if(p1n&p2n){
        p1 <- c(p12$x[1],p12$y[1])
        p2 <- c(p12$x[2],p12$y[2])
      } else if(p1n){
        p1 <- c(p12$x[1],p12$y[1])
      } else {
        p2 <- c(p12$x[1],p12$y[1])
      }
    }
    
    if(!is.null(ncores)) Model["ncores"] <- ncores

    px <- seq(p1[1],p2[1],length.out=nx)
    py <- seq(p1[2],p2[2],length.out=nx)
    xm <- sqrt((px[1]-px)^2+(py[1]-py)^2)
    zm <- Heights
    xz <- matrix(0,nrow=length(xm),ncol=length(zm))
    dimnames(xz) <- list(xm,zm)
    dm <- length(xm)
    Sindex <- expand.grid(seq.int(nx),seq.int(nz))
    SensorNames <- paste0("Sensor_",apply(Sindex,1,paste,collapse="_"))

    # browser()
    Sensors <- genSensors(data.frame(name = SensorNames, x= px[Sindex[[1]]], y = py[Sindex[[1]]], z = zm[Sindex[[2]]]))

    Interval[, c('Sensor Names (sep = ",")', 'Source Names (sep = ",")')] <- ""
    Sources <- Source
    Sources[,4] <- as.numeric(as.factor(paste(Sources[,1],Sources[,4])))
    Sources[,1] <- "Source"
    InL <- genInputList(Interval=Interval, Sensors = Sensors, Model = Model, Tolerances = Tolerances, Sources = Source)

    if(quiet){
      blackhole <- capture.output(x <- runbLS(InL,Cat.Path=Cat.Path,asDT=TRUE))
    } else {
      x <- runbLS(InL,Cat.Path=Cat.Path,asDT=TRUE)
    }
    if(depo <- !is.null(v_d)){
      y <- copy(x)
      x <- deposition(y,v_d)
      if(quiet){
        blackhole <- capture.output(x <- deposition(y,v_d))
      } else {
        x <- deposition(y,v_d)
      }
    }
    if(nrow(Interval)>1){
      warning("Averaging several intervals!")
      # xz[] <- eval(parse(text=paste0("x[,mean(",stats,"),keyby=Sensor][SensorNames,V1]")))
      xz[] <- x[, .(V1 = eval(parse(text=paste0("mean(",stats,")")))),
      keyby=Sensor][SensorNames, V1]
    } else {
      setkey(x,Sensor)  
      xz[] <- x[SensorNames, .(V1 = eval(parse(text=stats)))][,V1]
    }
  } else {
    xm <- Source$Img$x
    zm <- Source$Img$y
    xz <- Source$Img$z
    SensorNames <- Source$SensorNames
    px <- Source$Pts$x
    py <- Source$Pts$y
    p1 <- c(px[1],py[1])
    p2 <- c(rev(px)[1],rev(py)[1])
    depo <- Source$depo
    vDep <- Source$v_d
    x <- Source$bLSout
    Sensors <- attr(x,"ModelInput")$Sensors
    Src <- attr(x,"ModelInput")$Sources
    Interval <- attr(x,"ModelInput")$Interval
    WDmean <- (360 + atan2(sum(sin(Interval[,8]/180*pi)*Interval[,2])/sum(Interval[,2]),sum(cos(Interval[,8]/180*pi)*Interval[,2])/sum(Interval[,2]))/pi*180) %% 360
    if(!missing(v_d)&&!all.equal(vDep,v_d)){
      depo <- TRUE
      y <- copy(x)
      x <- deposition(y,v_d)
      if(quiet){
        blackhole <- capture.output(x <- deposition(y,v_d))
      } else {
        x <- deposition(y,v_d)
      }
    } else {
      v_d <- vDep
      if(depo){
        x <- Source$bLSdep
        y <- Source$bLSout
      }
    }
    if(missing(stats)){
      stats <- Source$stats
    } else {
      if(nrow(Interval)>1){
        warning("Averaging several intervals!")
        # xz[] <- eval(parse(text=paste0("x[,mean(",stats,"),keyby=Sensor][SensorNames,V1]")))
        xz[] <- x[, .(V1 = eval(parse(text=paste0("mean(",stats,")")))),
        keyby=Sensor][SensorNames, V1]
      } else {
        xz[] <- x[SensorNames, .(V1 = eval(parse(text=stats)))][,V1]
      }
    }
    Source <- Src
  }

  srX <- range(Source[,2])
  srY <- range(Source[,3])
  
  contour(xm,zm,xz,levels=levels(xz),add=add,labels=labels,...)
  
  if(fill){
    .filled.contour(xm,zm,xz, levels=c(levels(xz),max(xz)), col=rev(heat.colors(length(levels(xz)))))
    # .filled.contour(xm,zm,xz, levels=c(levels(xz),max(xz)), col=heat.colors(length(levels(xz))))
    contour(xm,zm,xz,levels=levels(xz),add=TRUE,labels=labels,...)
  }

  zGround <- max(Interval[,"Zo [m]"]+Interval[,"d [m]"])
  abline(h=zGround)
  abline(h=min(Interval[,"Zo [m]"]+Interval[,"d [m]"]),lty=2)
  if(addSource){
    srca <- attr(bLSmodelR:::procSources(Source),"SourceList")
    for(i in seq_along(srca)){
      srcb <- srca[[i]]
      upoly <- unique(srcb[,3])
      for(p in upoly){
        srcc <- srcb[srcb[,3] == p,]
        xSource <- numeric(0)
        ySource <- numeric(0)
        for(j in 1:(NROW(srcc)-1)){
          Pxy <- line.intersection(p1, p2, srcc[j,1:2], srcc[j + 1,1:2], 
            check.inside = "both")
          if(!is.na(Pxy[1])){
            xSource <- c(xSource,Pxy[[1]])
            ySource <- c(ySource,Pxy[[2]])
          }
        }
        distSource <- sqrt((px[1]-xSource)^2+(py[1]-ySource)^2)
        dSource <- sort(distSource)
        lSource <- length(dSource) 
        if(lSource>=1){
          for(j in seq(1,lSource,2)){
            if((j+1)<=lSource){
              lines(dSource[j+c(0,1)],c(zGround,zGround),col=SourceCol,lwd=3)
            }       
          }
          points(dSource,rep(zGround,length(dSource)),pch=20,col=SourceCol,lwd=3)
        }
      }
    }
  }
  if(Overview){
    if(Rug){
      rug(xm,0.02,lwd=2)    
      rug(zm,0.02,lwd=2,side=2)
    } 
    opar <- par()
    figNew <- opar$fig[c(2,4)]/max(opar$fin)*opar$fin[2:1]
    par(new=TRUE,fig=c(0.2*figNew[1],0.4*figNew[1],0.6*figNew[2],0.8*figNew[2]),mar=c(0,0,0,0),cex=0.4)
    rx <- range(range(px)+diff(range(px))*0.1*c(-1,1),srX)
    ry <- range(range(py)+diff(range(py))*0.1*c(-1,1),srY)
    # siteMap(Source,xlim=rx,ylim=ry)
    plot(rx,ry,asp=1,type="n",ylab="",xlab="")
    usr <- par()$usr
    rect(usr[1],usr[3],usr[2],usr[4],col="white",border=1)
    grid()
    box()
    siteMap(Source,add=TRUE)
    # pg <- srca[,2:4]
    # upg <- unique(pg[,3])
    # for(k in upg){
    #   pindex <- which(pg[,3]==k)
    #   pgxy <- pg[pindex,1:2]
    #   if(sum(abs(pgxy[1,]-pgxy[NROW(pgxy),]))>1E-3)pgxy <- rbind(pgxy,pgxy[1,])
    #   lines(pgxy,col=SourceCol)
    # }
    lines(c(p1[1],p2[1]),c(p1[2],p2[2]))
    if(Rug)points(px,py,pch=20)
    text(p1[1],p1[2],"0",pos=3,offset=0.2)
    text(p2[1],p2[2],sprintf("%1.1f",rev(xm)[1]),pos=3,offset=0.2)
    if(!is.null(Sensors)){
      points(Sensors[,2:3],pch=3,lwd=2)
      LS <- unique(Sensors[,5]) %w/o% ""
      for(i in LS){
        lines(Sensors[Sensors[,5]%in%i,2:3],lty=3)
      }
    }
    if(addWR)addWindrose(WDmean,WRpos)
    if(addSB)addScaleBar(SBpos)
    op <- options()
    options(warn=-1)
    par(opar)
    options(op)
  }
  if(addSB)addScaleBar(1)
  if(depo){
    tmp <- copy(y)
    y <- copy(x)
    x <- tmp
  }
  return(invisible(structure(
    list(
      Img=list(x=xm,y=zm,z=xz),
      stats=stats,
      SensorNames=SensorNames,
      Pts=list(x=px,y=py),
      bLSout=x,
      depo=depo,
      bLSdep=if(depo)y else NULL,
      v_d=v_d
      )
    ,class=c("contourXZ","list")
    )))
}



line.intersection <- function(l1_p1, l1_p2, l2_p1, l2_p2, 
  check.inside = c("none","both","l1","l2")){
  x1 <- l1_p1[[1]]
  x2 <- l1_p2[[1]]
  x3 <- l2_p1[[1]]
  x4 <- l2_p2[[1]]
  y1 <- l1_p1[[2]]
  y2 <- l1_p2[[2]]
  y3 <- l2_p1[[2]]
  y4 <- l2_p2[[2]]
  d <- (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4) 
  if(d == 0){
    warning("lines are parallel")
    return(NA)
  }

  t <- ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / d
  u <- -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / d
  tinside <- t >= 0 & t <= 1
  uinside <- u >= 0 & u <= 1
  switch(check.inside[1]
    ,"both" = {
      if(!tinside || !uinside){
        return(NA)
      }
    }
    ,"l1" = {
      if(!tinside){
        return(NA)
      }
    }
    ,"l2" = {
      if(!uinside){
        return(NA)
      }
    }
  )
  if(tinside){
    c(x1 + t * (x2 - x1), y1 + t * (y2 - y1))
  } else {
    c(x3 + u * (x4 - x3), y3 + u * (y4 - y3))
  }

}
# line.intersection(c(0,0), c(1,0), c(1, 0.5), c(1, 1))
# line.intersection(c(0,0), c(1,0), c(1, 0.5), c(1, 1), check.inside = "l1")
# line.intersection(c(0,0), c(1,0), c(1, 0.5), c(1, 1), check.inside = "l2")
# line.intersection(c(0,0), c(1,0), c(1, 0.5), c(1, 1), check.inside = "both")


# line.intersection(c(0,0), c(1,0), c(1, 0), c(1, 1))
# line.intersection(c(0.75,0), c(1,0), c(0.5, -1), c(0.5, 1))
# line.intersection(c(0.75,0), c(1,0), c(0.5, -1), c(0.5, 1), check.inside = "l1")

# line.intersection(c(0,0), c(1,1), c(0, 1), c(1, 0))



if(FALSE){

  library(bLSmodelR)
  library(ibts)
  library(RgoogleMaps)

  Pfad <- "~/repos/4_Projects/3_FerEVS"

  PfadRSaves <- paste0(Pfad,"/Auswertung/RSaves")
  CatPath <- "~/Y-Drive/Marcel Bühler/Test"

  load(file=paste0(PfadRSaves,"/FerEVS_MK2_bLS_Geometry.RData"))
  load(paste0(PfadRSaves,"/All_MK2_GF_v1.RData"))
  Sonic1 <- readRDS(file=paste0(PfadRSaves,"/FerEVS_MK2_Sonic1_30min_pf_vanDik2004.rds"))
  Sonic2 <- readRDS(file=paste0(PfadRSaves,"/FerEVS_MK2_Sonic2_30min_pf_vanDik2004.rds"))
  
  All_sonic <- merge(All,Sonic1[,c("z_sonic","sWu","d")])
  All_sonic <- merge(All_sonic,Sonic2[,c("z_sonic","sWu","d")])

  indS1 <- which(All_sonic[,"WD"] > 315 | All_sonic[,"WD"] < 135)
  indS2 <- which(All_sonic[,"WD"] < 315 & All_sonic[,"WD"] > 135)

  All_sonic$z_sWu <- All_sonic$sWu <- All_sonic$d <- NA
  All_sonic[indS1]$z_sWu <- All_sonic[indS1]$z_sonic.x
  All_sonic[indS1]$sWu <- All_sonic[indS1]$sWu.x
  All_sonic[indS1]$d <- All_sonic[indS1]$d.x
  All_sonic[indS2]$z_sWu <- All_sonic[indS2]$z_sonic.y
  All_sonic[indS2]$sWu <- All_sonic[indS2]$sWu.y
  All_sonic[indS2]$d <- All_sonic[indS2]$d.y

  Sonic17 <- All_sonic[indAll17,c("Ustar","L","Zo","sUu","sVu","z_sWu" ,"WD","d")]
  SonicTest <- genInterval(as.data.frame(Sonic17))

  EVS_Map <- ReadMapTile(file.path(Pfad,"Figures/EVS_GoogleMaps_new.png"))

  source("~/repos/3_Scripts/3_Koordinatentransformation/WGS84_CH1903.R")

  ncores = 2

  Snc <- SonicTest[1:10,]
  # ncores = 4
  # Snc <- SonicTest[c(1,10),]
  Snc$N0 <- 2E3

  source("/home/christoph/repos/2_Packages/1_bLSmodelR/Functions Aufraeumen/ToDo/contourXZ.r")

  plot(Sources,Sensors)
  # locator(2)
  p1 <- c(711568.7, 260835.6)
  p2 <- c(711822.1, 261053.1)
  lines(c(p1[1],p2[1]),c(p1[2],p2[2]))

  XZ_out_17 <- contourXZ(Sources[Sources[,1] == "EVS",], p1 = p1, p2 = p2,
    nz = 20, nx = 100,
    Interval = Snc[1,], Cat.Path = CatPath
    )
    # levels = c(0.05, 0.1, 0.2, 0.5, 1, 2, 5),
    # fill = TRUE,ncores=ncores)

  contourXZ(XZ_out_17, fill = TRUE, levels = function(x)c(0.05, 0.1, 0.2, 0.5, 1, 2, 5))


  XZ_out_17 <- contourXZ(Sources[Sources[,1] == "EVS",], 
    p1 = colMeans(Sources[Sources[,1] == "EVS",2:3]), p2 = "WD",
    nz = 20, nx = 100, ncores = 2, pDist = c(50,-100),
    Interval = Snc[1,], Cat.Path = CatPath
    )

}
