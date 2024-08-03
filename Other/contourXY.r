
contourXY <- function (
  Source, z.meas, stats = "CE", xlim = c(-100, 100), ylim = c(-100, 100),
  origin = c(0, 0), MyMap = NULL, StaticMapArgs = NULL, alpha = 0.3,
  fill = FALSE, col = NULL, lty = NULL, lwd = NULL, coord2wgs84 = NULL,
  Ustar = NULL, L = NULL, Zo = NULL, Suu = NULL,
  Svu = NULL, Swu = NULL, z_sWu = NULL, WD = NULL, d = NULL, axs = c("r", "i"),
  showLegend = TRUE, dx = NULL, dy = NULL, levels = function (x) 
  quantile(range(x, na.rm = TRUE), c(0.01, 0.05, 0.1, 0.5, 0.9)),
  siteMapArgs = NULL, useSTRtree = TRUE, avoidGEOS = FALSE, sigNums = 1,
  use.sym = FALSE, use.avg = FALSE, ncores = NULL, Cat.Path = NULL, add = FALSE,
  Tolerances = NULL, Interval = NULL, Model = NULL, MaxFetch = NULL, N0 = NULL,
  lpos = NULL, leg.bg.col = NULL, addSource = FALSE, addWR = !add,
  WRpos = ceiling(WDmean / 90), WRfrac = 20, WRscale = 1,
  SourceCol = "darkgreen", addSB = !add, SBpos = 3, labels = sprintf("%1.2f",
    levels(xy)), v_d = NULL, quiet = FALSE, SourceName = NULL, ...) {


    # TODO:
    #       Fix without map
    #       Fix map XY not metric!!!

  if(!is.null(ncores) && ncores > 1) data.table::setDTthreads(ncores)
  if (!inherits(Source, "contourXY")) {

    if (inherits(Source, "InputList")) {
      if (is.null(Source$Sources)) stop("At least one source must be provided!")
      # Interval:
      if (is.null(Interval)) {
        Interval <- Source$Interval
      }
      # Model
      if (is.null(Model)) Model <- Source$Model
      # Tolerances
      if (is.null(Tolerances)) Tolerances <- Source$Tolerances
      # Sources
      Source <- Source$Sources
    } else if (inherits(Source, "Sources")) {
      # Interval:
      if(is.null(Interval)){
        if(length(wn <- which(sapply(list(Ustar,L,Zo,Suu,Svu,Swu,z_sWu,WD,d),is.null))))stop(paste("Insufficient windstatistics and windprofiles!\nPlease supply:\n",paste(c("- Ustar","- L","- Zo","- Suu","- Svu","- Swu","- z_sWu","- WD","- d")[wn],collapse="\n"),"\n",sep=""))
        Interval <- genInterval(Data=data.frame(Ustar,L,Zo,Suu,Svu,Swu,z_sWu,WD,d,"","",stringsAsFactors=FALSE))
      }
      # Model
      if(is.null(Model))Model <- genModel()
      # Tolerances
      if(is.null(Tolerances))Tolerances <- genTolerances()
    } else {
      stop("Argument Source must be one of the following classes: contourXY, InputList or Sources!")
    }
    
    if(length(wn <- which(!sapply(list(Ustar,L,Zo,Suu,Svu,Swu,z_sWu,WD,d,N0,MaxFetch),is.null)))){
      IntervalIn <- do.call(genInterval, list(Ustar = Ustar, L = L, Zo = Zo,
        sUu = Suu, sVu = Svu, sWu = Swu, z_sWu = z_sWu, WD = WD, d = d, N0 = N0,
        MaxFetch = MaxFetch)[wn])
      Interval[, wn] <- IntervalIn[, wn]
    }
    WDmean <- (360 + atan2(sum(sin(Interval[,"WD [deg N]"]/180*pi)*Interval[,"Ustar [m/s]"])/sum(Interval[,"Ustar [m/s]"]),sum(cos(Interval["WD [deg N]"]/180*pi)*Interval[,"Ustar [m/s]"])/sum(Interval[,"Ustar [m/s]"]))/pi*180) %% 360
    
    if(!is.null(ncores)) Model["ncores"] <- ncores

    if(!is.null(MyMap)){
      if(!("TrueProj" %in% names(StaticMapArgs))||!StaticMapArgs$TrueProj){
        if(is.null(coord2wgs84)){
          stop("Please provide a coordinate transformation to WGS84 (coord2wgs84)")
        } else {
          coord2wgs84Original <- coord2wgs84
          coord2wgs84 <- function(x,y,mymap=MyMap,...){
            xy <- coord2wgs84Original(x,y,...)
            LatLon2XY.centered(mymap,xy$lat,xy$lon)
          }
        }
      } else {
        stop("'TrueProj' == FALSE is not allowed")
      }
    } else {
      xlim <- xlim + origin[1]
      ylim <- ylim + origin[2]      
    }

    if(!add){
      if(is.null(MyMap)){
        do.call(siteMap,c(list(Source,add=add,xlim=xlim,ylim=ylim),siteMapArgs))
      } else {
        do.call(PlotOnStaticMap,c(list(MyMap=MyMap),StaticMapArgs))
      }
    }
    usr <- par()$usr
    xlim <- usr[1:2]
    ylim <- usr[3:4]
    if(is.null(dx)){
      dx <- floor(diff(xlim)/64)
    }
    if(is.null(dy)){
      dy <- floor(diff(ylim)/64)
    }

    xlimOriginal <- xlim
    ylimOriginal <- ylim

    if(!is.null(coord2wgs84)){
      lim1 <-  optim(colMeans(Source[,2:3]),function(x)sum((unlist(coord2wgs84(x=x[1],y=x[2])) - c(xlim[1],ylim[1]))^2))$par
      lim2 <-  optim(colMeans(Source[,2:3]),function(x)sum((unlist(coord2wgs84(x=x[1],y=x[2])) - c(xlim[2],ylim[2]))^2))$par
      xlim <- c(lim1[1],lim2[1])
      ylim <- c(lim1[2],lim2[2])
    }

    px <- seq(xlim[1]-dx/2,xlim[2]+dx/2,dx)
    py <- seq(ylim[1]-dy/2,ylim[2]+dy/2,dy)
    xm <- px[-1]-dx/2
    ym <- py[-1]-dy/2
    xy <- matrix(0,nrow=length(xm),ncol=length(ym))
    dimnames(xy) <- list(xm,ym)
    dm <- length(xm)
    Sindex <- expand.grid(seq.int(length(xm)),seq.int(length(ym)))
    SensorNames <- paste0("Sensor_",apply(Sindex,1,paste,collapse="_"))

    # create Sensors
    Sensors <- structure(
        data.frame(
            'Sensor Name' = SensorNames, 
            'Sensor ID' = '1',
            'Node' = 1L,
            'x-Coord (m)' = px[Sindex[[1]]],
            'y-Coord (m)' = py[Sindex[[2]]],
            'Sensor Height (m)' = z.meas,
            'Distance between Point-Sensors (m)' = 0,
            'Number of Point-Sensors' = 1L,
            stringsAsFactors = FALSE,
            check.names = FALSE
        ),
        class = c('Sensors', 'data.frame'),
        Version = '4.2+'
        )

    
    Interval[,c(12,13)] <- ""
    Sources <- Source
    if(!is.null(SourceName)){
      Sources <- Sources[Sources[,1] %chin% SourceName,]
      # attr(Sources,"SourceList") <- attr(Sources,"SourceList")[SourceName]
    }
    Sources[,4] <- as.numeric(as.factor(paste(Sources[,1],Sources[,4])))
    Sources[,1] <- "Source"
    # names(attr(Sources,"SourceList")) <- "Source"

    # check distance origin to source
    ydist_os <- max(abs(outer(ylim, Sources[, 3], '-')))
    xdist_os <- max(abs(outer(xlim, Sources[, 2], '-')))
    if (any(c(xdist_os, ydist_os) > 2e3)) {
        stop('Maximum distance xlim/ylim to Source > 2000 m. Did you set origin correctly?')
    }

    S1 <- Sensors[1,]
    InL <- genInputList(Interval=Interval, Sensors = Sensors, Model = Model, Tolerances = Tolerances, Sources = Source)
    ModelInput <- genInputList(Interval=Interval, Sensors = S1, Model = Model, Tolerances = Tolerances, Sources = Source)
    Slist <- attr(bLSmodelR:::procSources(InL[["Sources"]]),"SourceList")  
    SourceNames <- names(Slist)
    for(i in SourceNames){
      Slist[[i]] <- cbind(Feld=i,Slist[[i]],stringsAsFactors=FALSE)
    }
    Scalc <- rbindlist(Slist)
    setnames(Scalc,c("Plot","x","y","pid"))
    setkey(Scalc,Plot,pid)
    rm(Slist)
    Srange <- Scalc[,rbind(
      cbind(x=max(x)-min(xm),y=max(y)-min(ym))
      ,cbind(x=min(x)-max(xm),y=min(y)-max(ym))
      )]

    if (fill) {
        require(maptools)
        require(sp)
        require(rgeos)
    }

    if (parl <- ncores > 1) {
      on.exit({                                       
          sfStop()                                                      
      })
      sfInit(parallel = TRUE, cpus = ncores, type = "SOCK")         
      cl <- sfGetCluster()
      sfClusterSetupRNG(seed = sample.int(1e+09, 6, TRUE))                                                                 
      gwd <- getwd()                   
      # sfExport("gwd", "C.Path") 
      sfLibrary(bLSmodelR)
      invisible(sfClusterCall(setwd, gwd))
      invisible(clusterCall(cl, data.table::setDTthreads, 1L))
    }

    # wtf?!!! why is SensorNames reordered after this paragraph?!?
    IntervalRun <- prepareIntervals(ModelInput,Cat.Path,TRUE)
    bLSmodelR:::.calcCatalogs(IntervalRun,ModelInput,Cat.Path,parl)
    setDT(Sensors)
    setnames(Sensors,c("x-Coord (m)","y-Coord (m)"),c("x","y"))
    setkey(Sensors,"Sensor Name")
    Out <- IntervalRun[rep(1,length(SensorNames)),.(rn)][
      ,":="(Sensor=SensorNames,
        CE=0,CE_se=NA_real_,CE_lo=NA_real_,CE_hi=NA_real_,
        uCE=0,uCE_se=NA_real_,uCE_lo=NA_real_,uCE_hi=NA_real_,
        vCE=0,vCE_se=NA_real_,vCE_lo=NA_real_,vCE_hi=NA_real_,
        wCE=0,wCE_se=NA_real_,wCE_lo=NA_real_,wCE_hi=NA_real_,
        N_TD=0,TD_Time_avg=NA_real_,TD_Time_max=NA_real_,Max_Dist=NA_real_,UCE=0)][,rn:=NULL]
    setkey(Out,Sensor)

    if (parl) {
        # get indides
        nr <- nrow(Interval)
        ns <- length(SensorNames)
        # distribute cores to rows
        nc <- max(1, floor(ncores / nr))
        nresid <- max(0, ncores - nr * nc)
        # distribute residual cores to rows
        nall <- rep(nc, nr) + c(rep(1, nresid), rep(0, nr - nresid))
        # create indices
        ind_parl <- lapply(nall, function(x) splitIndices(ns, x))
        # prepare input
        InList <- unlist(lapply(seq.int(nr), function(i) {
            lapply(ind_parl[[i]], function(x, a, b, d) {
                sens <- b[x]
                list(a, sens, d[sens])
            }, a = IntervalRun[i, ], b = SensorNames, d = Out)
        }), recursive = FALSE)
        if(!quiet)cat("Calculating tiles in parallel...\n")
        out <- unlist(clusterApply(cl, InList, calcDep, Cat.Path, use.avg, 
                use.sym, quiet = TRUE, Srange, Scalc, Sensors), recursive = FALSE)
        if(!quiet)cat("done.\n")
        sfStop()
        on.exit()
        Out <- rbindlist(out)
    } else {
      Out <- calcDep(list(IntervalRun, SensorNames, Out), Cat.Path, use.avg, use.sym, quiet, 
        Srange, Scalc, Sensors)[[1]]
    }
    setkey(Out,Sensor)

    # check if all 0
    if (Out[, all(CE == 0)]) {
        plot(InL$Sensors, Sources, sensors.text.args = list(labels = ''))
        stop('No plume visible in selected xylim ranges!')
    }

    if(!quiet)cat("\n")
    if(nrow(Interval)>1){
      warning("Averaging several intervals!")
      Out <- Out[,lapply(.SD,base::mean,na.rm=TRUE),keyby=Sensor]
    }
    if(use.sym){
      Out[,N_TD:=N_TD/2]
    }

    # wtf?!!! why is SensorNames reordered here?!?
    xy[] <- Out[paste0("Sensor_",apply(Sindex,1,paste,collapse="_")), .(V1 = eval(parse(text=stats)))][,V1]

    # xy[] <- eval(parse(text=paste0("Out[SensorNames,",stats,"]")))
    # if(depo <- !is.null(v_d)){
    #   y <- copy(x)
    #   x <- deposition(y,v_d)
    #   if(quiet){
    #     blackhole <- capture.output(x <- deposition(y,v_d))
    #   } else {
    #     x <- deposition(y,v_d)
    #   }
    # }

  } else {
    if (fill) {
        require(maptools)
        require(sp)
        require(rgeos)
    }
    xm <- Source$Img$x
    dx <- diff(xm[1:2])
    ym <- Source$Img$y
    dy <- diff(ym[1:2])
    xy <- Source$Img$z
    xlim <- if(missing(xlim))Source$limits$xlim
    ylim <- if(missing(ylim))Source$limits$ylim
    InL <- Source$bLSin
    Sensors <- InL$Sensors
    Src <- InL$Sources
    Interval <- InL$Interval
    Out <- Source$bLSout
    WDmean <- (360 + atan2(sum(sin(Interval[,"WD [deg N]"]/180*pi)*Interval[,"Ustar [m/s]"])/sum(Interval[,"Ustar [m/s]"]),sum(cos(Interval[,"WD [deg N]"]/180*pi)*Interval[,"Ustar [m/s]"])/sum(Interval[,"Ustar [m/s]"]))/pi*180) %% 360
    if(missing(stats)){
      stats <- Source$stats
    } else {
        Sindex <- expand.grid(seq.int(length(xm)),seq.int(length(ym)))
        SensorNames <- paste0("Sensor_",apply(Sindex,1,paste,collapse="_"))
        xy[] <- Out[SensorNames, .(V1 = eval(parse(text=stats)))][,V1]
    }
    coord2wgs84 <- Source$coord2wgs84
    MyMap <- Source$MyMap
    Source <- Src
    if(!add){
      if(is.null(MyMap)){
        do.call(siteMap,c(list(Source,add=add,xlim=xlim,ylim=ylim),siteMapArgs))
        abline(h=pretty(xlim),lty="dotted",col="lightgrey")
        abline(v=pretty(ylim),lty="dotted",col="lightgrey")
      } else {
        do.call(RgoogleMaps::PlotOnStaticMap,
          c(list(MyMap = MyMap), StaticMapArgs))
      }
      usr <- par()$usr
      xlim <- usr[1:2]
      ylim <- usr[3:4]
    }
    xlimOriginal <- xlim
    ylimOriginal <- ylim
    if(!is.null(coord2wgs84)){
      lim1 <-  optim(colMeans(Source[,2:3]),function(x)sum((unlist(coord2wgs84(x=x[1],y=x[2])) - c(xlim[1],ylim[1]))^2))$par
      lim2 <-  optim(colMeans(Source[,2:3]),function(x)sum((unlist(coord2wgs84(x=x[1],y=x[2])) - c(xlim[2],ylim[2]))^2))$par
      xlim <- c(lim1[1],lim2[1])
      ylim <- c(lim1[2],lim2[2])
    }
  }

  xy[is.na(xy)] <- 0

  if(!is.function(levels)){
    sub_levels <- eval(substitute(levels))
    levels <- function(x)sub_levels
  }

  brks <- levels(xy)
  cl <- contourLines(c(xm[1]-dx,xm,xm[length(xm)]+dx),c(ym[1]-dy,ym,ym[length(ym)]+dy),rbind(0,cbind(0,xy,0),0),levels=brks)
  clFac <- as.factor(sapply(cl,"[[","level"))
  uclFac <- unique(clFac)

  # transform x/y:
  if(!is.null(coord2wgs84)){
    cl <- lapply(cl,function(x){
      # a <- do.call(coord2wgs84,list(x=x$x,y=x$y))
      a <- coord2wgs84(x=x$x,y=x$y)
      list(level=x$level,x=a[[1]],y=a[[2]])
    })
  }
  if(addSource){
    for(i in attr(Source,"SourceList")){
      pg <- i
      if(!is.null(coord2wgs84)){
        pg[1:2] <- coord2wgs84(x=pg[,1],y=pg[,2])
      }
      upg <- unique(pg[,3])
      for(j in upg){
        ind <- which(pg[,3]==j)
        polygon(pg[ind,1],pg[ind,2],col="#006837")
      }
    }
  }
  if(fill){
    clP <- lapply(cl,function(x)sp::Polygon(cbind(x$x,x$y)))
    clPs <- lapply(uclFac,function(x,y,z,a,b)maptools::checkPolygonsHoles(sp::Polygons(y[z==x],as.character(x)),useSTRtree=a,avoidGEOS=b),y=clP,z=clFac,a=useSTRtree,b=avoidGEOS)
    clSp <- sp::SpatialPolygons(clPs, as.integer(uclFac))
    
    # Loesung:
    mlt <- switch(axs[1],
      "i"=1.08,
      1
      )
    xlim2 <- xlimOriginal - diff(xlimOriginal)*(1-1/mlt)/2;ylim2 <- ylimOriginal - diff(ylimOriginal)*(1-1/mlt)/2
    c1 <- sp::Polygon(cbind(c(xlim2[1],xlim2[2],xlim2[2],xlim2[1], xlim2[1]),c(ylim2[1],ylim2[1],ylim2[2],ylim2[2],ylim2[1])))
    c2 <- sp::Polygons(list(c1), "sclip")
    Pclip <- sp::SpatialPolygons(list(c2))

    SpPclip <- rgeos::gIntersection(clSp, Pclip, byid = TRUE)
    
    alphachar <- as.hexmode(round(alpha*255))

    if(is.null(col)){
      col <- ConcPalette(length(brks))
    } else if(is.function(col)){
      col <- col(length(brks))
    }
    if(is.null(lty))lty <- 1
    if(is.null(lwd))lwd <- 1
    fillcol <- paste0(col,alphachar)
    plot(SpPclip,col=fillcol,border=col,lty=lty,lwd=lwd,add=TRUE)
  } else {
    if(is.null(col)){
      col <- rev(grey.colors(length(brks),start=0,end=0.6))
    } else if(is.function(col)){
      col <- col(length(brks))
    }
    if(is.null(lty)){
      lty <- 1:8
      lty <- rev(rep(lty,length(brks))[seq_along(brks)])
    }
    if(is.null(lwd))lwd <- 1
    col <- rep(col,length(brks))[seq_along(brks)]
    lty <- rep(lty,length(brks))[seq_along(brks)]
    lwd <- rep(lwd,length(brks))[seq_along(brks)]
    for(i in seq_along(cl)){
      ind <- as.integer(clFac[i])
      lines(cl[[i]]$x,cl[[i]]$y,col=col[ind],lty=lty[ind],lwd=lwd[ind])
    } 
  }
  box()


  scale <- if(is.null(MyMap)) 1 else MyMap
  if(addSB)addScaleBar(SBpos,scale)
  if(addWR)addWindrose(WDmean,pos=WRpos,frac=WRfrac,scF=WRscale)

  if(showLegend){ 
    if(is.null(lpos))lpos <- c("bottomleft","topleft","topright","bottomright")[(floor(WDmean/90)-2)%%4+1]
    ltex <- signif(rev(brks),sigNums)
    if(is.null(leg.bg.col)&&!is.null(MyMap))leg.bg.col <- "white"
    if(fill){
      ltex <- paste0(">",ltex)
      legend(lpos,legend=ltex,xjust=0,fill=rev(fillcol),border=rev(col),bg=leg.bg.col)
    } else {
      legend(lpos,legend=ltex,xjust=0,col=rev(col),lty=rev(lty),lwd=rev(lwd),bg=leg.bg.col)
    }
  }
  # if(depo){
  #   tmp <- copy(y)
  #   y <- copy(x)
  #   x <- tmp
  # }
  return(invisible(structure(
    list(Img=list(x=xm,y=ym,z=xy),stats=stats,limits=list(xlim=xlimOriginal,ylim=ylimOriginal),bLSout=Out,bLSin=InL,MyMap=MyMap,coord2wgs84=coord2wgs84)
    ,class=c("contourXY","list")
    )))
}


calcDep <- function(inlist, Cat.PathIn, use.avg, use.sym,
  quiet, SrangeIn, ScalcIn, SensorsIn){
  # get input from list
  IntIn <- inlist[[1]]
  SensorNamesIn <- inlist[[2]]
  outOut <- vector("list", nrow(IntIn))
  Cat.Name <- ""
  Cat.N0 <- 0
  for (i in seq.int(nrow(IntIn))){
    N0 <- IntIn[i, N0]
    C0 <- data.table(Traj=1L:IntIn[i,N0],Conc=0)
    setkey(C0,"Traj")
    outOut[[i]] <- setkey(copy(inlist[[3]]), Sensor)
    if(Cat.Name!=IntIn[i,Cat.Name] || Cat.N0 != N0){
      Cat.Name <- IntIn[i,Cat.Name]
      Ctlg <- readCatalog(paste(Cat.PathIn,Cat.Name,sep="/"))
      Cat.Ustar <- attr(Ctlg, "Ustar")
      Cat.WD <- IntIn[i,WD]
      uvw <- uvw0(Ctlg)
      if(use.avg){
        Ctlg[,wTD:=1/(mean(1/wTD))]
        # Ctlg[,wTD:=mean(wTD)]
      }
      # subset?
      if(attr(Ctlg,"N0")>N0){
        seedN <- sample.int(1E9,1)
        env <- globalenv()
        oseed <- env$.Random.seed
        set.seed(seedN,kind="L'Ecuyer-CMRG")      
        takeSub <- sort(sample.int(attr(Ctlg,"N0"),N0))
        if (is.null(oseed)) {
            rm(list = ".Random.seed", envir = env)
        } else {
            assign(".Random.seed", value = oseed, envir = env)
        }
        attCat <- attributes(Ctlg)
        indexNew <- 1:N0
        names(indexNew) <- takeSub
        Ctlg <- Ctlg[Traj_ID %in% takeSub,]
        Ctlg[,":="(Traj_ID=indexNew[as.character(Traj_ID)])]
        for(ac in (names(attCat) %w/o% c("names","row.names",".internal.selfref","uvw0")))setattr(Ctlg,ac,attCat[[ac]])
        uvw <- uvw[takeSub,]
        setattr(Ctlg,"uvw0",uvw)
        HeaderSub <- paste(attr(Ctlg,"header"),"\n*** Subset of ",N0," Trajectories ***\n\n",sep="")
        class(HeaderSub) <- c("TDhead","character")
        setattr(Ctlg,"N0",N0)
        setattr(Ctlg,"header",HeaderSub)
      }
      Cat.N0 <- attr(Ctlg,"N0")
      if(use.sym){
        attCat <- attributes(Ctlg)
        Ctlg[,wTD:=wTD*2]
        Ctlg <- rbind(Ctlg,
          Ctlg[,.(
            Traj_ID=Traj_ID,
            Time=Time,
            x=x,
            y=-y,
            wTD=wTD
            )]
          )
        for(ac in (names(attCat) %w/o% c("names","row.names",".internal.selfref","uvw0")))setattr(Ctlg,ac,attCat[[ac]])
        uvw[,2] <- 0 
        setattr(Ctlg,"uvw0",uvw)
        setkey(Ctlg,"Traj_ID")
      }
      rotateCatalog(Ctlg,Cat.WD)
    } else if(Cat.WD!=IntIn[i,WD]){
      rotateCatalog(Ctlg,Cat.WD,back=TRUE)
      Cat.WD <- IntIn[i,WD]
      rotateCatalog(Ctlg,Cat.WD)
    }
    if(Cat.Ustar!=IntIn[i,Ustar]){
      ## korrigiere Catalog wTD uvw0 U0:
      initializeCatalog(IntIn[i,],Catalog=Ctlg)
      uvw <- uvw0(Ctlg)
      Cat.Ustar <- attr(Ctlg, "Ustar")
    }
    Catalog <- copy(Ctlg)

    if(!quiet)cat("\nGet TD inside source areas:\n")
    bLSmodelR:::tagNear(Catalog,SrangeIn)
    Catalog[,inside0:=inside]
    Catalog[,rn:=.I]
    if(Catalog[,any(inside)]){
      for(j in SensorNamesIn){
        if(!quiet)cat("Calculate Tile",which(j==SensorNamesIn),"/",length(SensorNamesIn),"\r")
        SourceAreaRelative <- copy(ScalcIn)[,":="(x=x-SensorsIn[j,x],y=y-SensorsIn[j,y])]
        Catalog[,inside:=inside0]
        bLSmodelR:::tagNear(Catalog,SourceAreaRelative)
        Catalog[,inside1:=inside]

        if(Catalog[,any(inside)]){
          # tag Inside Source
          TDinside <- SourceAreaRelative[,
          {
            bLSmodelR:::tagNear(Catalog[,inside:=inside1],.(x=x,y=y))
            cbind(ID=Catalog[(inside),rn],bLSmodelR::pnt.in.poly(Catalog[(inside),cbind(x,y)],cbind(x,y)))
          },by=pid][,sum(pip),by=ID]
          if(any(TDinside[,as.logical(V1)])){
            setkey(Catalog,rn)
            Catalog[TDinside,inside:=as.logical(V1)]

            # calc CE
            Ci <- Catalog[(inside), sum(2 / wTD), by = Traj_ID][
              C0, on = c(Traj_ID = "Traj")][is.na(V1), V1 := 0]
            Ci[, ":="(
              CE = V1,
              uCE = V1 * uvw[, "u0"],
              vCE = V1 * uvw[, "v0"],
              wCE = V1 * uvw[, "w0"]
              )]

            # Max_Dist etc.
            setkey(Catalog,Traj_ID)
            Cat <- Catalog[Catalog[(inside),.(minTime=min(Time)),by=Traj_ID]][Time>=minTime,]
            outOut[[i]][j,":="(
              CE=Ci[,mean(CE)],CE_se=Ci[,sqrt(var(CE)/N0)],
              uCE=Ci[,mean(uCE)],uCE_se=Ci[,sqrt(var(uCE)/N0)],
              vCE=Ci[,mean(vCE)],vCE_se=Ci[,sqrt(var(vCE)/N0)],
              wCE=Ci[,mean(wCE)],wCE_se=Ci[,sqrt(var(wCE)/N0)],
              Max_Dist = max(Max_Dist,rotateCatalog(Cat,IntIn[i,WD],back=TRUE)[,-min(x)],na.rm=TRUE),
              N_TD = N_TD + Catalog[,sum(inside)],
              TD_Time_avg = sum(TD_Time_avg,Catalog[(inside),-mean(Time)],na.rm=TRUE),
              TD_Time_max = max(TD_Time_max,Catalog[(inside),-min(Time)],na.rm=TRUE)
              )]
          }
        }

      }
    }
    outOut[[i]][,":="(
      CE_lo=CE + qt(0.025,N0-1)*CE_se,CE_hi=CE + qt(0.975,N0-1)*CE_se,
      uCE_lo=uCE + qt(0.025,N0-1)*uCE_se,uCE_hi=uCE + qt(0.975,N0-1)*uCE_se,
      vCE_lo=vCE + qt(0.025,N0-1)*vCE_se,vCE_hi=vCE + qt(0.975,N0-1)*vCE_se,
      wCE_lo=wCE + qt(0.025,N0-1)*wCE_se,wCE_hi=wCE + qt(0.975,N0-1)*wCE_se
    )]
  }
  outOut
}



# -> transform to xy (oder ist schon xy?)

# rotate into WD
# get MaxFetch (delta_x) & delta_y <-> xlim, ylim
# calc bLS (via runbLS?), oder bLSmodelR:::.calcCatalogs & 
#   bLSmodelR:::.calcCE -> parallel?!
# transform

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

  source("/home/christoph/repos/2_Packages/1_bLSmodelR/Functions Aufraeumen/ToDo/contourXY.r")

  # XY_out_17 <- contourXY(Sources[Sources[,1] == "EVS",], z.meas = 1.54, MyMap=EVS_Map,
  #   coord2wgs84 = function(x,y){
  #     # browser()
  #     list(lat=CH.to.WGS.lat(x,y),lon=CH.to.WGS.lng(x,y))},
  #   Interval = Snc, Cat.Path = CatPath, dx = 5, dy = 5,
  #   fill = TRUE,ncores=ncores)

  x11(width = 10, height = 10)
  par(mfrow=c(3,4))
  for(i in 1:10){
    XY_out_17 <- contourXY(Sources[Sources[,1] == "EVS",], z.meas = 1.54, MyMap=EVS_Map,
      coord2wgs84 = function(x,y){
        # browser()
        list(lat=CH.to.WGS.lat(x,y),lon=CH.to.WGS.lng(x,y))},
      Interval = Snc[i,], Cat.Path = CatPath, dx = 5, dy = 5,
      levels = c(0.05, 0.1, 0.2, 0.5, 1, 2, 5),
      fill = TRUE,ncores=ncores)
  }
  XY_out_17 <- contourXY(Sources[Sources[,1] == "EVS",], z.meas = 1.54, MyMap=EVS_Map,
    coord2wgs84 = function(x,y){
      # browser()
      list(lat=CH.to.WGS.lat(x,y),lon=CH.to.WGS.lng(x,y))},
    Interval = Snc, Cat.Path = CatPath, dx = 5, dy = 5,
    levels = c(0.05, 0.1, 0.2, 0.5, 1, 2, 5),
    fill = TRUE,ncores=ncores)

  contourXY(XY_out_17, fill = TRUE, levels = c(0.05, 0.1, 0.2, 0.5, 1, 2, 5))
  # debug(initializeCatalog)
  # undebug(initializeCatalog)

# -> browser() bei Mittelung einbauen (und bei rotateCatalog?)

# -> bLSmodelR::rotateCatalog stimmt überhaupt nicht. Nochmals in Ruhe überdenken!




}
