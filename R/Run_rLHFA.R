#' run rLakeHeatFluxAnalyzer
#'
#' @param LakeName Name of the Lake. All forcing and config files must have this name
#' @param Folder Folder the forcing and config files are stored in
#' @param skipLoad Use own config file (must then be placed in Folder)
#' @export
#' @examples Run_LHFA(LakeName="Esthwaite",Folder="../data",skipLoad=TRUE)

Run_LHFA <- function(LakeName,Folder,skipLoad = FALSE){
#----Author: Jordan S Read 2009 ----
  #----Modified by R. Iestyn Woolway ----
  #---- R translation by Johannes Feldbauer -----

  graphics.off()
  cat("\f")

  lhfa_version <- '1.1.2'

  #done <- FALSE
  if (!skipLoad){
    build_config(LakeName,Folder)
  }

  cat(paste0('Running LHFA version ', lhfa_version, '\n\n'))

  # -- variables --
  matSec <- 86400         # number of seconds in a day
  smplTs <- 50            # samplerate test length
  dateTl <- 0.00001       # tolerance for equal day (day fraction)
  cat(paste0('Reading ', LakeName, '.hfx file'))
  Cfg <- OpenCfg( LakeName,Folder )
  # -- end variables --

  if (!any(c(Cfg$plotYes, Cfg$writeYes))){
    stop(paste0('User must specify either to write results,',
                ' plot results, or both . (see ', LakeName, '.hfx)'))
  }

  cat('completed\n\n')

  cat('****Building program structure****\n')
  pltMods <- NULL
  if (Cfg$plotYes){
    cat(paste0('Checking for ', LakeName, '.plt file*'))
    pltFileName  <- paste0(Folder, '/', LakeName, '.plt')
    oper <- file.exists(pltFileName)
    if (!oper){
      cat('not found [*is optional]\n\n')
    } else {

      pltMods <- read.table(pltFileName)
      cat('completed\n')
    }
  }
  OC <- OutputConstructor(Cfg$outPuts,Cfg$pltMods)
  fNms <- names(OC$TT)
  for (j in 1:length(fNms)){
    if (OC$TT[[fNms[j]]]){
      cat(paste0(fNms[j], '\n'))
    }
  }
  cat('****completed****\n\n')

  if (OC$TT$openWtr){
    cat(paste0('Reading ', LakeName, '.wtr file \n'))
    wtrFileName  <- paste0(Folder, '/', LakeName, '.wtr')
    oper <- file.exists(wtrFileName)
    if (!oper){
      stop(paste0(LakeName, '.wtr file not found \n'))
    }

    wt <- gFileOpen(wtrFileName,treatAsWtr = TRUE)
    #headers <- textscan(wt$depths,'#s','Delimiter','\t')
    # headers <- headers{1}(2:end) # get rid of dateTime
    depths <- wt$depths
    # for d =1:length(depths)
    # txt <- headers{d}
    # depths(d) <- str2double(txt(5:end))
    # end
    mnDep <- min(depths)
    mnIdx <- which.min(depths)
    if (length(depths) > 1){
      wt$dat <- wt$dat[,mnIdx]
    }
    if (mnDep>0){
      cat(paste0(' ', mnDep, ' m is the shallowest depth in the .wtr file',
                 ' which will be used to represent surface water temperatures. \n'))
    }

    # remove nans
    idx <- is.na(wt$dat)
    wt$dat <- wt$dat[!idx]
    wt$dates <- wt$dates[!idx]
    cat('completed\n\n')

    #*** find samplerate of raw data
    if(length(wt$dates) < smplTs){
      tLen <- length(wt$dates)
    } else {
      tLen <- smplTs
    }
    steps <- rep(NA,tLen-1)
    for (i in 1:(tLen-1)){
      steps[i] <- wt$dates[i+1] - wt$dates[i]
    }
    if(min(steps)==0){
      matRs <- mean(steps)*matSec
    } else {
     matRs <- min(steps)*matSec #current sample rate of raw data in seconds
    }

    if ((Cfg$outRs - matRs)>dateTl){
      OC$TT$dwnSmple <- TRUE   #down sample if necessary
    }

  # *** down sampling ***
    if (OC$TT$dwnSmple){
      cat(paste0('Down sampling ', LakeName, '.wtr data'))
      DSwt <- DownSample_TS(wt$dates,Cfg$outRs,wt$dat)
      DS_wtrD <- DSwt$DS_dates
      DS_wtr <- DSwt$DS_varArray
      wt$dates <- DS_wtrD
      wt$dat <- DS_wtr
      varL <- length(wt$dates)
      rm (DS_wtrD, DS_wtr,DSwt)
      cat('completed\n\n')
    }
  }

  if (OC$TT$openWnd){
    cat(paste0('Reading ', LakeName, '.wnd file \n'))
    wndFileName  <- paste0(Folder, '/', LakeName, '.wnd')
    oper <- file.exists(wndFileName)
    if (!oper){
      stop(paste0(LakeName, '.wnd file not found \n'))
    }
    wnd <- gFileOpen(wndFileName)
    # set values between min and max values
    wnd$dat[wnd$dat < Cfg$wndMn] <- Cfg$wndMn
    wnd$dat[wnd$dat > Cfg$wndMx] <- Cfg$wndMx

    # remove nans
    idx <- is.na(wnd$dat)
    wnd$dat <- wnd$dat[!idx]
    wnd$dates <- wnd$dates[!idx]
    cat('completed\n\n')

    #*** find samplerate of raw data
    if (length(wnd$dates) < smplTs){
      tLen <- length(wnd$dates)
    } else{
      tLen <- smplTs
    }
    steps <- rep(NA,tLen-1)
    for (i in 1:(tLen-1)){
      steps[i] <- wnd$dates[i+1]-wnd$dates[i]
    }
    if (min(steps)==0){
      matRs <- mean(steps)*matSec
    } else {
      matRs <- min(steps)*matSec #current sample rate of raw data in seconds
    }
    rm( idx, steps, tLen)
    if ((Cfg$outRs - matRs)>dateTl){
      OC$TT$dwnSmple <- TRUE   #down sample if necessary
    }

    # *** down sampling ***
    if (OC$TT$dwnSmple){
      cat(paste0('Down sampling ', LakeName, '.wnd data \n'))
      DS_wnd <- DownSample_TS(wnd$dates,Cfg$outRs,wnd$dat)
      wnd$dates <- DS_wnd$DS_dates
      wnd$dat <- DS_wnd$DS_varArray
      varL <- length(wnd$dates)
      rm(DS_wnd)
      cat('completed\n\n')
    }
    dates <- wnd$dates
  }

  if (OC$TT$openSW){
    if (file.exists(paste0(Folder, '/', LakeName, '.sw'))){
      cat(paste0('Reading ', LakeName, '.sw file'))
      swFileName  <- paste0(Folder, '/', LakeName, '.sw')
      sw <- gFileOpen(swFileName)
      # remove negative values
      sw$dat[sw$dat < 0] <- 0

      # remove nans
      idx <- is.na(sw$dat)
      sw$dat <- sw$dat[!idx]
      sw$dates <- sw$dates[!idx]
      cat('completed\n\n')

  #*** find samplerate of raw data
      if (length(sw$dates) < smplTs){
        tLen <- length(sw$dates)
      } else {
        tLen <- smplTs
      }
      steps <- rep(NA,tLen-1)
      for (i in 1:(tLen-1)){
        steps[i] <- sw$dates[i+1] - sw$dates[i]
      }
      if (min(steps)==0){
        matRs <- mean(steps)*matSec
      } else {
        matRs <- min(steps)*matSec #current sample rate of raw data in seconds
      }
      rm (idx, steps, tLen)
      if((Cfg$outRs - matRs)>dateTl){
        OC$TT$dwnSmple <- TRUE   #down sample if necessary
      }

  # *** down sampling ***
      if (OC$TT$dwnSmple){
        cat(paste0('Down sampling ', LakeName, '.sw data \n'))
        DS_sw <- DownSample_TS(sw$dates,Cfg$outRs,sw$dat)
        sw$dates <- DS_sw$DS_dates
        sw$dat <- DS_sw$DS_varArray
        varL <- length(sw$dates)
        rm (DS_sw)
        cat('completed\n\n')
      }

    } else
      if (file.exists(paste0(Folder, '/', LakeName, '.par'))){
        cat(paste0('Reading ', LakeName, '.par file \n'))
        swFileName  <- paste0(Folder, '/', LakeName, '.par')
        sw <- gFileOpen(swFileName)
        sw$dat[sw$dat < 0] <- 0
        # remove nans
        idx <- is.na(sw$dat)
        sw$dat <- sw$dat[!idx]
        sw$dates <- sw$dates[!idx]
        parMult <- 0.4957
        sw$dat <- sw$dat*parMult
        cat('completed\n\n')

  #*** find samplerate of raw data
        if (length(sw$dat) < smplTs){
          tLen <- length(sw$dates)
        } else {
          tLen <- smplTs
        }
        steps <- rep(NA,tLen-1)
        for (i in 1:(tLen-1)){
          steps[i] <- sw$dates[i+1]-sw$dates[i]
        }
        if (min(steps)==0){
          matRs <- mean(steps)*matSec
        } else {
          matRs <- min(steps)*matSec #current sample rate of raw data in seconds
        }
        rm (idx, steps, tLen)
        if ((Cfg$outRs - matRs)>dateTl){
         OC$TT$dwnSmple <- TRUE   #down sample if necessary
        }

    # *** down sampling ***
        if (OC$TT$dwnSmple){
          cat(paste0('Down sampling ', LakeName, '.sw data \n'))
          DS_sw <- DownSample_TS(sw$dates,Cfg$outRs,sw$dat)
          sw$dates <- DS_sw$DS_dates
          sw$dat <- DS_sw$DS_varArray
          varL <- length(sw$dates)
          rm (DS_sw)
          cat('completed\n\n')
        }

     } else {
      stop(paste0(LakeName, '.sw nor .par file not found \n'))
    }

    dates <- sw$dates
  }

  if (OC$TT$openAirT){
    cat(paste0('Reading ', LakeName, '.airT file \n'))
    airTFileName  <- paste0(Folder, '/', LakeName, '.airT')
    oper <- file.exists(airTFileName)
    if (!oper){
      stop(paste0(LakeName, '.airT file not found'))
    }

    airT <- gFileOpen(airTFileName)

    # remove nans
    idx <- is.na(airT$dat)
    airT$dat <- airT$dat[!idx]
    airT$dates <- airT$dates[!idx]
    cat('completed\n\n')

    #*** find samplerate of raw data
    if (length(airT$dates) < smplTs){
      tLen <- length(airT$dates)
    } else {
      tLen <- smplTs
    }
    steps <- rep(NA,tLen-1)
    for (i in 1:(tLen-1)){
      steps[i] <- airT$dates[i+1]-airT$dates[i]
    }
    if (min(steps)==0){
      matRs <- mean(steps)*matSec
    } else {
      matRs <- min(steps)*matSec #current sample rate of raw data in seconds
    }
    rm( idx,steps, tLen)
    if ((Cfg$outRs - matRs)>dateTl){
      OC$TT$dwnSmple <- TRUE   #down sample if necessary
    }

    # *** down sampling ***
    if (OC$TT$dwnSmple){
      cat(paste0('Down sampling ', LakeName, '.airT data \n'))
      DS_airT <- DownSample_TS(airT$dates,Cfg$outRs,airT$dat)
      airT$dates <- DS_airT$DS_dates
      airT$dat <- DS_airT$DS_varArray
      varL <- length(airT$dates)
      rm( DS_airT)
      cat('completed\n\n')
    }
  }

  if (OC$TT$openRH){
    cat(paste0('Reading ', LakeName, '.rh file \n'))
    rhFileName  <- paste0(Folder, '/', LakeName, '.rh')
    oper <- file.exists(rhFileName)
    if (!oper){
      stop(paste0(LakeName, '.rh file not found'))
    }

    rh <- gFileOpen(rhFileName)
    rh$dat[rh$dat < 0] <- 0
    rh$dat[rh$dat > 100] <- 100

    # remove nans
    idx <- is.na(rh$dat)
    rh$dat <- rh$dat[!idx]
    rh$dates <- rh$dates[!idx]
    cat('completed\n\n')

    #*** find samplerate of raw data
    if (length(rh$dates) < smplTs){
      tLen <- length(rh$dates)
    } else {
      tLen <- smplTs
    }
    steps <- rep(NA,tLen-1)
    for (i in 1:(tLen-1)){
      steps[i] <- rh$dates[i+1]-rh$dates[i]
    }
    if(min(steps)==0){
      matRs <- mean(steps)*matSec
    } else {
      matRs <- min(steps)*matSec #current sample rate of raw data in seconds
    }
    rm( idx,steps, tLen)
    if ((Cfg$outRs - matRs)>dateTl){
      OC$TT$dwnSmple <- TRUE   #down sample if necessary
    }

    # *** down sampling ***
    if( OC$TT$dwnSmple){
      cat(paste0('Down sampling ', LakeName, '.rh data \n'))
      DS_rh <- DownSample_TS(rh$dates,Cfg$outRs,rh$dat)
      rh$dates <- DS_rh$DS_dates
      rh$dat <- DS_rh$DS_varArray
      varL <- length(rh$dates)
      rm(DS_rh)
      cat('completed\n\n')
    }
  }

  # make sure dates match
  if (!OC$TT$openSW){
    if (OC$TT$openWtr && OC$TT$openWnd && OC$TT$openRH && OC$TT$openAirT){
      idx <- intersect(intersect(intersect(wt$dates,wnd$dates),rh$dates),airT$dates)
      dates <- as.POSIXct(idx,origin = "1970-01-01")
      varL <- length(dates)
      wt$dat <- wt$dat[is.element(wt$dates,idx)]
      wnd$dat <- wnd$dat[is.element(wnd$dates,idx)]
      rh$dat <- rh$dat[is.element(rh$dat,idx)]
      airT$dat <- airT$dat[is.element(airT$dates,idx)]
      wt$dates <- as.POSIXct(idx,origin = "1970-01-01")
      wnd$dates <- as.POSIXct(idx,origin = "1970-01-01")
      rh$dates <- as.POSIXct(idx,origin = "1970-01-01")
      airT$dates <- as.POSIXct(idx,origin = "1970-01-01")
    }
  }

  if (!OC$TT$openWnd){
    if (OC$TT$openRH && OC$TT$openAirT && OC$TT$openWtr && OC$TT$openSW){
      idx <- intersect(intersect(intersect(rh$dates,airT$dates),wt$dates),sw$dates)
      dates <- as.POSIXct(idx,origin = "1970-01-01")
      varL <- length(dates)
      rh$dat <- rh$dat[is.element(rh$dates,idx)]
      airT$dat <- airT$dat[is.element(airT$dates,idx)]
      wt$dat <- wt$dat[is.element(wt$dates,idx)]
      sw$dat <- sw$dat[is.element(sw$dates,idx)]
      wt$dates <- as.POSIXct(idx,origin = "1970-01-01")
      sw$dates <- as.POSIXct(idx,origin = "1970-01-01")
      rh$dates <- as.POSIXct(idx,origin = "1970-01-01")
      airT$dates <- as.POSIXct(idx,origin = "1970-01-01")
    }
  }

  if (OC$TT$openRH && OC$TT$openAirT && OC$TT$openWtr && OC$TT$openSW && OC$TT$openWnd){
    idx <- intersect(intersect(intersect(intersect(rh$dates,airT$dates),wt$dates),
                               sw$dates),wnd$dates)
    dates <- as.POSIXct(idx,origin = "1970-01-01")
    varL <- length(dates)
    rh$dat <- rh$dat[is.element(rh$dates,idx)]
    airT$dat <- airT$dat[is.element(airT$dates,idx)]
    wt$dat <- wt$dat[is.element(wt$dates,idx)]
    sw$dat <- sw$dat[is.element(sw$dates,idx)]
    wnd$dat <- wnd$dat[is.element(wnd$dates,idx)]
    wt$dates <- as.POSIXct(idx,origin = "1970-01-01")
    sw$dates <- as.POSIXct(idx,origin = "1970-01-01")
    rh$dates <- as.POSIXct(idx,origin = "1970-01-01")
    airT$dates <- as.POSIXct(idx,origin = "1970-01-01")
    wnd$dates <- as.POSIXct(idx,origin = "1970-01-01")
  }

  if (!OC$TT$openRH && !OC$TT$openAirT && !OC$TT$openSW && !OC$TT$openWnd){
    if (OC$TT$openWtr){
      dates <- wt$dates
      varL <- length(dates)
    }
  }

  # look for long-wave radiation data
  if (OC$TT$openLWnet){
    if (file.exists(paste0(Folder, '/', LakeName, '.lwnet'))){
      cat(paste0('Reading ', LakeName, '.lwnet file \n'))
      lwnetFileName <- paste0(Folder, '/', LakeName, '.lwnet')

      lwnet <- gFileOpen(lwnetFileName)
      idx <- is.na(lwnet$dat)
      lwnet$dat <- lwnet$dat[!idx]
      lwnet$dates <- lwnet$dates[!idx]
      cat('completed\n\n')

      #*** find samplerate of raw data
      if (length(lwnet$dates) < smplTs){
        tLen <- length(lwnet$dates)
      } else {
        tLen <- smplTs
      }
      steps <- rep(NA,tLen-1)
      for (i in 1:(tLen-1)){
        steps[i] <- lwnet$dates[i+1]-lwnet$dates[i]
      }
      if (min(steps)==0){
        matRs <- mean(steps)*matSec
      } else {
        matRs <- min(steps)*matSec #current sample rate of raw data in seconds
      }
      rm(idx,steps, tLen)
      if ((Cfg$outRs - matRs)>dateTl){
        OC$TT$dwnSmple <- TRUE   #down sample if necessary
      }

      # *** down sampling ***
      if (OC$TT$dwnSmple){
        cat(paste0('Down sampling ', LakeName, '.lwnet data \n'))
        DS_lwnet <- DownSample_TS(lwnet$dates,Cfg$outRs,lwnet$dat)
        lwnet$dates <- DS_lwnet$DS_dates
        lwnet$dat <- DS_lwnet$DS_varArray
        varL <- length(lwnet$dates)
        rm(DS_lwnet)
        cat('completed\n\n')
      }

    } else {
        if (file.exists(paste0(Folder, '/', LakeName, '.lw'))){
          cat(paste0(LakeName, '.lwnet files not found, looking for .lw file\n\n'))
          cat(paste0('Reading ', LakeName, '.lw file \n'))
          lwFileName  <- paste(Folder, '/', LakeName, '.lw')

          lw <- gFileOpen(lwFileName)
          idx <- is.na(lw$dat)
          lw$dat <- lw$dat[!idx]
          lw$dates <- lw$dates[!idx]
          cat('completed\n\n')

          #*** find samplerate of raw data
          if (length(lw$dates) < smplTs){
            tLen <- length(lw$dates)
          } else {
            tLen <- smplTs
          }
          steps <- rep(NA,tLen-1)
          for (i in 1:(tLen-1)){
            steps[i] <- lw$dates[i+1]-lw$dates[i]
          }
          if (min(steps)==0){
            matRs <- mean(steps)*matSec
          } else {
            matRs <- min(steps)*matSec #current sample rate of raw data in seconds
          }
          rm(idx, steps, tLen)
          if ((Cfg$outRs - matRs)>dateTl){
            OC$TT$dwnSmple <- TRUE   #down sample if necessary
          }

          # *** down sampling ***
          if (OC$TT$dwnSmple){
            cat(paste0('Down sampling ', LakeName, '.lw data \n'))
            DS_lw <- DownSample_TS(lw$dates,Cfg$outRs,lw$dat)
            lw$dates <- DS_lwD$DS_dates
            lw <- DS_lw$DS_varArray
            varL <- length(lw$dates)
            rm(DS_lw)
            cat('completed\n\n')
          }

          # find when wtr and lw dates intersect
          idx <- intersect(wt$dates,lw$dates)
          lw$dat <- lw$dat[is.element(lw$dates,idx)]
          wt$dat <- wt$dat[is.element(wt$dates,idx)]
          wt$dates <- as.POSIXct(idx,origin = "1970-01-01")

          Tk <- wt$dat + 273.13 # .wtr already called at this point
          emiss <- 0.972
          S_B <- 5.67E-8
          LWo <- S_B*emiss*Tk^4

          # define lwnet
          lwnet <- list()
          lwnet$dat <- lw$dat - LWo
          lwnet$dates <- as.POSIXct(idx,origin = "1970-01-01")
        } else {
          cat(paste0('',LakeName,
                     '.lwnet and .lw file missing, using .airT, .rh, .sw, .wtr instead\n\n'))
          press <- 101325*(1 - 2.25577E-5*Cfg$alt)^5.25588 # Pa
          press <- press/100 # mb
          Qlnet <- calc_lwnet(dates,Cfg$lat,press,airT$dat,rh$dat,sw$dat,wt$dat)$lwnet
          lwnet <- list()
          lwnet$dat <- Qlnet
          lwnet$dates <- dates
        }
    }
  }

  # re-adjust dates depending on lw and lwnet files
  if (OC$TT$openRH && OC$TT$openAirT && OC$TT$openWtr && OC$TT$openSW && OC$TT$openLWnet){
    idx <- intersect(intersect(intersect(intersect(rh$dates,airT$dates),
                                         wt$dates),sw$dates),lwnet$dates)
    dates <- as.POSIXct(idx,origin="1970-01-01")
    varL <- length(dates)
    rh$dat <- rh$dat[is.element(rh$dates,idx)]
    airT$dat <- airT$dat[is.element(airT$dates,idx)]
    wt$dat <- wt$dat[is.element(wt$dates,idx)]
    sw$dat <- sw$dat[is.element(sw$dates,idx)]
    lwnet$dat <- lwnet$dat[is.element(lwnet$dates,idx)]
    wt$dates <- as.POSIXct(idx,origin="1970-01-01")
    sw$dates <- as.POSIXct(idx,origin="1970-01-01")
    rh$dates <- as.POSIXct(idx,origin="1970-01-01")
    airT$dates <- as.POSIXct(idx,origin="1970-01-01")
    lwnet$dates <- as.POSIXct(idx,origin="1970-01-01")
  }

  if( OC$TT$openRH && OC$TT$openAirT && OC$TT$openWtr && OC$TT$openSW &&
      OC$TT$openLWnet && OC$TT$openWnd){
    idx <- intersect(intersect(intersect(intersect(intersect(rh$dates,airT$dates),wt$dates),
                                         sw$dates),lwnet$dates),wnd$dates)
    dates <- idx
    varL <- length(dates)
    rh$dat <- rh$dat[is.element(rh$dates,idx)]
    airT <- airT[is.element(airT$dates,idx)]
    wt$dat <- wt$dat[is.element(wt$dates,idx)]
    sw$dat <- sw$dat[is.element(sw$dates,idx)]
    lwnet$dat <- lwnet$dat[is.element(lwnet$dates,idx)]
    wnd$dat <- wnd$dat[is.element(wnd$dates,idx)]
    wt$dates <- as.POSIXct(idx,origin="1970-01-01")
    sw$dates <- as.POSIXct(idx,origin="1970-01-01")
    rh$dates <- as.POSIXct(idx,origin="1970-01-01")
    airT$dates <- as.POSIXct(idx,origin="1970-01-01")
    lwnet$dates <- as.POSIXct(idx,origin="1970-01-01")
    wnd$dates <- as.POSIXct(idx,origin="1970-01-01")
  }

  # if data isn't downsampled, define times and variable length
  if (!OC$TT$dwnSmple){
    if (OC$TT$openWtr){
      dates <- wt$dates
      varL <- length(dates)
    }
    if (OC$TT$openAirT){
      dates <- airT$dates
      varL <- length(dates)
    }
    if (OC$TT$openRH){
      dates <- rh$dates
      varL <- length(dates)
    }
    if (OC$TT$openSW){
      dates <- sw$dates
      varL <- length(dates)
    }
    if (OC$TT$openWnd){
      dates <- wnd$dates
      varL <- length(dates)
    }
    if (OC$TT$openLWnet){
      dates <- lwnet$dates
      varL <- length(dates)
    }
  }

  #****-----varL is the length of output files as of here-------*****
  # water temperature
  if (OC$TT$wrt_wTemp){
    OC$OC$writeTable$wTemp <- wt$dat
  }

  # calculate surface fluxes
  if (OC$TT$senslatYes || OC$TT$QtotYes){
    mm <- sens_latent(wt$dat,wnd$dat,airT$dat,rh$dat,Cfg$wndH,Cfg$htH,Cfg$hqH,Cfg$alt,Cfg$lat)
  }

  # atmospheric stability
  if (OC$TT$wrt_obu){
    zL1 <- Cfg$wndH/mm[,35]
    zL1[zL1 > 15] <- 15
    zL1[zL1 < -15] <- -15
    OC$OC$writeTable$obu <- zL1
  }

  # momentum flux
  if (OC$TT$wrt_tau){
    OC$OC$writeTable$tau <- mm[,1]
  }

  # sensible heat flux
  if (OC$TT$wrt_Qh || OC$TT$QtotYes){
    Qh <- mm[,3]
    if (OC$TT$wrt_Qh){
      OC$OC$writeTable$Qh <- Qh
    }
  }

  # latent heat flux
  if (OC$TT$wrt_Qe || OC$TT$QtotYes){
    Qe <- mm[,2]
    if (OC$TT$wrt_Qe){
      OC$OC$writeTable$Qe <- Qe
    }
  }

  # air shear velocity
  if (OC$TT$wrt_uSt_a){
    OC$OC$writeTable$uSt_a <- mm[,4]
  }

  # air shear velocity (neutral)
  if (OC$TT$wrt_uSt_aN){
    mm2 <- neutral_transfer_coeff(Uz = wnd$dat,hu = Cfg$wndH)
    OC$OC$writeTable$uSt_aN <- mm2[,1]
  }

  # wind speed at 10 m
  if (OC$TT$wrt_u10){
    OC$OC$writeTable$u10 <- mm[,7]
  }

  # wind speed at 10 m (neutral)
  if (OC$TT$wrt_u10N){
    mm2 <- neutral_transfer_coeff(wnd$dat,Cfg$wndH)
    OC$OC$writeTable$u10N <- mm2[,2]
  }

  # air temperature at 10 m
  if (OC$TT$wrt_t10){
    OC$OC$writeTable$t10 <- mm[,8]
  }

  # relative humidity at 10 m
  if (OC$TT$wrt_rh10){
    OC$OC$writeTable$rh10 <- mm[,10]
  }

  # transfer coefficient for momentum
  if (OC$TT$wrt_C_D){
    OC$OC$writeTable$C_D <- mm[,14]
  }

  # transfer coefficient for heat
  if (OC$TT$wrt_C_E){
    OC$OC$writeTable$C_E <- mm[,15]
  }

  # transfer coefficient for humidity
  if (OC$TT$wrt_C_H){
    OC$OC$writeTable$C_H <- mm[,16]
  }

  # transfer coefficient for momentum at 10 m
  if (OC$TT$wrt_C_D10){
    OC$OC$writeTable$C_D10 <- mm[,17]
  }

  # transfer coefficient for heat at 10 m
  if (OC$TT$wrt_C_E10){
    OC$OC$writeTable$C_E10 <- mm[,18]
  }

  # transfer coefficient for humidity at 10 m
  if (OC$TT$wrt_C_H10){
    OC$OC$writeTable$C_H10 <- mm[,19]
  }

  # neutral transfer coefficient for momentum
  if (OC$TT$wrt_C_D10N){
    mm2 <- neutral_transfer_coeff(wnd$dat,Cfg$wndH)
    OC$OC$writeTable$C_D10N <- mm2[,3]
  }

  # neutral drag coefficient for heat
  if (OC$TT$wrt_C_E10N){
    mm2 <- neutral_transfer_coeff(wnd$dat,Cfg$wndH)
    OC$OC$writeTable$C_E10N <- mm2[,4]
  }

  # neutral drag coefficient for humidity
  if (OC$TT$wrt_C_H10N){
    mm2 <- neutral_transfer_coeff(wnd$dat,Cfg$wndH)
    OC$OC$writeTable$C_H10N <- mm2[,5]
  }

  # neutral transfer coefficient for momentum
  if (OC$TT$wrt_C_DN){
    mm2 <- neutral_transfer_coeff(wnd$dat,Cfg$wndH)
    OC$OC$writeTable$C_DN <- mm2[,6]
  }

  # neutral drag coefficient for heat
  if (OC$TT$wrt_C_EN){
    mm2 <- neutral_transfer_coeff(wnd$dat,Cfg$wndH)
    OC$OC$writeTable$C_EN <- mm2[,7]
  }

  # neutral drag coefficient for humidity
  if (OC$TT$wrt_C_HN){
    mm2 <- neutral_transfer_coeff(wnd$dat,Cfg$wndH)
    OC$OC$writeTable$C_HN <- mm2[,8]
  }

  # evaporation
  if (OC$TT$wrt_Evap){
    OC$OC$writeTable$Evap <- mm[,27]
  }

  # net long wave heat flux
  if (OC$TT$wrt_Qlnet){
    OC$OC$writeTable$Qlnet <- lwnet$dat
  }

  # incoming long wave heat flux
  if (OC$TT$wrt_Qlin){
    press <- 101325*(1 - 2.25577E-5*Cfg$alt)^5.25588 # Pa
    press <- press/100 # mb
    lwm <- calc_lwnet(dates,Cfg$lat,press,airT$dat,rh$dat,sw$dat,wt$dat)$lw
    OC$OC$writeTable$Qlin <- lwm
  }

  # outgoing long wave heat flux
  if (OC$TT$wrt_Qlout){
    Tk <- wt$dat + 273.13
    emiss <- 0.972
    S_B <- 5.67E-8
    LWo <- S_B*emiss*Tk^4
    OC$OC$writeTable$Qlout <- LWo
  }

  # reflected short wave radiaiton
  if (OC$TT$wrt_Qsr || OC$TT$QtotYes){
    if (Cfg$outRs >= 86400){
      #datev <- datevec(dates)
      #datev[,4] <- 12
      dates2 <- as.POSIXct(paste0(format(dates,"%Y-%m-%d")," 12:00"))
      sw_alb <- sw_albedo(dates2,Cfg$lat)
    } else {
      sw_alb <- sw_albedo(dates,Cfg$lat)
    }
    Qsr <- sw$dat*sw_alb # reflected short wave radiation
    if (OC$TT$wrt_Qsr){
      OC$OC$writeTable$Qsr <- Qsr
    }
  }

  # total surface heat flux
  if (OC$TT$wrt_Qtot){
    # Qtot <- sw - Qsr - Qe - Qh + lwnet
    Qtot <- sw$dat - Qsr - Qe - Qh + lwm - LWo
    OC$OC$writeTable$Qtot <- Qtot
  }

  # short wave radiation
  if (OC$TT$wrt_Qs){
    OC$OC$writeTable$Qs <- sw$dat
  }

  # short wave radiation
  if (OC$TT$wrt_Qsin){
    if (Cfg$outRs >= 86400) {# quick fix for >daily averages
      ## needs translation ##
      #datev <- datevec(dates)
      #datev(:,4) <- 12
      #dates2 <- datenum(datev)
      sw_alb <- sw_albedo(dates2,Cfg$lat)
    } else {
      sw_alb <- sw_albedo(dates,Cfg$lat)
    }
    Qsr <- sw$dat*sw_alb # reflected short wave radiation
    OC$OC$writeTable$Qsin <- sw$dat - Qsr
  }

  # air density 10 m
  if (OC$TT$wrt_rhoa10){
    airT10 <- mm[,8]
    rh10 <- mm[,10]
    press <- 101325*(1 - 2.25577e-5*Cfg$alt)^5.25588 # Pa
    press <- press/100 # mb
    e_s <- 6.11*exp(17.27*airT10/(237.3 + airT10)) # saturated vapour pressure at ta, mb
    e_a <- rh10*e_s/100 # vapour pressure, mb
    q_z <- 0.622*e_a/press # specific humidity, kg kg-1
    R_a <- 287*(1 + 0.608*q_z)
    rhoa10 <- 100*press/(R_a*(airT10 + 273.16))
    OC$OC$writeTable$rhoa10 <- rhoa10
  }

  # water density
  if (OC$TT$wrt_rhow){
    rhow <- 1000*(1-1.9549*0.00001*abs(wt$dat-3.84)^1.68)
    OC$OC$writeTable$rhow <- rhow
  }

  # air density
  if (OC$TT$wrt_rhoa){
    press <- 101325*(1 - 2.25577e-5*Cfg$alt)^5.25588 # Pa
    press <- press/100 # mb
    e_s <- 6.11*exp(17.27*airT$dat/(237.3 + airT$dat)) # saturated vapour pressure at ta, mb
    e_a <- rh$dat*e_s/100 # vapour pressure, mb
    q_z <- 0.622*e_a/press # specific humidity, kg kg-1
    R_a <- 287*(1 + 0.608*q_z)
    rhoa <- 100*press/(R_a*(airT$dat + 273.16))
    OC$OC$writeTable$rhoa <- rhoa
  }

  # build plot array
  if (Cfg$plotYes){
    cat('Plotting results \n')
    plotLA_results(OC$writeTable,OC$plotTable,dates,LakeName,Folder)
    cat('completed\n\n')
  }

  # build file array
  writeNames <- list()
  cnt <- 1
  for (k in 1:length(OC$outputOptions)){
    if (length(OC$OC$writeTable[[OC$outputOptions[k]]])>0){
      writeNames[[cnt]] <- OC$outputOptions[k]
      cnt <- cnt+1
    }
  }

  # write to file
  if (Cfg$writeYes){
    cat('Writing results to file \n')
  }
  ##herehere##
  if (Cfg$writeYes && length(writeNames)>0){
    outputFile <- paste0(Folder, '/', LakeName, '_results.txt')
    outFile <- file(outputFile)
    if (isOpen(outFile)){
      stop(paste0(Folder, '/', LakeName, '_results.csv file in use, please close'))
    }
    close(outFile)
    #wrt <- function(writer,ap=TRUE){cat(writer,file = outputFile,append = ap)} # build a subfunction that writes
    # the contents of the input "writer"
    # to the file everytime wrt is called
    #wrt('DateTime',ap=FALSE)
    coln <- 'DateTime'
    for (i in 1:(cnt-1)){
      #wrt(paste0(OC$delimO," ", writeNames[[i]], OC$plotTable[[writeNames[[i]]]]["YLabel"],
      #          collapse = ""))
      coln <- c(coln,paste0(writeNames[[i]], OC$plotTable[[writeNames[[i]]]]["YLabel"],
                 collapse = ""))
    }
    #wrt('\r\n')
    #for (j in 1:varL){
      #wrt(format(dates[j],OC$dateOutput)) #change 'dateOutput'
      # in the 'OutputConstructor.m' file
      data_out <- data.frame(format(dates,OC$dateOutput))
      for (i in 1:length(writeNames)){
        #wrt(paste0(OC$delimO," ", as.numeric(OC$OC$writeTable[[writeNames[[i]]]][j]),
        #            collapse=""))
        data_out <- cbind(data_out,as.numeric(OC$OC$writeTable[[writeNames[[i]]]]))
      }
      #wrt('\r\n')
    #}
  colnames(data_out) <- coln
  data_out[,2:length(writeNames)] <- round(data_out[,2:length(writeNames)],7)
  write.table(data_out,outputFile,sep = OC$delimO,row.names = FALSE,quote = FALSE)
  }
  if (Cfg$writeYes){
    cat('completed\n\n')
  }
  cat('Lake Heat Flux Analyzer is complete')
  #profile off
  #profile viewer
}
