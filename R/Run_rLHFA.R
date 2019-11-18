#' run rLakeHeatFluxAnalyzer
#'
#' @param LakeName Name of the Lake. All forcing and config files must have this name
#' @param Folder Folder the forcing and config files are stored in
#' @param skipLoad Use own config file (must then be placed in Folder)
#' @export

Run_LHFA <- function(LakeName,Folder,skipLoad=FALSE)
#----Author: Jordan S Read 2009 ----
  #----Modified by R. Iestyn Woolway ----
  #---- R translation by Johannes Feldbauer -----

  lhfa_version <- '1.1.2'


  rm(list = ls())
  graphics.off()
  cat("\f")

  done <- FALSE
  if (!skipLoad){
    build_config(LakeName,Folder)
  }

  cat(paste0('Running LHFA version ', lhfa_version, '\n\n'))

  # -- variables --
  matSec <- 86400         # number of seconds in a day
  smplTs <- 50            # samplerate test length
  dateTl <- 0.00001       # tolerance for equal day (day fraction)
  cat(pasteo('Reading ', LakeName, '.hfx file'))
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
    cat(['Checking for ' LakeName '.plt file*'])
    pltFileName  <- paste0(Folder, '/', LakeName, '.plt')
    oper <- file.exists(pltFileName)
    if (!oper){
      cat('not found [*is optional]\n\n')
    } else {

      pltMods <- read.table(pltFileName)
      cat('completed\n')
    }
  }
  OC <- OutputConstructor(outPuts,pltMods)
  fNms <- names(OC$TT)
  for (j in 1:length(fNms)){
    if (OC$TT[[fNms[j]]]){
      cat(paste0(fNms[j], '\n'))
    }
  }
  cat('****completed****\n\n')

  if (OC$TT$openWtr){
    cat(paste0('Reading ', LakeName, '.wtr file'))
    wtrFileName  <- paste0(Folder, '/', LakeName, '.wtr')
    oper <- file.exists(wtrFileName)
    if (!oper){
      stop(paste0(LakeName, '.wtr file not found'))
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
                 ' which will be used to represent surface water temperatures'))
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
##+++++++++++++++++++++ here I am +++++++++++++++++++++++++++##
  # *** down sampling ***
    if (OC$TT$dwnSmple){
      cat(['Down sampling ' LakeName '.wtr data'])
      [DS_wtrD,DS_wtr] <- DownSample_TS(wt$dates,Cfg$outRs,wt$dat)
      wt$dates <- DS_wtrD
      wt$dat <- DS_wtr
      varL <- length(wt$dates)
      clear DS_wtrD DS_wtr
      cat('completed\n\n')
    }
  }

  if OC$TT$openWnd
  cat(['Reading ' LakeName '.wnd file'])
  wndFileName  <- [Folder '/' LakeName '.wnd']
  oper <- fopen(wndFileName)
  if eq(oper,-1)
  error([LakeName '.wnd file not found'])
  end
  fclose all
  [wndD,wnd] <- gFileOpen(wndFileName)
  wnd(wnd < Cfg$wndMn) <- Cfg$wndMn
  wnd(wnd > Cfg$wndMx) <- Cfg$wndMx

  # remove nans
  idx <- isnan(wnd)
  wnd(idx) <- []
  wndD(idx) <- []
  cat('completed\n\n')

  #*** find samplerate of raw data
  if length(wndD) < smplTs
  tLen <- length(wndD)
  else
    tLen <- smplTs
  end
  steps <- NaN(1,tLen-1)
  for i <- 1:tLen-1
  steps(i) <- wndD(i+1)-wndD(i)
  end
  if eq(min(steps),0)
  matRs <- mean(steps)*matSec
  else
    matRs <- min(steps)*matSec #current sample rate of raw data in seconds
  end
  clear vals ind numMx numCont steps tLen
  if (Cfg$outRs - matRs)>dateTl
  OC$TT$dwnSmple <- TRUE   #down sample if necessary
  end

  # *** down sampling ***
    if OC$TT$dwnSmple
  cat(['Down sampling ' LakeName '.wnd data'])
  [DS_wndD,DS_wnd] <- DownSample_TS(wndD,Cfg$outRs,wnd)
  wndD <- DS_wndD
  wnd <- DS_wnd
  varL <- length(wndD)
  clear DS_wndD DS_wnd
  cat('completed\n\n')
  end
  dates <- wndD
  end

  if OC$TT$openSW
  if exist([Folder '/' LakeName '.sw']) > 0
  cat(['Reading ' LakeName '.sw file'])
  swFileName  <- [Folder '/' LakeName '.sw']
  fclose all
  [swD,sw] <- gFileOpen(swFileName)
  sw(sw < 0) <- 0

  # remove nans
  idx <- isnan(sw)
  sw(idx) <- []
  swD(idx) <- []
  cat('completed\n\n')

  #*** find samplerate of raw data
  if length(swD) < smplTs
  tLen <- length(swD)
  else
    tLen <- smplTs
  end
  steps <- NaN(1,tLen-1)
  for i <- 1:tLen-1
  steps(i) <- swD(i+1)-swD(i)
  end
  if eq(min(steps),0)
  matRs <- mean(steps)*matSec
  else
    matRs <- min(steps)*matSec #current sample rate of raw data in seconds
  end
  clear vals ind numMx numCont steps tLen
  if (Cfg$outRs - matRs)>dateTl
  OC$TT$dwnSmple <- TRUE   #down sample if necessary
  end

  # *** down sampling ***
    if OC$TT$dwnSmple
  cat(['Down sampling ' LakeName '.sw data'])
  [DS_swD,DS_sw] <- DownSample_TS(swD,Cfg$outRs,sw)
  swD <- DS_swD
  sw <- DS_sw
  varL <- length(swD)
  clear DS_swD DS_sw
  cat('completed\n\n')
  end

  else if exist([Folder '/' LakeName '.par']) > 0
  cat(['Reading ' LakeName '.par file'])
  swFileName  <- [Folder '/' LakeName '.par']
  fclose all
  [swD,sw] <- gFileOpen(swFileName)
  sw(sw < 0) <- 0
  # remove nans
  idx <- isnan(sw)
  sw(idx) <- []
  swD(idx) <- []
  parMult <- 0.4957
  sw <- sw.*parMult
  cat('completed\n\n')

  #*** find samplerate of raw data
  if length(swD) < smplTs
  tLen <- length(swD)
  else
    tLen <- smplTs
  end
  steps <- NaN(1,tLen-1)
  for i <- 1:tLen-1
  steps(i) <- swD(i+1)-swD(i)
  end
  if eq(min(steps),0)
  matRs <- mean(steps)*matSec
  else
    matRs <- min(steps)*matSec #current sample rate of raw data in seconds
  end
  clear vals ind numMx numCont steps tLen
  if (Cfg$outRs - matRs)>dateTl
  OC$TT$dwnSmple <- TRUE   #down sample if necessary
  end

  # *** down sampling ***
    if OC$TT$dwnSmple
  cat(['Down sampling ' LakeName '.sw data'])
  [DS_swD,DS_sw] <- DownSample_TS(swD,Cfg$outRs,sw)
  swD <- DS_swD
  sw <- DS_sw
  varL <- length(swD)
  clear DS_swD DS_sw
  cat('completed\n\n')
  end

  else
    error([LakeName '.sw nor .par file not found'])
  end
  end
  dates <- swD
  end

  if OC$TT$openAirT
  cat(['Reading ' LakeName '.airT file'])
  airTFileName  <- [Folder '/' LakeName '.airT']
  oper <- fopen(airTFileName)
  if eq(oper,-1)
  error([LakeName '.airT file not found'])
  end
  fclose all
  [airTD,airT] <- gFileOpen(airTFileName)

  # remove nans
  idx <- isnan(airT)
  airT(idx) <- []
  airTD(idx) <- []
  cat('completed\n\n')

  #*** find samplerate of raw data
  if length(airTD) < smplTs
  tLen <- length(airTD)
  else
    tLen <- smplTs
  end
  steps <- NaN(1,tLen-1)
  for i <- 1:tLen-1
  steps(i) <- airTD(i+1)-airTD(i)
  end
  if eq(min(steps),0)
  matRs <- mean(steps)*matSec
  else
    matRs <- min(steps)*matSec #current sample rate of raw data in seconds
  end
  clear vals ind numMx numCont steps tLen
  if (Cfg$outRs - matRs)>dateTl
  OC$TT$dwnSmple <- TRUE   #down sample if necessary
  end

  # *** down sampling ***
    if OC$TT$dwnSmple
  cat(['Down sampling ' LakeName '.airT data'])
  [DS_airTD,DS_airT] <- DownSample_TS(airTD,Cfg$outRs,airT)
  airTD <- DS_airTD
  airT <- DS_airT
  varL <- length(airTD)
  clear DS_airTD DS_airT
  cat('completed\n\n')
  end
  end

  if OC$TT$openRH
  cat(['Reading ' LakeName '.rh file'])
  rhFileName  <- [Folder '/' LakeName '.rh']
  oper <- fopen(rhFileName)
  if eq(oper,-1)
  error([LakeName '.rh file not found'])
  end
  fclose all
  [rhD,rh] <- gFileOpen(rhFileName)
  rh(rh < 0) <- 0
  rh(rh > 100) <- 100

  # remove nans
  idx <- isnan(rh)
  rh(idx) <- []
  rhD(idx) <- []
  cat('completed\n\n')

  #*** find samplerate of raw data
  if length(rhD) < smplTs
  tLen <- length(rhD)
  else
    tLen <- smplTs
  end
  steps <- NaN(1,tLen-1)
  for i <- 1:tLen-1
  steps(i) <- rhD(i+1)-rhD(i)
  end
  if eq(min(steps),0)
  matRs <- mean(steps)*matSec
  else
    matRs <- min(steps)*matSec #current sample rate of raw data in seconds
  end
  clear vals ind numMx numCont steps tLen
  if (Cfg$outRs - matRs)>dateTl
  OC$TT$dwnSmple <- TRUE   #down sample if necessary
  end

  # *** down sampling ***
    if OC$TT$dwnSmple
  cat(['Down sampling ' LakeName '.rh data'])
  [DS_rhD,DS_rh] <- DownSample_TS(rhD,Cfg$outRs,rh)
  rhD <- DS_rhD
  rh <- DS_rh
  varL <- length(rhD)
  clear DS_rhD DS_rh
  cat('completed\n\n')
  end
  end

  # make sure dates match
  if !OC$TT$openSW
  if OC$TT$openWtr && OC$TT$openWnd && OC$TT$openRH && OC$TT$openAirT
  idx <- intersect(intersect(intersect(wt$dates,wndD),rhD),airTD)
  dates <- idx
  varL <- length(dates)
  wt$dat <- wt$dat(ismember(wt$dates,idx))
  wnd <- wnd(ismember(wndD,idx))
  rh <- rh(ismember(rhD,idx))
  airT <- airT(ismember(airTD,idx))
  wt$dates <- idx
  wndD <- idx
  rhD <- idx
  airTD <- idx
  end
  end

  if !OC$TT$openWnd
  if OC$TT$openRH && OC$TT$openAirT && OC$TT$openWtr && OC$TT$openSW
  idx <- intersect(intersect(intersect(rhD,airTD),wt$dates),swD)
  dates <- idx
  varL <- length(dates)
  rh <- rh(ismember(rhD,idx))
  airT <- airT(ismember(airTD,idx))
  wt$dat <- wt$dat(ismember(wt$dates,idx))
  sw <- sw(ismember(swD,idx))
  wt$dates <- idx
  swD <- idx
  rhD <- idx
  airTD <- idx
  end
  end

  if OC$TT$openRH && OC$TT$openAirT && OC$TT$openWtr && OC$TT$openSW && OC$TT$openWnd
  idx <- intersect(intersect(intersect(intersect(rhD,airTD),wt$dates),swD),wndD)
  dates <- idx
  varL <- length(dates)
  rh <- rh(ismember(rhD,idx))
  airT <- airT(ismember(airTD,idx))
  wt$dat <- wt$dat(ismember(wt$dates,idx))
  sw <- sw(ismember(swD,idx))
  wnd <- wnd(ismember(wndD,idx))
  wt$dates <- idx
  swD <- idx
  rhD <- idx
  airTD <- idx
  wndD <- idx
  end

  if !OC$TT$openRH && !OC$TT$openAirT && !OC$TT$openSW && !OC$TT$openWnd
  if OC$TT$openWtr
  dates <- wt$dates
  varL <- length(dates)
  end
  end

  # look for long-wave radiation data
  if OC$TT$openLWnet
  if exist([Folder '/' LakeName '.lwnet']) > 0
  cat(['Reading ' LakeName '.lwnet file'])
  lwnetFileName <- [Folder '/' LakeName '.lwnet']
  fclose all
  [lwnetD,lwnet] <- gFileOpen(lwnetFileName)
  idx <- isnan(lwnet)
  lwnet(idx) <- []
  lwnetD(idx) <- []
  cat('completed\n\n')

  #*** find samplerate of raw data
  if length(lwnetD) < smplTs
  tLen <- length(lwnetD)
  else
    tLen <- smplTs
  end
  steps <- NaN(1,tLen-1)
  for i <- 1:tLen-1
  steps(i) <- lwnetD(i+1)-lwnetD(i)
  end
  if eq(min(steps),0)
  matRs <- mean(steps)*matSec
  else
    matRs <- min(steps)*matSec #current sample rate of raw data in seconds
  end
  clear vals ind numMx numCont steps tLen
  if (Cfg$outRs - matRs)>dateTl
  OC$TT$dwnSmple <- TRUE   #down sample if necessary
  end

  # *** down sampling ***
    if OC$TT$dwnSmple
  cat(['Down sampling ' LakeName '.lwnet data'])
  [DS_lwnetD,DS_lwnet] <- DownSample_TS(lwnetD,Cfg$outRs,lwnet)
  lwnetD <- DS_lwnetD
  lwnet <- DS_lwnet
  varL <- length(lwnetD)
  clear DS_lwnetD DS_lwnet
  cat('completed\n\n')
  end

  else if exist([Folder '/' LakeName '.lw']) > 0
  cat([LakeName '.lwnet files not found, looking for .lw file\n\n'])
  cat(['Reading ' LakeName '.lw file'])
  lwFileName  <- [Folder '/' LakeName '.lw']
  fclose all
  [lwD,lw] <- gFileOpen(lwFileName)
  idx <- isnan(lw)
  lw(idx) <- []
  lwD(idx) <- []
  cat('completed\n\n')

  #*** find samplerate of raw data
  if length(lwD) < smplTs
  tLen <- length(lwD)
  else
    tLen <- smplTs
  end
  steps <- NaN(1,tLen-1)
  for i <- 1:tLen-1
  steps(i) <- lwD(i+1)-lwD(i)
  end
  if eq(min(steps),0)
  matRs <- mean(steps)*matSec
  else
    matRs <- min(steps)*matSec #current sample rate of raw data in seconds
  end
  clear vals ind numMx numCont steps tLen
  if (Cfg$outRs - matRs)>dateTl
  OC$TT$dwnSmple <- TRUE   #down sample if necessary
  end

  # *** down sampling ***
    if OC$TT$dwnSmple
  cat(['Down sampling ' LakeName '.lw data'])
  [DS_lwD,DS_lw] <- DownSample_TS(lwD,Cfg$outRs,lw)
  lwD <- DS_lwD
  lw <- DS_lw
  varL <- length(lwD)
  clear DS_lwD DS_lw
  cat('completed\n\n')
  end

  # find when wtr and lw dates intersect
  idx <- intersect(wt$dates,lwD)
  lw <- lw(ismember(lwD,idx))
  wt$dat <- wt$dat(ismember(wt$dates,idx))
  wt$dates <- idx

  Tk <- wt$dat + 273.13 # .wtr already called at this point
  emiss <- 0.972
  S_B <- 5.67E-8
  LWo <- S_B*emiss*Tk.^4

  # define lwnet
  lwnet <- lw - LWo
  lwnet <- lwnet
  lwnetD <- idx
  else
    cat(['' LakeName '.lwnet and .lw file missing, using .airT, .rh, .sw, .wtr instead\n\n'])

  press <- 101325.*(1 - 2.25577e-5.*Cfg$alt).^5.25588 # Pa
  press <- press./100 # mb
  [!,!,Qlnet] <- calc_lwnet(dates,Cfg$lat,press,airT,rh,sw,wt$dat)
  lwnet <- Qlnet
  lwnetD <- dates
  end
  end
  end

  # re-adjust dates depending on lw and lwnet files
  if OC$TT$openRH && OC$TT$openAirT && OC$TT$openWtr && OC$TT$openSW && OC$TT$openLWnet
  idx <- intersect(intersect(intersect(intersect(rhD,airTD),wt$dates),swD),lwnetD)
  dates <- idx
  varL <- length(dates)
  rh <- rh(ismember(rhD,idx))
  airT <- airT(ismember(airTD,idx))
  wt$dat <- wt$dat(ismember(wt$dates,idx))
  sw <- sw(ismember(swD,idx))
  lwnet <- lwnet(ismember(lwnetD,idx))
  wt$dates <- idx
  swD <- idx
  rhD <- idx
  airTD <- idx
  lwnetD <- idx
  end

  if OC$TT$openRH && OC$TT$openAirT && OC$TT$openWtr && OC$TT$openSW && OC$TT$openLWnet && OC$TT$openWnd
  idx <- intersect(intersect(intersect(intersect(intersect(rhD,airTD),wt$dates),swD),lwnetD),wndD)
  dates <- idx
  varL <- length(dates)
  rh <- rh(ismember(rhD,idx))
  airT <- airT(ismember(airTD,idx))
  wt$dat <- wt$dat(ismember(wt$dates,idx))
  sw <- sw(ismember(swD,idx))
  lwnet <- lwnet(ismember(lwnetD,idx))
  wnd <- wnd(ismember(wndD,idx))
  wt$dates <- idx
  swD <- idx
  rhD <- idx
  airTD <- idx
  lwnetD <- idx
  wndD <- idx
  end

  # if data isn't downsampled, define times and variable length
  if !OC$TT$dwnSmple
  if OC$TT$openWtr
  dates <- wt$dates
  varL <- length(dates)
  end
  if OC$TT$openAirT
  dates <- airTD
  varL <- length(dates)
  end
  if OC$TT$openRH
  dates <- rhD
  varL <- length(dates)
  end
  if OC$TT$openSW
  dates <- swD
  varL <- length(dates)
  end
  if OC$TT$openWnd
  dates <- wndD
  varL <- length(dates)
  end
  if OC$TT$openLWnet
  dates <- lwnetD
  varL <- length(dates)
  end
  end

  #****-----varL is the length of output files as of here-------*****
  # water temperature
  if OC$TT$wrt_wTemp
  OC$OC$writeTable$wTemp <- wt$dat
  end

  # calculate surface fluxes
  if OC$TT$senslatYes || OC$TT$QtotYes
  mm <- sens_latent(wt$dat,wnd,airT,rh,Cfg$wndH,Cfg$htH,Cfg$hqH,Cfg$alt,Cfg$lat)
  end

  # atmospheric stability
  if OC$TT$wrt_obu
  zL1 <- Cfg$wndH./mm(:,35)
  zL1(zL1 > 15) <- 15
  zL1(zL1 < -15) <- -15
  OC$OC$writeTable$obu <- zL1
  end

  # momentum flux
  if OC$TT$wrt_tau
  OC$OC$writeTable$tau <- mm(:,1)
  end

  # sensible heat flux
  if OC$TT$wrt_Qh || OC$TT$QtotYes
  Qh <- mm(:,3)
  if OC$TT$wrt_Qh
  OC$OC$writeTable$Qh <- Qh
  end
  end

  # latent heat flux
  if OC$TT$wrt_Qe || OC$TT$QtotYes
  Qe <- mm(:,2)
  if OC$TT$wrt_Qe
  OC$OC$writeTable$Qe <- Qe
  end
  end

  # air shear velocity
  if OC$TT$wrt_uSt_a
  OC$OC$writeTable$uSt_a <- mm(:,4)
  end

  # air shear velocity (neutral)
  if OC$TT$wrt_uSt_aN
  mm2 <- neutral_transfer_coeff(wnd,Cfg$wndH)
  OC$OC$writeTable$uSt_aN <- mm2(:,1)
  end

  # wind speed at 10 m
  if OC$TT$wrt_u10
  OC$OC$writeTable$u10 <- mm(:,7)
  end

  # wind speed at 10 m (neutral)
  if OC$TT$wrt_u10N
  mm2 <- neutral_transfer_coeff(wnd,Cfg$wndH)
  OC$OC$writeTable$u10N <- mm2(:,2)
  end

  # air temperature at 10 m
  if OC$TT$wrt_t10
  OC$OC$writeTable$t10 <- mm(:,8)
  end

  # relative humidity at 10 m
  if OC$TT$wrt_rh10
  OC$OC$writeTable$rh10 <- mm(:,10)
  end

  # transfer coefficient for momentum
  if OC$TT$wrt_C_D
  OC$OC$writeTable$C_D <- mm(:,14)
  end

  # transfer coefficient for heat
  if OC$TT$wrt_C_E
  OC$OC$writeTable$C_E <- mm(:,15)
  end

  # transfer coefficient for humidity
  if OC$TT$wrt_C_H
  OC$OC$writeTable$C_H <- mm(:,16)
  end

  # transfer coefficient for momentum at 10 m
  if OC$TT$wrt_C_D10
  OC$OC$writeTable$C_D10 <- mm(:,17)
  end

  # transfer coefficient for heat at 10 m
  if OC$TT$wrt_C_E10
  OC$OC$writeTable$C_E10 <- mm(:,18)
  end

  # transfer coefficient for humidity at 10 m
  if OC$TT$wrt_C_H10
  OC$OC$writeTable$C_H10 <- mm(:,19)
  end

  # neutral transfer coefficient for momentum
  if OC$TT$wrt_C_D10N
  mm2 <- neutral_transfer_coeff(wnd,Cfg$wndH)
  OC$OC$writeTable$C_D10N <- mm2(:,3)
  end

  # neutral drag coefficient for heat
  if OC$TT$wrt_C_E10N
  mm2 <- neutral_transfer_coeff(wnd,Cfg$wndH)
  OC$OC$writeTable$C_E10N <- mm2(:,4)
  end

  # neutral drag coefficient for humidity
  if OC$TT$wrt_C_H10N
  mm2 <- neutral_transfer_coeff(wnd,Cfg$wndH)
  OC$OC$writeTable$C_H10N <- mm2(:,5)
  end

  # neutral transfer coefficient for momentum
  if OC$TT$wrt_C_DN
  mm2 <- neutral_transfer_coeff(wnd,Cfg$wndH)
  OC$OC$writeTable$C_DN <- mm2(:,6)
  end

  # neutral drag coefficient for heat
  if OC$TT$wrt_C_EN
  mm2 <- neutral_transfer_coeff(wnd,Cfg$wndH)
  OC$OC$writeTable$C_EN <- mm2(:,7)
  end

  # neutral drag coefficient for humidity
  if OC$TT$wrt_C_HN
  mm2 <- neutral_transfer_coeff(wnd,Cfg$wndH)
  OC$OC$writeTable$C_HN <- mm2(:,8)
  end

  # evaporation
  if OC$TT$wrt_Evap
  OC$OC$writeTable$Evap <- mm(:,27)
  end

  # net long wave heat flux
  if OC$TT$wrt_Qlnet
  OC$OC$writeTable$Qlnet <- lwnet
  end

  # incoming long wave heat flux
  if OC$TT$wrt_Qlin
  press <- 101325.*(1 - 2.25577e-5.*Cfg$alt).^5.25588 # Pa
  press <- press./100 # mb
  [lw,!,!] <- calc_lwnet(dates,Cfg$lat,press,airT,rh,sw,wt$dat)
  OC$OC$writeTable$Qlin <- lw
  end

  # outgoing long wave heat flux
  if OC$TT$wrt_Qlout
  Tk <- wt$dat + 273.13
  emiss <- 0.972
  S_B <- 5.67E-8
  LWo <- S_B*emiss*Tk.^4
  OC$OC$writeTable$Qlout <- LWo
  end

  # reflected short wave radiaiton
  if OC$TT$wrt_Qsr || OC$TT$QtotYes
  if Cfg$outRs >= 86400
  datev <- datevec(dates)
  datev(:,4) <- 12
  dates2 <- datenum(datev)
  sw_alb <- sw_albedo(dates2,Cfg$lat)
  else
  sw_alb <- sw_albedo(dates,Cfg$lat)
  end
  Qsr <- sw.*sw_alb # reflected short wave radiation
  if OC$TT$wrt_Qsr
  OC$OC$writeTable$Qsr <- Qsr
  end
  end

  # total surface heat flux
  if OC$TT$wrt_Qtot
  # Qtot <- sw - Qsr - Qe - Qh + lwnet
  Qtot <- sw - Qsr - Qe - Qh + lw - LWo
  OC$OC$writeTable$Qtot <- Qtot
  end

  # short wave radiation
  if OC$TT$wrt_Qs
  OC$OC$writeTable$Qs <- sw
  end

  # short wave radiation
  if OC$TT$wrt_Qsin
  if Cfg$outRs >= 86400 # quick fix for >daily averages
  datev <- datevec(dates)
  datev(:,4) <- 12
  dates2 <- datenum(datev)
  sw_alb <- sw_albedo(dates2,Cfg$lat)
  else
  sw_alb <- sw_albedo(dates,Cfg$lat)
  end
  Qsr <- sw.*sw_alb # reflected short wave radiation
  OC$OC$writeTable$Qsin <- sw - Qsr
  end

  # air density 10 m
  if OC$TT$wrt_rhoa10
  airT10 <- mm(:,8)
  rh10 <- mm(:,10)
  press <- 101325.*(1 - 2.25577e-5.*Cfg$alt).^5.25588 # Pa
  press <- press./100 # mb
  e_s <- 6.11.*exp(17.27.*airT10./(237.3 + airT10)) # saturated vapour pressure at ta, mb
  e_a <- rh10.*e_s./100 # vapour pressure, mb
  q_z <- 0.622.*e_a./press # specific humidity, kg kg-1
  R_a <- 287.*(1 + 0.608.*q_z)
  rhoa10 <- 100*press./(R_a.*(airT10 + 273.16))
  OC$OC$writeTable$rhoa10 <- rhoa10
  end

  # water density
  if OC$TT$wrt_rhow
  rhow <- 1000*(1-1.9549*0.00001*abs(wt$dat-3.84).^1.68)
  OC$OC$writeTable$rhow <- rhow
  end

  # air density
  if OC$TT$wrt_rhoa
  press <- 101325.*(1 - 2.25577e-5.*Cfg$alt).^5.25588 # Pa
  press <- press./100 # mb
  e_s <- 6.11.*exp(17.27.*airT./(237.3 + airT)) # saturated vapour pressure at ta, mb
  e_a <- rh.*e_s./100 # vapour pressure, mb
  q_z <- 0.622.*e_a./press # specific humidity, kg kg-1
  R_a <- 287.*(1 + 0.608.*q_z)
  rhoa <- 100*press./(R_a.*(airT + 273.16))
  OC$OC$writeTable$rhoa <- rhoa
  end

  # build plot array
  if Cfg$plotYes
  cat('Plotting results')
  plotLA_results(OC$writeTable,OC$plotTable,dates,LakeName,Folder)
  cat('completed\n\n')
  end

  # build file array
  writeNames <- {}
  cnt <- 1
  for k <- 1:length(OC$outputOptions)
  if !islogical(OC$OC$writeTable$(char(OC$outputOptions{k})))
  writeNames{cnt} <- OC$outputOptions{k}
  cnt <- cnt+1
  end
  end

  # write to file
  if Cfg$writeYes
  cat('Writing results to file')
  end

  if Cfg$writeYes && !isempty(writeNames)
  outputFile <- [Folder '/' LakeName '_results.txt']
  outFile <- fopen(outputFile,'w')
  if eq(outFile,-1)
  error([Folder '/' LakeName '_results.csv file in use, please close'])
  end
  wrt <- @(writer)cat(outFile,writer) # build a subfunction that writes
  # the contents of the input "writer"
  # to the file everytime wrt is called
  wrt('DateTime')
  for i <- 1:cnt-1
  wrt([OC$delimO writeNames{i} OC$plotTable[[writeNames[i]]][YLabel][]])
  end
  wrt('\r\n')
  for j <- 1:varL
  wrt(datestr(dates(j),OC$dateOutput)) #change 'dateOutput'
  # in the 'OutputConstructor.m' file
  for i <- 1:length(writeNames)
  wrt([OC$delimO num2str(OC$OC$writeTable$(char(writeNames{i}))(j))])
  end
  wrt('\r\n')
  end
  fclose all
  end
  if Cfg$writeYes
  cat('completed\n\n')
  end
  disp('Lake Heat Flux Analyzer is complete')
  #profile off
  #profile viewer
  end
