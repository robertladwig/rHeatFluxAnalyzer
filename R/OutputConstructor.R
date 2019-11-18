#' Function to construct input data dependency
#'
#' @param outputNames Names of the variables that should be outputed
#' @param pltMods Modifications to the plot settings
#' @return Returns a list with several tables and lists for further calculations


OutputConstructor <- function(outputNames,pltMods){

  outputOptions <-  c('tau','Qh','Qe','C_DN','C_EN','C_HN','C_D10N','C_E10N','C_H10N','C_D',
                    'C_E','C_H','C_D10','C_E10','C_H10','u10','u10N','t10','rh10','Qlnet',
                    'Qlin','Qlout','uSt_a','uSt_aN','Evap','wTemp','Qsr','obu','Qtot','Qs',
                    'Qsin','rhoa10','rhow','rhoa')

  # writable outputs
  programOptions= c('openWtr','openWnd','openWnd','openRH','openAirT','openSW',
                    'openLWnet','dwnSmple','senslatYes','QtotYes')
  # program flow specifiers

  # Error Check
  if (length(outputNames)!=sum(outputNames%in%outputOptions)){
  error(paste0('output "',outputNames[!outputNames%in%outputOptions],
        '" not recognized'))
  }
  # Other defaults
  dateInput <-  'yyyy-mm-dd HH:MM'
  delimI <- '\t'

  dateOutput  <- 'yyyy-mm-dd HH:MM'
  delimO <- '\t'

  # Figure defaults
  isStringMod <- c('figUnits','figType','fontName','figRes')
  figUnits    <- 'inches'
  figWidth    <- 6   # relative to fig_units
  figHeight   <- 3   # relative to fig_units
  leftMargin  <- .75 # relative to fig_units
  rightMargin <- .1  # relative to fig_units
  topMargin   <- .4  # relative to fig_units
  botMargin   <- .4  # relative to fig_units
  figType     <- 'png'
  figRes      <- '150' # dots per inch (not relative to units?)
  fontName    <- 'Arial'
  fontSize    <- 12
  heatMapMin  <- 0
  heatMapMax  <- 30


  plt <- list(figUnits=figUnits,figWidth=figWidth,figHeight=figHeight,
              leftMargin=leftMargin,rightMargin=rightMargin,topMargin=topMargin,
              botMargin=botMargin,figType=figType,figRes=figRes,
              fontName=fontName,fontSize=fontSize,heatMapMin=heatMapMin,
              heatMapMax=heatMapMax)

  if (length(pltMods)>0){
    # use plot mods to modify plotting defaults
    fN <- names(pltMods)
    for (n in 1:length(fN)){
      if (!any(fN[n]==isStringMod)){
        if (pltMods[[fN[n]]]!=plt[[fN[n]]]){
          cat(paste0('>>User plot modification replacing ', fN[n], '=',
                     as.character(plt[[fN[n]]]), ' with ', pltMods[[fN[n]]],  '<<<\n'))
          plt[[fN[n]]] <- as.numeric(pltMods[[fN[n]]])
        }
      }
        if (pltMods[[fN[n]]]!=plt[[fN[n]]]){
            cat(paste0('>>User plot modification replacing ', fN[n], '=',
                       plt[[fN[n]]], ' with ', pltMods[[fN[n]]],  '<<<\n'))
            plt[[fN[n]]] <- as.character(pltMods[[fN[n]]])
          }
    }
  }


  fig_Defaults <- list(Units=plt$figUnits,
                       Color='w',
                       PaperUnits=plt$figUnits,
                       PaperPosition=c(0, 0, plt$figWidth, plt$figHeight),
                       Position=c(1, 1, plt$figWidth, plt$figHeight),
                       PaperPositionMode='manual',
                       PaperSize=c(plt$figWidth, plt$figHeight))

  print_Defaults <- list(format=c('-d', plt$figType),
                         res=c('-r', plt$figRes),
                         toClose=TRUE)

  position <- c(plt$leftMargin/plt$figWidth, plt$botMargin/plt$figHeight,
                (plt$figWidth-plt$leftMargin-plt$rightMargin)/plt$figWidth,
                (plt$figHeight-plt$topMargin-plt$botMargin)/plt$figHeight)

  axes_Defaults <- list(FontName=plt$fontName,
                        FontSize=plt$fontSize,
                        Layer='top',
                        Position=position,
                        Box='on',
                        YLabel='',
                        Title='',
                        YDir='normal',
                        YScale='linear',
                        CLim=c(0, 1))
  # build listures
  plotTable <- list(FigD=fig_Defaults,
                    PrintD=print_Defaults)
  truthTable <- list()
  writeTable <- list()
  outputConlist <- list()
  for (j in 1:length(programOptions)){
    truthTable[[programOptions[j]]] <- FALSE
  }
  for (j in 1:length(outputOptions)){
    truthTable[[paste0('wrt_', outputOptions[j])]] <- FALSE
    writeTable[[outputOptions[j]]] <- FALSE
    plotTable[[outputOptions[j]]] <- axes_Defaults
  }
  for (j in 1:length(outputOptions)){
    outputConlist[[outputOptions[j]]] <- truthTable
  }
  RunNeed <- list()
  RunAxes <- list()
  # Water Temperature
  name <- 'wTemp'
  RunNeed[[name]] <- c('openWtr','wrt_wTemp')
  RunAxes[[name]] <- list(YLabel=' (^{o} C)',
                          Title='Surface water temperature')

  # Momentum flux
  name <- 'tau'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','wrt_tau',
                       'senslatYes')
  RunAxes[[name]] <- list(YLabel=' (N m^{-2})',
                          Title='Surface flux of momentum')

  # Sensible heat flux
  name <- 'Qh'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','wrt_Qh',
                       'senslatYes')
  RunAxes[[name]] <- list(YLabel=' (W m^{-2})',
                          Title='Sensible heat flux')

  # Latent heat flux
  name <- 'Qe'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','wrt_Qe',
                       'senslatYes')
  RunAxes[[name]] <- list(YLabel=' (W m^{-2})',
                          Title='Latent heat flux')

  # Transfer coefficient for momentum (neutral)
  name <- 'C_DN'
  RunNeed[[name]] <- c('openWnd','wrt_C_DN')
  RunAxes[[name]] <- list(YLabel='',
                          Title='Transfer coefficient for momentum (neutral)')

  # Transfer coefficient for heat (neutral)
  name <- 'C_EN'
  RunNeed[[name]] <- c('openWnd','wrt_C_EN')
  RunAxes[[name]] <- list(YLabel='',
                          Title='Transfer coefficient for heat (neutral)')

  # Transfer coefficient for humidity (neutral)
  name <- 'C_HN'
  RunNeed[[name]] <- c('openWnd','wrt_C_HN')
  RunAxes[[name]] <- list(YLabel='',
                          Title='Transfer coefficient for humidity (neutral)')

  # Transfer coefficient for momentum at 10 m (neutral)
  name <- 'C_D10N'
  RunNeed[[name]] <- c('openWnd','wrt_C_D10N')
  RunAxes[[name]] <- list(YLabel='',
                          Title='Transfer coefficient for momentum at 10 m(neutral)')

  # Transfer coefficient for heat at 10 m (neutral)
  name <- 'C_E10N'
  RunNeed[[name]] <- c('openWnd','wrt_C_E10N')
  RunAxes[[name]] <- list(YLabel='',
                          Title='Transfer coefficient for heat at 10 m(neutral)')

  # Transfer coefficient for humidity at 10 m (neutral)
  name <- 'C_H10N'
  RunNeed[[name]] <- c('openWnd','wrt_C_H10N')
  RunAxes[[name]] <- list(YLabel='',
                          Title='Transfer coefficient for humidity at 10 m(neutral)')

  # Transfer coefficient for momentum
  name <- 'C_D'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','wrt_C_D',
                       'senslatYes')
  RunAxes[[name]] <- list(YLabel='',
                          Title='Transfer coefficient for momentum')

  # Transfer coefficient for heat
  name <- 'C_E'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','wrt_C_E',
                       'senslatYes')
  RunAxes[[name]] <- list(YLabel='',
                          Title='Transfer coefficient for heat')

  # Transfer coefficient for humidity
  name <- 'C_H'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','wrt_C_H',
                       'senslatYes')
  RunAxes[[name]] <- list(YLabel='',
                          Title='Transfer coefficient for humidity')

  # Transfer coefficient for momentum at 10 m
  name <- 'C_D10'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','wrt_C_D10',
                       'senslatYes')
  RunAxes[[name]] <- list(YLabel='',
                          Title='Transfer coefficient for momentum at 10 m')

  # Transfer coefficient for heat at 10 m
  name <- 'C_E10'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','wrt_C_E10',
                       'senslatYes')
  RunAxes[[name]] <- list(YLabel='',
                          Title='Transfer coefficient for heat at 10 m')

  # Transfer coefficient for humidity at 10 m
  name <- 'C_H10'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','wrt_C_H10',
                       'senslatYes')
  RunAxes[[name]] <- list(YLabel='',
                          Title='Transfer coefficient for humidity at 10 m')

  # Wind speed at 10 m
  name <- 'u10'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','wrt_u10',
                       'senslatYes')
  RunAxes[[name]] <- list(YLabel=' (m s^{-1})',
                          Title='Wind speed at 10 m')

  # Wind speed at 10 m (neutral)
  name <- 'u10N'
  RunNeed[[name]] <- c('openWnd','wrt_u10N')
  RunAxes[[name]] <- list(YLabel=' (m s^{-1})',
                          Title='Wind speed at 10 m (neutral)')

  # Relative humidity at 10 m
  name <- 'rh10'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','wrt_rh10',
                       'senslatYes')
  RunAxes[[name]] <- list(YLabel=as.character(37),
                          Title='Relative humidity at 10 m')

  # Air temperature at 10 m
  name <- 't10'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','wrt_t10',
                       'senslatYes')
  RunAxes[[name]] <- list(YLabel=' (^o{C})',
                          Title='Air temperature at 10 m')

  # Long wave heat flux
  name <- 'Qlnet'
  RunNeed[[name]] <- c('openRH','openAirT','openWtr','openSW',
                       'openLWnet','wrt_Qlnet')
  RunAxes[[name]] <- list(YLabel=' (W m^{-2})',
                          Title='Net long wave heat radiation')

  # Incoming long wave radiation
  name <- 'Qlin'
  RunNeed[[name]] <- c('openRH','openAirT','openWtr','openSW',
                       'open_LWnet','wrt_Qlin')
  RunAxes[[name]] <- list(YLabel=' (W m^{-2})',
                          Title='Incoming long wave radiation')

  # Outgoing long wave heat flux
  name <- 'Qlout'
  RunNeed[[name]] <- c('openWtr','wrt_Qlout')
  RunAxes[[name]] <- list(YLabel=' (W m^{-2})',
                          Title='Outgoing long wave radiation')

  # Air shear velocity
  name <- 'uSt_a'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','wrt_uSt_a',
                       'senslatYes')
  RunAxes[[name]] <- list(YLabel=' (m s^{-1})',
                          Title='Air shear velocity')

  # Air shear velocity neutral
  name <- 'uSt_aN'
  RunNeed[[name]] <- c('openWnd','wrt_uSt_aN')
  RunAxes[[name]] <- list(YLabel=' (m s^{-1})',
                          Title='Air shear velocity (neutral)')

  # Evaporation
  name <- 'Evap'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','wrt_Evap',
                       'senslatYes')
  RunAxes[[name]] <- list(YLabel=' (mm day^{-1})',
                          Title='Evaporation')

  # Reflected short wave radiation
  name <- 'Qsr'
  RunNeed[[name]] <- c('openSW','wrt_Qsr')
  RunAxes[[name]] <- list(YLabel=' (W m^{-2})',
                    Title='Reflected short wave radiation')

  # Monin-Obukhov length scale
  name <- 'obu'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','wrt_obu',
                       'senslatYes')
  RunAxes[[name]] <- list(YLabel=' zLw^{-1}',
                          Title='Atmospheric stability')

  # total heat flux
  name <- 'Qtot'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','openSW',
                       'openLWnet','wrt_Qtot','QtotYes')
  RunAxes[[name]] <- list(YLabel=' (W m^{-2})',
                          Title='Total surface heat flux')

  # short wave radiation
  name <- 'Qs'
  RunNeed[[name]] <- c('openSW','wrt_Qs')
  RunAxes[[name]] <- list(YLabel=' (W m^{-2})',
                          Title='Short wave radiation')

  # net incoming short wave radiation
  name <- 'Qsin'
  RunNeed[[name]] <- c('openSW','wrt_Qsin')
  RunAxes[[name]] <- list(YLabel=' (W m^{-2})',
                          Title='Net incoming short wave radiation')

  # air density at 10m
  name <- 'rhoa10'
  RunNeed[[name]] <- c('openWtr','openWnd','openRH','openAirT','wrt_rhoa10',
                       'senslatYes')
  RunAxes[[name]] <- list(YLabel=' (kg m^{-3})',
                          Title='Air density (10 m)')

  # water density
  name <- 'rhow'
  RunNeed[[name]] <- c('openWtr','wrt_rhow')
  RunAxes[[name]] <- list(YLabel=' (kg m^{-3})',
                          Title='Water density')

  # air density
  name <- 'rhoa'
  RunNeed[[name]] <- c('openAirT','openRH','wrt_rhoa')
  RunAxes[[name]] <- list(YLabel=' (kg m^{-3})',
                          Title='Air density')

  # Fill Conlist
  fNms <- names(outputConlist)
  for (k in 1:length(fNms)){
    strN <- RunNeed[[fNms[k]]]
    for (j in 1:length(strN)){
      outputConlist[[fNms[k]]][[strN[j]]] <- TRUE
    }
    axN <- names(RunAxes[[fNms[k]]])
    for (j in 1:length(axN)){
      plotTable[[fNms[k]]][[axN[j]]] <- RunAxes[[fNms[k]]][[axN[j]]]
    }
  }

  fNms <- names(truthTable)
  for (k in 1:length(fNms)){
    for (j in 1:length(outputNames)){
      if (outputConlist[[outputNames[j]]][[fNms[k]]]){
        truthTable[[fNms[k]]] <- TRUE
      }
    }
  }

  fNms <- names(writeTable)
  for (k in 1:length(fNms)){
    for (j in 1:length(outputNames)){
      if (outputConlist[[outputNames[j]]][[paste0('wrt_', fNms[k])]]){
        writeTable[[fNms[k]]] <- TRUE
      }
    }
  }
  return(list(TT=truthTable,
              outputOptions=outputOptions,
              writeTable=writeTable,
              plotTable=plotTable,
              dateInput=dateInput,
              dateOutput=dateOutput,
              delimI=delimI,
              delimO=delimO))
}

