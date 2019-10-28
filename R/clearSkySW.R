
#' Calculate the clear sky short wave radiation
#'
#' @name clearSkySW
#' @param time Vector of times as POSIXct
#' @param lat latitude for which to calculate clear sky radiation
#' @param pressure Air pressure in (hPa) if available
#' @param temperature Air temperature in (°C) if available
#' @param RH relative humidity in percentage if available
#' @return A numerical vecot of clear sky short wave radiation
#' @export
#' @references Crawford and Duchon, 1991

clearSkySW <- function(time,lat,pressure=NULL,temperature=NULL,RH=NULL){
  
  
  
  if (is.null(temperature)){
    
    temperature <- 20
    
  }
  
  if (is.null(pressure)){
    
    pressure <- 1013
    
  }
  
  if (is.null(RH)){
    
    RH       <- .7
    
  } else{
    
    # this should only be done if RH is given in %
    
    RH <- RH*0.01
    
  }
  
  ## resize vectors
  
  #time <- reshape(time,length(time),1)
  
  if (length(time)>length(lat)){
    
    lat  <- rep(1,length(time))*lat[1]
    
  }
  
  if (length(time)>length(pressure)){
    
    pressure  <- rep(1,length(time))*pressure[1]
    
  }
  
  if (length(time)>length(temperature)){
    
    temperature  <- rep(1,length(time))*temperature[1]
    
  }
  
  if (length(time)>length(RH)){
    
    RH  <- rep(1,length(time))*RH[1]
    
  }
  
  
  
  # # integrating lat and long into a timezone calculation...
  
  # URL <- ['http://www.askgeo.com/api/334001/jk6n0p5jknsbe5gq0dh1su5bsr/'...
  
  #     'timezone.json?points=' num2str(lat) '#2C' num2str(long)]
  
  # str <- urlread(URL)
  
  # currentStI <- strfind(str,'"currentOffsetMs"')
  
  # latitudeStI <- strfind(str,',"latitude"')
  
  # GMT_Off <- str2double(str(currentStI+18:latitudeStI)) # in ms, UTC
  
  # hr_Off <- GMT_Off/3600/1000     # hour offset from GMT
  
  # calculate SW from clear sky
  
  # from Crawford and Duchon 1991 *** key reference here ***
  
  t_noon <- 12.5          # solar noon (actually will depend on timezone)
  
  #tm  <- datevec(time)
  
  #n  <- floor(time-datenum(tm[,1],1,0))  # Julian day
  
  n <- as.numeric(format(time,"%j")) + (as.numeric(format(time,"%H"))*3600 +
                                          
                                          as.numeric(format(time,"%M"))*60 +
                                          
                                          as.numeric(format(time,"%S")))/86400
  
  
  
  tm  <-  as.numeric(format(time,"%H"))           # time (hours)
  
  
  
  cosN <- 1+0.034*cos(2*pi*(n-1)/365)
  
  I0 <- 1353*cosN*cosN     # Meyers & Dale 1983
  
  
  
  H <- (pi/12)*(t_noon-tm)   # t_noon is solar noon, t is the local solar time
  
  # e.g. t = 12.5 and H = 0 at local solar noon
  
  
  
  sigma <- 180/pi*(0.006918-.399912*cos(2*pi/365*(n-1))+0.070257*
                     
                     sin(2*pi/365*(n-1))-.006758*cos(4*pi/365*(n-1))+
                     
                     0.000907*sin(4*pi/365*(n-1))-0.002697*cos(6*pi/365*(n-1))+
                     
                     0.00148*sin(6*pi/365*(n-1)))
  
  
  
  # - - cosine of solar zenith angle - -
  
  sin1 <- sin(lat*2*pi/360)
  
  sin2 <- sin(sigma*2*pi/360)
  
  cos1 <- cos(lat*2*pi/360)
  
  cos2 <- cos(sigma*2*pi/360)
  
  cos3 <- cos(H)
  
  cosZ <- sin1*sin2+cos1*cos2*cos3
  
  
  
  p <- pressure                       # in millibars
  
  m1 <- 1224*cosZ*cosZ+1
  
  m <- 35*m1^-.5                     # optical air mass at p = 1013 mb
  
  
  
  Tr1 <- m*(0.000949*p+0.051)
  
  T_rT_pg <- 1.021-0.084*Tr1^.5      # Atwater & Brown 1974
  
  
  
  Es <- SatVaporFromTemp(temperature)*1000
  
  E  <- (RH*Es)
  
  Td1 <- 243.5*log(E/6.112)
  
  Td2 <- 17.67-log(E/6.112)
  
  T_d <- Td1/Td2                     # dewpoint (degrees c)
  
  T_d <- T_d*9/5+32                   # dewpoint (degrees F)
  
  #lat=46
  
  G   <- getSmithGamma(lat,time)      # empirical constant dependent upon
  
  # time of the year and latitude (Smith
  
  # 1966 tale 1
  
  T_a <- m
  
  for (i in 1:length(m)){
    
    T_a[i] <- 0.935^m[i]            # Meyers & Dale 1993
    
  }
  
  
  
  mu <- exp(0.1133-log(G+1)+0.0393*T_d) # precipitable water
  
  
  
  mu1 <- mu*m
  
  T_w <- 1-.077*mu1^0.3
  
  
  
  ISW <- I0*cosZ*T_rT_pg*T_w*T_a
  
  useI <- ISW<0
  
  clrSW <- rep(0,length(time))
  
  clrSW[useI] <- ISW[useI]
  
  return(clrSW)
  
}

