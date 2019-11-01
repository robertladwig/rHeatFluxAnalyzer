##---------------------- sw_albedo.m ---------------------------------------------------------------

#' albedo for short wave radiation for given latitude and time of year
#'
#' @name sw_alb
#' @param lat latitude in degrees North.
#' @param time date as a POSIXct object in UTC time zone(???)
#' @return numeric vector of albedo for short wave radiation for given latitude and time of year
#' @export
#' @references science and dude, 2001 - some title
#' @seealso test
#' @examples
#'  lat <- 54
#'  time <- as.POSIXct("2001-10-10",tz="UTC")
#'  sw_alb(lat,time)

sw_alb <-  function(time,lat){





  ## convert time to day of the year

  # using base R

  doy <- as.numeric(format(time,"%j")) + (as.numeric(format(time,"%H"))*3600 +

                                            as.numeric(format(time,"%M"))*60 +

                                            as.numeric(format(time,"%S")))/86400

  # using lubridate

  doy2 <- as.numeric(format(time,"%j")) + (hour(time)*3600 + minute(time)*60 +

                                             second(time))/86400



  # Assign Parameter Values #

  Tropic <-  23.45

  YrDays <-  365.0

  DclDay <-  284.0

  DgCrcl <-  360.0

  DgToRd <-  57.29577951

  RefInd <-  1.33 # This is a mild function of T and S. http://scubageek.com/articles/wwwh2o.html



  Lat <-  lat/DgToRd

  z <-  DclDay + floor(doy)

  x <-  DgCrcl*z/YrDays/DgToRd

  y <-  sin(x)

  Decl <-  Tropic/DgToRd*y



  # Hour-angle calculation where time is hours from midnight  #

  HrAng <-  ((doy-floor(doy))*24 -12)*15.0/DgToRd



  # Zenith angle calculation #

  Zenith <-  acos(sin(Decl)*sin(Lat)+cos(Decl)*cos(Lat)*cos(HrAng))



  # Angle of Refraction calculation based on Snell's Law #

  RefAng <-  asin(sin(Zenith)/RefInd) # Angle of refraction



  # Albedo Calculation #

  A1 <-  tan(Zenith-RefAng)^2

  A2 <-  tan(Zenith+RefAng)^2

  A3 <-  sin(Zenith-RefAng)^2

  A4 <-  sin(Zenith+RefAng)^2

  sw_alb <-  0.5 * (A1/A2 + A3/A4)



  # Set albedo to 1 if greater than 1 #

  sw_alb[sw_alb > 1] <-  1

  return(sw_alb)

}





##------------------------ neutral_transfer_coeff.m ------------------------------------------------

#' Calculate the neutral transfer coefficient
#'
#' @name neutral_transfer_coeff
#' @param Uz vector of wind speed in (m/s)
#' @param hu height of the wind speed measurement in m
#' @return A data.frame containing ustar, u10, C_D10N, C_E10N, C_H10N, and C_DN
#' @export
#' @examples
#' neutral_transfer_coeff(6,10)


neutral_transfer_coeff <- function (Uz,hu){



  # Uz: wind speed, m/s.

  # hu: height of wind measurement, m.



  # define constants

  const_vonKarman <-  0.41 # von Karman constant

  const_Charnock <-  0.013 # charnock constant

  const_Gravity <-  9.81 # gravitational acceleration, m/s2



  # kinematic viscosity, m2 s-1

  KinV <-  1.5e-5



  # estimate initial values of u* and zo

  Uz[Uz == 0] <-  0.001 # need wind to be > 0, otherwise equation breaks down see Amorocho and Devries (1980)

  ustarN <-  Uz*sqrt(0.00104+0.0015/(1+exp((-Uz+12.5)/1.56)))

  zo <-  (const_Charnock*ustarN^2/const_Gravity) + (0.11*KinV/ustarN)

  zo_prev <-  zo*1.1

  for (i in 1:length(Uz)){

    while (abs((zo[i] - zo_prev[i]))/abs(zo_prev[i]) > 0.00001){

      ustarN[i] <-  const_vonKarman*Uz[i]/(log(hu/zo[i]))

      dummy <-  zo[i]

      zo[i] <- (const_Charnock*ustarN[i]^2/const_Gravity) + (0.11*KinV/ustarN[i])

      zo_prev[i] <-  dummy

    }

  }



  # calculate neutral transfer coefficients

  C_DN <- (ustarN^2)/(Uz^2)

  re <- ustarN*zo/KinV

  zot <- zo*exp(-2.67*(re)^(0.25) + 2.57)

  #zot <- real(zot)

  zoq <- zot

  C_HN <- const_vonKarman*sqrt(C_DN)/(log(hu/zot))

  C_EN <- C_HN



  # calculate neutral transfer coefficients at 10 m

  C_D10N <- (const_vonKarman/log(10/zo))*(const_vonKarman/log(10/zo))

  C_E10N <- (const_vonKarman*const_vonKarman)/(log(10/zo)*log(10/zoq))

  C_H10N <- C_E10N



  # calculate wind speed at 10 m following Amorocho and Devries (1980)

  if (hu != 10){

    u10 <- Uz/(1-sqrt(C_DN)/const_vonKarman*log(10/hu))

  } else {

    u10 <- Uz

  }



  # define output matrix

  mm <- data.frame(ustarN = ustarN, u10 = u10, C_D10N = C_D10N,

                   C_E10N = C_E10N, C_H10N = C_H10N, C_DN = C_DN,

                   C_EN = C_EN,C_HN = C_HN)

  #mm <- real(mm)

  return(mm)

}






##------------------------ calc_lwnet.m ------------------------------------------------------------


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


#' Calculate net long wave radiation
#'
#' @name calc_lwnet
#' @param time Vector of times as POSIXct
#' @param lat latitude for which to calculate clear sky radiation
#' @param press Air pressure in (hPa) if available
#' @param ta Air temperature in (°C) if available
#' @param rh relative humidity in percentage if available
#' @param sw observed showrt wave radiation
#' @param ts Water surface temperature in °C
#' @return A data.frame with lw, LWo, and lwnet
#' @export
#' @references Read et al. 2012


calc_lwnet <- function(time,lat,press,ta,rh,sw,ts){



  clrSW <- clearSkySW(time,lat,press,ta,rh*.01)

  clf <-  1-sw/clrSW



  # impose limits

  ltI <- clf<0

  clf[ltI] <- 0

  gtI <- clf>1

  clf[gtI] <- 1



  ## needs translation !!

  # when short-wave is zero make cc = day time average

  #dateV <- datevec(time)

  #[~,~,b] <- unique(dateV(:,1:3),'rows')

  b <- as.numeric(format(time,"%d"))

  daily_av <- aggregate(clf,by=list(b),mean, na.rm=TRUE)

  c2 <- clf

  for (ii in 1:length(unique(b))){

    c2[b == ii] <- daily_av[ii]

  }

  clf[sw == 0] <- c2[sw == 0]



  # if still NaN, assume cloud cover is zero (at night, same as Read et al. 2012)

  clf[is.nan(clf)]  <-  0



  #month <- datevec(Jday)

  month <- as.numeric(format(time,"%m"))

  lw <- Mdl_LW(ta,rh*0.01,clf,month)

  Tk <- ts + 273.13     # sample Ts to match LW

  emiss <- 0.972

  S_B <- 5.67E-8

  LWo <- S_B*emiss*Tk^4

  lwnet <- lw-LWo



  return(data.frame(lw = lw,LWo = LWo,lwnet = lwnet))



}

#' Calculate saturation vapor pressure from temperature
#'
#' @name SatVaporFromTemp
#' @param Temperature Water temperature
#' @return Vector of saturation vapor pressures in bar
#' @export
#' @references Empirical fit from P.R. Lowe, 1976. An Approximating Polynomial for the
#' Computation of Saturation Vapor Pressure. Journal of Applied Meteorology
#' @examples
#' SatVaporFromTemp(15)

SatVaporFromTemp <- function( Temperature ){

  #Empirical fit from P.R. Lowe, 1976. An Approximating Polynomial for the

  #Computation of Saturation Vapor Pressure. Journal of Applied Meteorology

  Temp <- Temperature

  a0 <- 6.107799961

  a1 <- 4.436518521E-1

  a2 <- 1.428945805E-2

  a3 <- 2.650648471E-4

  a4 <- 3.031240396E-6

  a5 <- 2.034080948E-8

  a6 <- 6.136820929E-11

  Pressure <- a0 + Temp*(a1+Temp*(a2+Temp*(a3+Temp*(a4+Temp*(a5*a6*Temp)))))

  Pressure <- Pressure/1000 #now in bar

  return(Pressure)

}

#' Function to get Smith Gamma
#'
#' @name getSmithGamma
#' @usage
#' @param lat Latitude
#' @param time time as POSIXct format
#' @return Smiths Gamma
#' @export
#' @examples
#'  lat <- 54
#'  time <- as.POSIXct("2001-10-10",tz="UTC")
#'  getSmithGamma(lat,time)

getSmithGamma <- function(lat,time){



  doy <- as.numeric(format(time,"%j")) + (as.numeric(format(time,"%H"))*3600 +

                                            as.numeric(format(time,"%M"))*60 +

                                            as.numeric(format(time,"%S")))/86400



  if(lat<0){

    doy <- doy+182    # convert for S Hemisphere

    lat <- abs(lat)

  }

  if(length(lat)<length(time)){

    lat <- rep(lat,length(time))

  }

  yr <- as.numeric(format(time,"%Y"))

  #yr <- yr[,1]

  #time <- time-datenum(yr,0,1)

  latDeg <- seq(5,85,by=10)

  season <- c(-10, 81, 173, 264, 355, 446)  # winter,spring,summer,fall,winter,spring

  ## need to check if matrix is correct!!

  gammaTable <- matrix(data = c(3.37, 2.85, 2.8,  2.64, 3.37, 2.85,

                                2.99, 3.02, 2.7,  2.93, 2.99, 3.02,

                                3.6,  3,    2.98, 2.93, 3.6,  2.98,

                                3.04, 3.11, 2.92, 2.94, 3.04, 3.11,

                                2.7,  2.95, 2.77, 2.71, 2.7,  2.95,

                                2.52, 3.07, 2.67, 2.93, 2.52, 3.07,

                                1.76, 2.69, 2.61, 2.61, 1.76, 2.69,

                                1.6,  1.67, 2.24, 2.63, 1.6,  1.67,

                                1.11, 1.44, 1.94, 2.02, 1.11, 1.44),9,6,byrow = TRUE)

  ## needs translation!! --> package pracma has the same function

  SmithGamma <- pracma::interp2(latDeg,season,t(gammaTable),lat,doy)

  return(SmithGamma)

}





##-------------------------------- Mdl_LW.m --------------------------------------------------------

#' Model long wave radiation
#'
#' @name Mdl_LW
#' @param airT air temperature in °C
#' @param RH relative humidity in %
#' @param clf ??
#' @param month month of the year
#' @return numerical vector of long wave radiation
#' @export
#' @examples
#' Mdl_LW(20,50,0.1,12)

Mdl_LW <- function(airT,RH,clf,month){

  S_B <- 5.67E-8

  vp <- VaporPressure(airT,RH)   # air vapor pressure (mbars)

  T_k <- 273.13+airT

  cl1 <- 1-clf

  cl2 <- 1.22+0.06*sin((month+2)*pi/6)

  cl3 <- vp/T_k

  T_k1 <- T_k^4

  LW <- T_k1*(clf+cl1*cl2*cl3^(1/7))*S_B

  return(LW)

}

#' calculate vapor pressure from relative humidity
#'
#' @name VaporPressure
#' @param Temperature Air temperature in °C
#' @param RelativeHumidity Relative humidity as 0-1
#' @return Numerical vector of vapor pressures
#' @export


VaporPressure <- function ( Temperature,RelativeHumidity ){

  #Temperature in °C, relative Humidity as decimal (1 > Rh > 0)

  #VaporPressure in mb

  satPressure <- SatVaporFromTemp( Temperature )*1000 #satP in mb

  VaporPressure <- RelativeHumidity*satPressure

  return(VaporPressure)

}




##-------------------------------- nanmean.m -------------------------------------------------------



## not sure if we even need this but:

#' calculate mean of a vector containing NaNs
#'
#' @name nanmean
#' @param inputdata data for which to calculate the mean
#' @return arithmetic mean of inputdata with NaNs removed. Same as mean(x,na.rm=TRUE)
#' @seealso mean()
#' @export
#' @examples
#' x <- c(1,3,1,4,NaN,4)
#' nanmean(x)

nanmean <- function( inputdata ){

  return(mean(inputdata,na.rm = TRUE))

}
