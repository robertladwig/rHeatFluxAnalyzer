##---------------------- sw_albedo.m ---------------------------------------------------------------

#' albedo for short wave radiation for given latitude and time of year
#'
#' @name sw_albedo
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

sw_albedo <-  function(time,lat){





  ## convert time to day of the year

  # using base R

  doy <- as.numeric(format(time,"%j")) + (as.numeric(format(time,"%H"))*3600 +

                                            as.numeric(format(time,"%M"))*60 +

                                            as.numeric(format(time,"%S")))/86400

  # using lubridate

  #doy2 <- as.numeric(format(time,"%j")) + (hour(time)*3600 + minute(time)*60 +
  #                                             second(time))/86400



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

