

#' calculate vapor pressure from relative humidity
#'
#' @name VaporPressure
#' @param Temperature Air temperature in degrees C
#' @param RelativeHumidity Relative humidity as 0-1
#' @return Numerical vector of vapor pressures
#' @export


VaporPressure <- function ( Temperature,RelativeHumidity ){

  #Temperature in C, relative Humidity as decimal (1 > Rh > 0)

  #VaporPressure in mb

  satPressure <- SatVaporFromTemp( Temperature )*1000 #satP in mb

  VaporPressure <- RelativeHumidity*satPressure

  return(VaporPressure)

}

