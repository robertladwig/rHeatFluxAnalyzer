#' Model long wave radiation
#'
#' @name Mdl_LW
#' @param airT air temperature in degrees C
#' @param RH relative humidity in \%
#' @param clf cloud cover fraction
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
