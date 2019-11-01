
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