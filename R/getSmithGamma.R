#' Function to get Smith Gamma
#'
#' @name getSmithGamma
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



  if(min(lat)<0){

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


