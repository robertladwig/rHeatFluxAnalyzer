
#' @title Calculate the neutral transfer coefficient
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

  Uz <- complex(real = Uz)

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

  zot <- Re(zot)

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

  mm <- apply(mm,2,Re)

  return(mm)

}



