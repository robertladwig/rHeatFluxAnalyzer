

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