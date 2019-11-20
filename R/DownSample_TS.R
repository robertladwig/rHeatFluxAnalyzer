#' Down sample time series
#'
#' @param Dates Vector of original data times
#' @param outRs Sampling inerval for output time steps
#' @param VarArray Array of original data
#' @export

DownSample_TS <- function(Dates,outRs,VarArray=NULL ){
  #----Author: Jordan S Read, 2009----
    #Dates in matDate format: (datenum('7/1/2009')) will be in matDate format
  #outRs: seconds. If outRs < step between dates, original values are
  #returned.
  #if some dates are not uniuque, values will be averaged***
    #

  if(length(VarArray)==0){
   DS_varArray <- NULL
  }

  matStep <- outRs#/86400 #matdates based on day unit, conver seconds/day
  u_matStep <- 1/matStep
  sortedDates <- Dates # assumes values passed are pre-sorted 

  tempDates <- floor(as.numeric(sortedDates)*u_matStep)*matStep

  uniStInd <- which(!duplicated(tempDates))
  DS_dates <- Dates[uniStInd]
  if (length(VarArray!=0)){
    DS_Length <- length(DS_dates)
    DS_varArray <- rep(NA,DS_Length)
    start <- 1
    for (j in 1:DS_Length){
      matchInd <- start:uniStInd[j]
      if (uniStInd[j]>length(VarArray) || length(VarArray[matchInd])==0){
        DS_varArray[j] <- NA
      } else {
        DS_varArray[j] <- mean(VarArray[matchInd])
      }
      start <- uniStInd[j]+1
    }
  }
  return(list(DS_dates=DS_dates,
              DS_varArray=DS_varArray ))
}
