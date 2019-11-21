
## not sure if we even need this but:

#' calculate mean of a vector containing NaNs
#'
#' @name nanmean
#' @export
#' @param inputdata data for which to calculate the mean
#' @return arithmetic mean of inputdata with NaNs removed. Same as mean(x,na.rm=TRUE)
#' @seealso mean()
#' @examples
#' x <- c(1,3,1,4,NaN,4)
#' nanmean(x)

nanmean <- function( inputdata ){
  
  return(mean(inputdata,na.rm = TRUE))
  
}
