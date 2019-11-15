#' function to read in config file for rLHFA
#'
#' @param LakeName Name of the Lake. Forcing and config file must have the same name
#' @param folder Folder where the config file is stored
#' @return Returns a list with options from the config file

OpenCfg <- function(LakeName,folder){
 # original Matlab Author: Jordan S Read 2009
 # translated to R from Johannes Feldbauer 2019


  fileName  <-  paste0(folder, '/', LakeName ,'.hfx')

  Test <- gsub("\t","",readLines(fileName))
  Test <- gsub("#.*","",Test)

  #TPuts <-  sub("\t","",sub(" ","",unlist(strsplit(x = readLines(fID,n = 3)[3],split = ","))))
  tOut <- unlist(strsplit(x = Test[3],split = ","))
  outPuts <- gsub(" ","",tOut)

  outRs <- as.numeric(Test[4])
  wndH  <- as.numeric(Test[5])
  htH   <- as.numeric(Test[6])
  hqH   <- as.numeric(Test[7])
  lat   <- as.numeric(Test[8])
  alt   <- as.numeric(Test[9])
  wndMx <- as.numeric(Test[10])
  wndMn <- as.numeric(Test[11])
  plotYes <- FALSE
  writeYes  <-  FALSE
  if (Test[12]%in%c("y","Y","yes","Yes")){
    plotYes <-  TRUE
  }
  if (Test[13]%in%c("y","Y","yes","Yes")){
    writeYes <-  TRUE
  }
  # function for getting number for single line entry. Ignores spaces
  return(list(outPuts=outPuts,outRs=outRs,wndH=wndH,hqH=hqH,htH=htH,lat=lat,
              alt=alt,wndMx=wndMx,wndMn=wndMn,plotYes=plotYes,writeYes=writeYes))
}
