# Make common functions for the development of the rLakeHeatFluxAnalyzer package

library(data.table)
library(stringr)

getExampleMeteo <- function(lakename="Esthwaite",path="../data"){

  airTemp=fread(paste0(path,"/",lakename,".airT"))
  relHum=fread(paste0(path,"/",lakename,".rh"))
  shortWave=fread(paste0(path,"/",lakename,".sw"))
  wind=fread(paste0(path,"/",lakename,".wnd"))

  metFile=Reduce(merge,list(wind,airTemp,relHum,shortWave))
  metFile[,dateTime:=as.POSIXct(dateTime)]
  return(metFile)
}

getExampleWaterT <- function(lakename="Esthwaite",path="../data"){
  waterTemp=fread(paste0(path,"/",lakename,".wtr"))
  waterTemp[,dateTime:=as.POSIXct(dateTime)]
  return(waterTemp)
}


getExampleConfig <- function(lakename="Esthwaite",path="../data"){
  filename=paste0(path,"/",lakename,".hfx")
  connection=file(filename,open = "r")
  configFile=readLines(connection)
  close(connection)

  return(configFile)
}

getExampleMatlabResults <- function(lakename="Esthwaite",path="../Matlab_ExampleFiles"){
  resultsFile=fread(paste0(path,"/",lakename,"_results.txt"))
  resultsFile[,DateTime:=as.POSIXct(DateTime)]
  return(resultsFile)
}

getValueFromConfig <- function(configFile,parameter){
  # Takes as argument the configFile as read by the "getExampleConfig" function

  full_string=configFile[str_which(configFile,parameter)]
  return(str_split(full_string,"\t\t")[[1]][1])
}
