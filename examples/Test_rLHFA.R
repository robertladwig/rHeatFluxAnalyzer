# Test rLHFA functions by comparing to the MATLAB output

library(ggplot2)

Sys.setenv(tz="GMT")

source("CommonFunctions.R") # Path to CommonFunctions.R

example_wtr=getExampleWaterT()
example_meteo=getExampleMeteo()
example_config=getExampleConfig()
example_results=getExampleMatlabResults()


setwd("../R") # setwd to the "R" folder on github

# Test the rLHFA
# Order needs to be like this, because some functions require other functions in the package

# Saturated vapor pressure from temperature
source("SatVaporFromTemp.R")
rLHFA_SatVapor=SatVaporFromTemp(Temperature = example_meteo$airT)
# No MATLAB output available


# Smith Gamma (?)
source("getSmithGamma.R")
rLHFA_smith_gamma=getSmithGamma(lat=as.numeric(getValueFromConfig(example_config,"latitude")),
                                time = example_meteo$dateTime)
# No MATLAB output available


# Clear-sky short-wave radiation
source("clearSkySW.R")
rLHFA_SkySW=clearSkySW(time=example_meteo$dateTime,lat=as.numeric(getValueFromConfig(example_config,"latitude")))
# No MATLAB output available

ggplot()+
  geom_line(aes(x=example_meteo$dateTime,y=rLHFA_SkySW))
max(rLHFA_SkySW)


source("VaporPressure.R")
rLHFA_VapPres=VaporPressure(Temperature = example_meteo$airT,RelativeHumidity = example_meteo$RH)
# No MATLAB output available


# Net long-wave radiation
source("Mdl_LW.R")
source("calc_lwnet.R")
rLHFA_lw=calc_lwnet(time=example_meteo$dateTime,lat=as.numeric(getValueFromConfig(example_config,"latitude")),
                       press=NULL, ta=example_meteo$airT,rh=example_meteo$RH,sw=example_meteo$sw,
                       ts=example_wtr[[2]])
rLHFA_lwnet=rLHFA_lw$lwnet
rLHFA_lwin=rLHFA_lw$lw
rLHFA_lwout=rLHFA_lw$LWo
MATLAB_lwnet=example_results$`Qlnet (W m^{-2})`
MATLAB_lwin=example_results$`Qlin (W m^{-2})`
MATLAB_lwout=example_results$`Qlout (W m^{-2})`

df_LW=rbindlist(list(data.table("time"=example_meteo$dateTime, "LW_net"=rLHFA_lwnet,
                                "LW_in"=rLHFA_lwin,"LW_out"=rLHFA_lwout,"version"="rLHFA"),
                     data.table("time"=example_meteo$dateTime, "LW_net"=MATLAB_lwnet,
                                "LW_in"=MATLAB_lwin,"LW_out"=MATLAB_lwout,"version"="MATLAB")))
ggplot(data=df_LW)+
  geom_line(aes(time,LW_net,colour=version))+
  labs(subtitle=paste("Sum of squares of the error:",sum(df_LW$LW_net[version=="rLHFA"]-df_LW$LW_net[version=="MATLAB"])))
ggplot(data=df_LW)+
  geom_line(aes(time,LW_in,colour=version))+
  labs(subtitle=paste("Sum of squares of the error:",sum(df_LW$LW_in[version=="rLHFA"]-df_LW$LW_in[version=="MATLAB"])))
ggplot(data=df_LW)+
  geom_line(aes(time,LW_out,colour=version))+
  labs(subtitle=paste("Sum of squares of the error:",sum(df_LW$LW_out[version=="rLHFA"]-df_LW$LW_out[version=="MATLAB"])))


ggplot(data=df_LW[time<as.POSIXct("2009-01-03")])+
  geom_line(aes(time,LW_in,colour=version))


