# Test rLHFA functions by comparing to the MATLAB output

library(ggplot2)

Sys.setenv(tz="GMT")

source("CommonFunctions.R") # Path to CommonFunctions.R

example_wtr=getExampleWaterT()
example_meteo=getExampleMeteo()
example_config=getExampleConfig()
example_results=getExampleMatlabResults()


setwd("../R") # setwd to the "R" folder on github

# Load all the scripts
for(i in list.files()){
  source(i)
}

# Test the rLHFA
# Order needs to be like this, because some functions require other functions in the package

# Saturated vapor pressure from temperature
rLHFA_SatVapor=SatVaporFromTemp(Temperature = example_meteo$airT)
# No MATLAB output available


# Smith Gamma (?)
rLHFA_smith_gamma=getSmithGamma(lat=as.numeric(getValueFromConfig(example_config,"latitude")),
                                time = example_meteo$dateTime)
# No MATLAB output available


# Clear-sky short-wave radiation
rLHFA_SkySW=clearSkySW(time=example_meteo$dateTime,lat=as.numeric(getValueFromConfig(example_config,"latitude")))
# No MATLAB output available

ggplot()+
  geom_line(aes(x=example_meteo$dateTime,y=rLHFA_SkySW))
max(rLHFA_SkySW)


rLHFA_VapPres=VaporPressure(Temperature = example_meteo$airT,RelativeHumidity = example_meteo$RH)
# No MATLAB output available


# Net long-wave radiation
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
  labs(subtitle=paste("Sum of squares of the error:",sum(abs(df_LW$LW_net[df_LW$version=="rLHFA"]-df_LW$LW_net[df_LW$version=="MATLAB"]))))
ggplot()+
  geom_line(aes(x=example_meteo$dateTime,y=df_LW$LW_net[df_LW$version=="rLHFA"]-df_LW$LW_net[df_LW$version=="MATLAB"]))+
  labs(title="Difference MATLAB and R output over time")

ggplot(data=df_LW)+
  geom_line(aes(time,LW_in,colour=version))+
  labs(subtitle=paste("Sum of squares of the error:",sum(abs(df_LW$LW_in[df_LW$version=="rLHFA"]-df_LW$LW_in[df_LW$version=="MATLAB"]))))
ggplot()+
  geom_line(aes(x=example_meteo$dateTime,y=df_LW$LW_in[df_LW$version=="rLHFA"]-df_LW$LW_in[df_LW$version=="MATLAB"]))+
  labs(title="Difference MATLAB and R output over time")

ggplot(data=df_LW)+
  geom_line(aes(time,LW_out,colour=version))+
  labs(subtitle=paste("Sum of squares of the error:",sum(abs(df_LW$LW_out[df_LW$version=="rLHFA"]-df_LW$LW_out[df_LW$version=="MATLAB"]))))
ggplot()+
  geom_line(aes(x=example_meteo$dateTime,y=df_LW$LW_out[df_LW$version=="rLHFA"]-df_LW$LW_out[df_LW$version=="MATLAB"]))+
  labs(title="Difference MATLAB and R output over time")


# Latent heat flux
rLHFA_Evap=sens_latent(example_wtr$temp1,example_meteo$wnd,example_meteo$airT,example_meteo$RH,
                       hu=2,ht=2,hq=2,alt=0,lat=54)$Evap
MATLAB_Evap=example_results$`Evap (mm day^{-1})`

df_Evap=rbindlist(list(data.table("time"=example_meteo$dateTime, "Evap"=rLHFA_Evap,"version"="rLHFA"),
                     data.table("time"=example_meteo$dateTime, "Evap"=MATLAB_Evap,"version"="MATLAB")))

ggplot(data=df_Evap)+
  geom_line(aes(time,Evap,colour=version))+
  labs(subtitle=paste("Sum of squares of the error:",sum(abs(df_Evap$Evap[df_Evap$version=="rLHFA"]-df_Evap$Evap[df_Evap$version=="MATLAB"]))))
ggplot()+
  geom_line(aes(x=example_meteo$dateTime,y=df_Evap$Evap[df_LW$version=="rLHFA"]-df_Evap$Evap[df_LW$version=="MATLAB"]))+
  labs(title="Difference MATLAB and R output over time")
