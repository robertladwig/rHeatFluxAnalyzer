setwd(dirname(rstudioapi::getSourceEditorContext()$path))


library(dygraphs)
library(xts)
library(rLakeHeatFluxAnalyzer)

graphics.off()
rm(list = ls())
cat("\f")

rm(list=ls())

# read in example data
airT <- read.delim("../data/Esthwaite.airT")
rh <- read.delim("../data/Esthwaite.rh")
sw <- read.delim("../data/Esthwaite.sw")
wnd <- read.delim("../data/Esthwaite.wnd")
wtr <- read.delim("../data/Esthwaite.wtr")


# set latitude
lat <- 54


# combine data in a data.frame
dat <- data.frame(airT,RH=rh$RH,sw=sw$sw,wnd=wnd$wnd,W_temp=wtr$temp1)
# convert time to POSIXct format
dat$dateTime <- as.POSIXct(dat$dateTime,tz="UTC")



# test sw_alb
alb <- sw_alb(time = dat$dateTime,lat = 54)

# plot result
dygraph(xts(alb,dat$dateTime))


plot(dat$dateTime,alb,'l')

abline(h=0.02,col=2)



# set hight of wind measurement
hu <- 10
# test neutral_transfer_coeff
test_neut <- neutral_transfer_coeff(Uz = wnd$wnd,hu = hu)
# plot results
plot(dat$dateTime,test_neut$C_HN,'l')



# test clearSkySW
test_csSW <- clearSkySW(time = dat$dateTime,lat = lat,pressure = 1024,temperature = dat$airT,RH = dat$RH)
# plot results
plot(dat$dateTime,test_csSW,'l')



cls <- 1 - sw$sw/test_csSW
# test  Mdl_LW
test_Mdl_LW <- Mdl_LW(dat$airT,dat$RH,cls,as.numeric(format(dat$dateTime,"%m")))
# plot result
dygraph(xts(x = test_Mdl_LW,order.by = dat$dateTime))

