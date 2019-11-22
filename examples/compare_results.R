rm(list = ls())
graphics.off()
cat("\f")

wnd <- gFileOpen("../data/Esthwaite.wnd")

res_m <- read.table("../Matlab_ExampleFiles/Esthwaite_results.txt",header = TRUE,sep = "\t")
res_m$DateTime <- as.POSIXct(res_m$DateTime)

res_r <- read.table("../data/Esthwaite_results.txt",header = TRUE,sep = "\t")
res_r$DateTime <- as.POSIXct(res_r$DateTime)

names <- read.table("../data/Esthwaite_results.txt",nrows = 1,stringsAsFactors = FALSE)
res_d <- res_m - res_r

summary(res_d)

for (i in 2:35) {
par(mfrow=c(2,1))
plot(res_m$DateTime,res_m[,i],'l',main = as.character(names[i]))
lines(res_r$DateTime,res_r[,i],col=2)
plot(res_m$DateTime,res_d[,i]/res_m[,i],main = "rel. difference",'l')
}
