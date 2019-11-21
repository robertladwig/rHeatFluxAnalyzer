rm(list = ls())
graphics.off()
cat("\f")

res_m <- read.table("../Matlab_ExampleFiles/Esthwaite_results.txt",header = TRUE,sep = "\t")
res_m$DateTime <- as.POSIXct(res_m$DateTime)

res_r <- read.table("../data/Esthwaite_results.txt",header = TRUE,sep = "\t")
res_r$DateTime <- as.POSIXct(res_r$DateTime)


