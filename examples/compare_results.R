rm(list = ls())
graphics.off()
cat("\f")

res_m <- read.table("../Matlab_ExampleFiles/Esthwaite_results.txt",header = TRUE,sep = "\t")
res_m$DateTime <- as.POSIXct(res_m$DateTime)

res_r <- read.table("../data/Esthwaite_results.txt",header = TRUE,sep = "\t")
res_r$DateTime <- as.POSIXct(res_r$DateTime)

names <- read.table("../data/Esthwaite_results.txt",nrows = 1,stringsAsFactors = FALSE)
res_d <- res_m - res_r

for (i in 1:34) {

plot(res_m$DateTime,res_m[,i],'l',main = as.character(names[i]))
lines(res_r$DateTime,res_r[,i],col=2)
}
