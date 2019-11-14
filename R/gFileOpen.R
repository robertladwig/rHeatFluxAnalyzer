gFileOpen <- function(fileName,treatAsWtr=FALSE){

  # author: Johannes Feldbauer 2019



  headers <- as.character(read.table(fileName,sep = "\t",nrows = 1,stringsAsFactors = FALSE))
  # format <-  '%s'
  # for (i in 2:length(regexp(headers,'\t','split'))){
  #   format <-  strcat(format,'%f')
  # }
  d  <- read.table(fileName,sep = "\t",header = TRUE)


  dat <-  d[,2:ncol(d)]
  dates <-  as.POSIXct(d[,1])

  if (treatAsWtr){
    heads <- unlist(regmatches(headers[2:length(headers)],
                               gregexpr("[0-9]+\\.*[0-9]*",headers[2:length(headers)])))
    depth <- as.numeric(heads)
    headers <- depth  # NOW are numbers...
  } else {
    headers <- NULL
    }

  return(list(dates=dates,dat=dat,depth=headers))
}
