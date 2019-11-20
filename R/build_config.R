build_config <- function(LakeName,directory){

done <- FALSE
spcFrac<- 2 #number of spaces <- char
num2delim <- 20
defSlt <- 'wTemp'
# see if file exists

fI <- file.exists(paste0(directory, '/', LakeName, '.hfx'))
if (fI){
  cat(paste0("File '",LakeName, ".hfx' allready exists.\n"))
  create <- readline(paste0("Create new file? Alternatively use ",LakeName, ".hfx"," (Y/N) \n"))
  if(!create%in%c("y","Y","yes","YES","Yes")){
    cat("Done \n")
  }
} else {
  # define possible outputs
  outputOptions <- matrix(c('Momentum flux','tau',
    'Sensible heat flux','Qh','Latent heat flux','Qe',
    'Transfer coefficient for momentum (neutral)','C_DN',
    'Transfer coefficient for heat (neutral)','C_EN',
    'Transfer coefficient for humidity (neutral)','C_HN',
    'Transfer coefficient for momentum at 10 m(neutral)','C_D10N',
    'Transfer coefficient for heat at 10 m (neutral)','C_E10N',
    'Transfer coefficient for humidity at 10 m (neutral)','C_H10N',
    'Transfer coefficient for momentum','C_D',
    'Transfer coefficient for heat','C_E',
    'Transfer coefficient for humidity','C_H',
    'Transfer coefficient for momentum at 10 m','C_D10',
    'Transfer coefficient for heat at 10 m','C_E10',
    'Transfer coefficient for humidity at 10 m','C_H10',
    'Wind speed at 10 m','u10',
    'Wind speed at 10 m (neutral)','u10N',
    'Air temperature at 10 m','t10',
    'Relative humidity at 10 m','rh10',
    'Net long wave radiation','Qlnet',
    'Incoming long wave radiation','Qlin',
    'Outgoing long wave radiation','Qlout',
    'Air shear velocity','uSt_a',
    'Air shear velocity (neutral)','uSt_aN',
    'Evaporation','Evap',
    'Water temperature','wTemp',
    'Reflected Short wave radiation','Qsr',
    'Atmospheric stability','obu',
    'Total surface heat flux','Qtot',
    'Short wave radiation','Qs',
    'Net incoming short wave radiation','Qsin',
    'Air density at 10m','rhoa10',
    'Water density','rhow',
    'Air density','rhoa'),34,2,byrow = TRUE)

  # defaults for parameters
  defaults <- c('86400','2','2','2','54','0','98','0.2','Y','Y')
  # long names for outputs
  names    <- c('output resolution (s)','wind height (m)',
    'temperature height (m)','humidity height (m)',
    'latitude (°N)','altitude (m)',
    'max wind speed (m/s)','min wind speed (m/s)',
    'plot figure (Y/N)','write results (Y/N)')

  use_def <- readline("Use default configurations? (Y/N)")

  if(use_def%in%c("y","Y","yes","YES","Yes")){
    user_input <- defaults
    output_choose <- 1:28
  } else {
    user_input <- rep("",10)
    for (i in 1:10) {

      cat(paste0(names[i],". (default value ",defaults[i],") \n"))
      user_input[i] <- readline("input value: ")
    }

    output_choose <- rep(FALSE,32)

    for (i in 1:32) {

      cat(paste0("Output ",outputOptions[i,1],"? \n"))
      temp <- readline("(Y/N)")
      output_choose[i] <- ifelse(temp=="Y",TRUE,FALSE)
    }

  }

  Output_user <- outputOptions[output_choose,2]

  file_out <- paste0("Configuration file for ",LakeName,"\n \n",
                     paste0(Output_user,collapse = ","),"\t\t#outputs \n",
                     paste0(apply(matrix(c(user_input,paste0(" \t\t # ",names)),2,10,byrow = TRUE),2,
                                  paste0,collapse=""),collapse = "\n"))
  cat(x = file_out,file=paste0(directory, '/', LakeName, '.hfx'))
  }
}
