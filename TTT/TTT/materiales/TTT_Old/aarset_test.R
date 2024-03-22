# Set working directory to file location
# -----------------------------------------------------------------------------
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


source("tttV10.R",        encoding="utf-8")
source("hypothesisV10.R", encoding="utf-8")

# These inputs are for each of the sizers test

#kernel  = "epanechnikov" # Which kernel do you want
kernel  = "gaussian"     
#kernel  = "triweight"    

#myMethod = "quadratic"
myMethod = "cubic"

variance = "bootstrap"
#variance = "moments"
#variance = "dirichlet"



weibullAditiveDataPath = file.path(paste(getwd(),"/data/weibullAditive.txt", sep = ""))
logDataPath            = file.path(paste(getwd(),"/data/logNormal.txt",      sep = ""))
ThreeWeibullsgDataPath = file.path(paste(getwd(),"/data/3weibull.txt",       sep = ""))


# Select whether we want a low time test (TRUE) or a complete test (FALSE) 

lowTest = TRUE

if(lowTest==TRUE){
  
  n       = 100             # Get how big is your data
  h       = 5              # How many windows do you have
  
  xgrid   = 50             # How many squares in the x axys
  ygrid   = h              # How many squares in the y axys
  
  
}else{
  
  n       = 1000           # Get how big is your data
  h       = 11             # How many windows do you have
  
  xgrid   = 401           # How many squares in the x axys
  ygrid   = h             # How many squares in the y axys
  
}

# These inputs are for hypothesis testing

repetitions   = 3                        # How many times do you want to repeat the tests (M=1000 probably, page 24 Table1.)
blueThreshold = 0.9                      # the proportion of blue pixels that should be in order to consider a test pass

passBlueTest  = rep(FALSE,repetitions)   # In here, we accumulate the result of each test

# For each test:
for (i in 1:repetitions) {
  
  # Xis
  # ---- Generate some random data
  
  #myData  = getRandomData(n,"weibull", 3, 1) 
  #myData  = getRandomData(n,"normal", 1, 0) 
  #myData  = getRandomData(n,"lognormal", 3, 1)
  #myData  = getRandomData(n,"exponential", 1) 
  #myData  = getRandomData(n,"addweibull", 3, 1)
  #myData  = getRandomData(n,"gamma", 3, 1) 
  #myData  = getRandomData(n,"gamma", 1/5, 5)
  
  # ---- Get data from file
  
  #myData = getDataFromFile(weibullAditiveDataPath)
  #myData = getDataFromFile(logDataPath)
  myData = getDataFromFile(ThreeWeibullsgDataPath)
 # myData = read.table('aarset.txt',h=T)
  
  # Normalize the data
  # myData = myData/mean(myData)
  
  
  
  # Run the TTT
  # ---- Expecific parameters
  
  # myResult <- ttt(myData, xgrid, ygrid, hOffset = 0.5,  kernel = "epanechnikov", myMethod = "quadratic", variance = "bootstrap",  ESSLimit = 5, saveCSV = TRUE, saveLog = FALSE)
  
  # ---- Default parameters
  
  myResult <- ttt(myData, xgrid, ygrid, saveCSV = TRUE, saveLog = FALSE)
  
  #TODO: All warnings are because there are missing points in the
  #      family plots or the phi vector. Due to having NA values
  #      when we do the cramer system divided by 0. Just need
  #      to silent all warnings when everything is complete
  #      and there is no other type of warning
  
  
  # Count the number of pixels
  myHypothesis <- hypothesisColors(myResult)
  
  
  # Save the result for later and repeat
  # -- Zero  Sizer
  #    (Save nothing from the zero SiZer)
  
  # -- First Sizer
  # ---- Red
  # ---- Purple
  # ---- Blue
  #    (Save only the blue information from the first SiZer)
  if(myHypothesis$percentage[6] > blueThreshold ){
    
    passBlueTest[i] = TRUE
    
  }
  
  # -- Second Sizer
  #    (Save nothing from the second SiZer)
  
}




# Repetir el experimento R veces y sacer el summary the cada hipotesis
summary(passBlueTest)
