# Set working directory to file location
# -----------------------------------------------------------------------------
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

source("tttV11.R",        encoding="utf-8")

#' For a given array of data, creates several plots describing ....
#'
#' @param dataPath        where is the data you want to use for the tests
#'
#' @param testRepetition  array with how many repetitions you want for each test.
#'                        i.e.: [10, 100, 1000]
#'                        
#' @param blueThreshold   float between 0 and 1. Each repetetion is considered suscessfull
#'                        if the proportion of blue pixels in the SiZer maps is greater or
#'                        equal to this number.
#'                        
#' @param orangeThreshold float between 0 and 1. Each repetetion is considered suscessfull
#'                        if the proportion of orange pixels in the SiZer maps is greater or
#'                        equal to this number.
#'                        
#'                        Note that the we only count colored pixels. If a pixels haven't
#'                        been colored due to being outside of the ESS, that pixel doesn't
#'                        count for the total of pixels inside the SiZer.
#'                        
#' @param xgrid   how many samples we take from the data. This also determine several plots, including the SiZer X lenght
#' 
#' @param ygrid   how many h we are going to generate. This also determine several plots, including the SiZer Y lenght
#'
#' @param kernel  which kernel do you want
#' 
#' @param myMethod which type of interpolation do you use:
#' 
#' @param variance how to calculate the variances.
#' 
#' @param quantileMethod You can choose which quantile method to use when
#'                       calculating the confident intervals.
#' 
#' @param ESSLimit Effective Sample Space.
#'                 
#' @param bootstrapSample If you choose to use the variance with the bootstrap
#'                        method, you can specify the number of sample use for
#'                        the bootstrap
#'                        
#' @return

runTest <- function(dataPath, testRepetitions, blueThreshold, orangeThreshold,
                    xgrid, ygrid, kernel, myMethod, variance, quantileMethod,
                    ESSLimit, bootstrapSample){

  # We are going to save the results in a CSV file in here:
  righNow           = gsub("-","",gsub(" ","_",gsub(":","",as_datetime(Sys.time()))))
  csvFileNameBlue   = paste(righNow,"_","blueTestResult.csv"  , sep = "")
  csvFileNameOrange = paste(righNow,"_","orangeTestResult.csv", sep = "")
  runFolder         = file.path(paste(getwd(),"/results/", sep = ""))
  csvPathBlue       = file.path(runFolder, csvFileNameBlue)
  csvPathOrange     = file.path(runFolder, csvFileNameOrange)
  
  # ---- Get data from file
  myData = getDataFromFile(dataPath)

  # Normalize the data
  # myData = myData/mean(myData)
  
  # TODO: Bootstrap the data???
  
  # Get total tests
  # This is the total set of test, so around 3. Is not the 1000s number.
  totalTests  = length(testRepetitions)
  allzeros    = rep(0,totalTests) # dummy variable
  
  # Create a list for each of the test
  # Each list contain another list with size the total of repetitions for each test
  # In each of those slots, we save the amount of success of tests
  # So each of these list do have 1000s of slots to put results inside.
  blueTestResults   = rep(list(), totalTests)
  orangeTestResults = rep(list(), totalTests)
  
  for (i in 1:totalTests){
    
    blueTestResults[[i]]   = rep(NA, testRepetitions[i])
    orangeTestResults[[i]] = rep(NA, testRepetitions[i])
    
  }
  
  
  # Create an empty dataframe with the results for each test
  # In here we will write the result of the hypothesis testing
  # AFTER all the tests are done.
  # Sample Size / Min / Max / Mean / Median / P75 / P95
  blueTestTable   = data.frame(sampleSize=testRepetitions, min=allzeros,
                               max=allzeros, mean = allzeros, median=allzeros,
                               p75=allzeros, p95=allzeros)
  
  orangeTestTable = data.frame(sampleSize=testRepetitions, min=allzeros,
                               max=allzeros, mean = allzeros, median=allzeros,
                               p75=allzeros, p95=allzeros)
  
  
  # Prepare some feedback for the user
  absoluteRepetitions = sum(testRepetitions)
  repetitionsDone     = 0
  
  # For each test
  for (i in 1:totalTests){
    
    # Repeat the test as many times as specified
    repetitions = testRepetitions[i]
    
    # In here, we accumulate the result of each test    
    # Init to FALSE at the beggining of each test
    passBlueTest    = rep(FALSE,repetitions)   
    passOrangeTest  = rep(FALSE,repetitions)
    
    for (j in 1:repetitions){
      
      
      # Tell the user in which step are you:
      
      print("-----------------------------")
      print("I'm in repetition:")
      print(j)
      print("Of test:")
      print(i)
      print("Completed: (%)")
      print(repetitionsDone*100/absoluteRepetitions)
      repetitionsDone = repetitionsDone + 1
      print("-----------------------------")
      
      # TODO: Bootstrap the sample???
      
      # Run the TTT
      myResult <- ttt(myData, xgrid = xgrid, ygrid = ygrid, kernel = kernel, myMethod = myMethod,
                      ESSLimit = ESSLimit, bootstrapSample = bootstrapSample,
                      savePlots = FALSE, saveCSV = FALSE, silent = TRUE)
      
      #TODO: All warnings are because there are missing points in the
      #      family plots or the phi vector. Due to having NA values
      #      when we do the cramer system divided by 0. Just need
      #      to silent all warnings when everything is complete
      #      and there is no other type of warning
      
      # Count the number of pixels
      myHypothesis <- summaryColors(myResult)
      
      #print(myHypothesis)
      
      # Save the result for later and repeat
      # -- Zero  Sizer
      #    (Save nothing from the zero SiZer)
      
      # -- First Sizer
      # ---- Red
      # ---- Purple
      # ---- Blue
      #    (Save only the blue information from the first SiZer)
      #if(myHypothesis$percentage[6] > blueThreshold ){
        
       # passBlueTest[j] = TRUE
        
      #}
      
      # -- Second Sizer
      # ---- Orange
      #if(myHypothesis$percentage[9] > orangeThreshold ){
        
       # passOrangeTest[j] = TRUE
        
      #}
      
      # ---- Cyan
      # ---- Green
      
      
      # Count pass statistics and write it in dataset
      # Put it into the appropiate blue list
      #totalPassed = sum(passBlueTest)
      #blueTestResults[[i]][j] = totalPassed
      
      # Repeat for other colors
      #totalPassed = sum(passOrangeTest)
      #orangeTestResults[[i]][j] = totalPassed
      
      blueTestResults[[i]][j]   = myHypothesis$percentage[6]
      orangeTestResults[[i]][j] = myHypothesis$percentage[9]
      
    }
    
  }
  
  #summaryBlue   = summary(blueTestResults)
  #summaryOrange = summary(orangeTestResults)
  
  #print(summaryBlue)
  #print(summaryOrange)
  
  print(blueTestResults)
  print(orangeTestResults)
  
  
  # Do the statistics and save it into the test tables
  
  
  write.csv2(blueTestTable,   file = csvPathBlue)
  write.csv2(orangeTestTable, file = csvPathOrange)
  
  #return(c(blueTestTable, orangeTestTable))
  
  return(c(blueTestResults, orangeTestResults))
  
}


# Saved data files
#weibullAditiveDataPath = file.path(paste(getwd(),"/data/weibullAditive.txt", sep = ""))
#logDataPath            = file.path(paste(getwd(),"/data/logNormal.txt",      sep = ""))
ThreeWeibullsgDataPath = file.path(paste(getwd(),"/data/3weibull.txt",       sep = ""))
GranugDataPath         = file.path(paste(getwd(),"/data/granulocytic1.txt",  sep = ""))
AarsetgDataPath        = file.path(paste(getwd(),"/data/aarset.txt",         sep = ""))

ThreeWeibullsgDataPathLittle = file.path(paste(getwd(),"/data/3weibullllittle.txt",    sep = "")) # This one is the same but only 20 datapoints
                                                                                        # Used for testing

dataPath = ThreeWeibullsgDataPathLittle

# Select inputs for the TTT
xgrid = 50
ygrid = 5

# ---- Which kernel do you want
#kernel  = "epanechnikov" 
kernel  = "gaussian"     
#kernel  = "triweight" 

# ---- Smoothing method
#myMethod = "quadratic"
myMethod = "cubic"

# ---- Quantile method 
quantileMethod = "normal"
#quantileMethod = "simultaneous"

# ---- Samples for the bootstrap
bootstrapSample = 100

ESSLimit = 5


# FOR INDIVIDUAL TTT TESTING

#myData = getDataFromFile(ThreeWeibullsgDataPath)
#tttResult <- ttt(myData, xgrid=xgrid, ygrid=ygrid, kernel=kernel, myMethod=myMethod, ESSLimit=ESSLimit, bootstrapSample=bootstrapSample, savePlots = TRUE, saveCSV = FALSE, silent = TRUE)

# FOR HYPOTHESIS TESTING

# Select inputs for the hypothesis testing
# ---- How many times do you want to repeat the tests
testRepetitions = c(10,20,30)                
# ---- The proportion of blue pixels that should be in order to consider a test pass
blueThreshold     = 0.2                      
orangeThreshold   = 0.15   


testResults = runTest(dataPath, testRepetitions, blueThreshold, orangeThreshold,
                      xgrid, ygrid, kernel, myMethod, variance, quantileMethod,
                      ESSLimit, bootstrapSample)

