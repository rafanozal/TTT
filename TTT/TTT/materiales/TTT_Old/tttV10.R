# Add the needed libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(lubridate)
library(latex2exp)
library(ggpubr)

source("toolsV9.R",    encoding="utf-8")


#' For a given array of data, creates several plots describing ....
#' TODO: Write a comprenhsive summary of TTT and what each plot do
#'
#' @param myData  array with the data. Doesn't need to be sorted.
#' 
#' @param xgrid   how many samples we take from the data. This also determine several plots, including the SiZer X lenght
#' 
#'                Default is 401. Beware that this could take some time to calculate.
#' 
#' @param ygrid   how many h we are going to generate. This also determine several plots, including the SiZer Y lenght
#' 
#'                Default is 11.
#' 
#' @param hMin    the minimum h you want to try.
#'  
#'                Default is 1/(ygrid-1).
#'                
#'                If you set a minimum bigger than the maximum, or smaller than
#'                0, the default will be used instead
#'                
#' @param hMax    the maximum h you want to try.
#' 
#'                Default is 1.
#'                
#'                If you set a maximum smaller than the minimum, or bigger than
#'                1, the default will be used instead
#'      
#'      
#' @param kernel  which kernel do you want:
#' 
#                 "gaussian" (DEFAULT)
#'                "epanechnikov" 
#'                "biweight"
#'                "triweight"
#' 
#' @param myMethod which type of interpolation do you use:
#' 
#'                "quadratic" (DEFAULT)
#'                "cubic"
#' 
#' @param variance how to calculate the variances (CPU intensive)
#' 
#'                "bootstrap" (DEFAULT)
#'                    - Does a bootstrap of 100 iteration over the data TODO: Explain better
#'                    
#'                "moments"
#'                    - Calculate the variance using classical central moments.
#'                      When the number of samples is big the factorial number does whacky things and you can get weird results
#'                      
#'                "beta"
#'                    - Calculate the variance using the beta distribution
#'                      TODO: Explain double beta using Dirichlet distribution
#' 
#' @param quantileMethod You can choose which quantile method to use when
#'                       calculating the confident intervals.
#'                       
#'                       "normal" (DEFAULT)
#'                       
#'                           - Classic Z-score based on a normal distribution.
#'                             For this method you can specify the alpha
#'                             parameter. For example, for alpha 0.05, you will
#'                             get a Z of 1.96
#'                            
#'                       "simultaneous"
#'                       
#'                          - It uses a Z-score different for each pixel,
#'                            depending on the h value for that pixel.
#'                            
#'                            Theta(delta) =
#'                            
#'                            2 null { (delta sqrt(3log(xgrid))) / 2h } - 1
#'                            
#'                            q = null^-1 ((1-alpha/2)^(1/Theta(delta)xgrid))
#'                            
#' 
#' @param alpha Level of significance for the confident intervals.
#' 
#'              Default is 0.05
#'  
#'              See the quantileMethod for more info
#' 
#' @param ESSLimit Effective Sample Space. (Default = 5)
#' 
#'                 How many numbers do you need to have around to be a valid result. 
#'                 
#' @param bootstrapSample If you choose to use the variance with the bootstrap
#'                        method, you can specify the number of sample use for
#'                        the bootstrap
#'                        
#'                 Default is 100
#'                 
#' @param saveCSV  Save all numbers used to generate all the plots into a CSV file (Default = TRUE)
#' 
#' @param saveLog  Save every single calculation into a TXT file. File grow exponential, and is around 20MB for 400 x 10 run. (Default = FALSE)
#'                 Use this only for debuggin.
#' 
#' @param blackAndWhite Save your plots with a black and white theme, so the journal don't whine about printing technicolors (Default = FALSE)
#'
#' @return Several plots will appear into the result folder. Which is timestamp to the moment you run everything.
#'         The function itself, return the raw data for all the plots in a dataframe with this form:
#'         
#'         p0 - A value between 0 and 1
#'         h  - The softener used for the kernel
#'         
#'         This two numbers correspond which each pixel in the SiZer plots. So each row of the dataframe represent the information for those pixels
#'         
#'         phiZero
#'         phiOne
#'         phiTwo
#'         
#'         The values for the function, first derivative, and second derivative.
#'         
#'         zeroVariance
#'         firstVariance
#'         secondVariance
#'         
#'         The variance for that given derivatives
#'         
#'         ESS - A value between 0 and infinity.
#'         
#'         The Effective Sample Space value for that pixel.
#'         
#'         LeftIntervalZero
#'         LeftIntervalFirst
#'         LeftIntervalSecond
#'         RightIntervalZero
#'         RightIntervalFirst
#'         RightIntervalSecond
#'         
#'         These are the limits on the left and the right, for the confident interval in each derivative.
#'         
#'         ZeroInZero
#'         ZeroInFirst        (average inside)
#'         ZeroInSecond
#'         
#'         ZeroSmallerZero
#'         ZeroSmallerFirst   (average under)
#'         ZeroSmallerSecond
#'         
#'         ZeroBiggerZero
#'         ZeroBiggerFirst    (average above)
#'         ZeroBiggerSecond
#'         
#'         This are all boolean values and tell you whether the average is
#'         inside the confident interval, under the lowest limit, or above
#'         the upper limit. For each of the derivatives.
#'         
#'         ColorCodeZero
#'         ColorCodeFirst
#'         ColorCodeSecond
#'         
#'         Which color correspond in the SiZer map
#'         
#'         distanceZero
#'         distanceFirst
#'         distanceSecond
#'         
#'         In the continuos SiZer maps, this tells you how many sigmas away is
#'         the derivative which respect the average
#'         
#'         PhiVector
#'         
#'         The phi vector value that is assigned to that p0 value. This column
#'         is redundant and it repeat itself each time the same p0 appear in a
#'         row.
#'         
#'         
#' @export
#'
#' @examples
#' 
#' This should take about 40 seconds if you run it with a CPU made from a toaster:
#' 
#' xgrid = 50
#' ygrid = 10
#' myData  = getRandomData(50,"gamma", 1/5, 5)
#'
#' ttt(myData, xgrid, ygrid, kernel = "epanechnikov", myMethod = "quadratic", variance = "bootstrap",  ESSLimit = 5, saveCSV = TRUE, saveLog = FALSE)
#' 
#' 
ttt <- function(myData, xgrid = 401, ygrid = 11, hMin = 0, hMax = 0,
                kernel = "gaussian", myMethod = "quadratic",
                variance = "bootstrap",  quantileMethod = "normal", 
                alpha = 0.05, ESSLimit = 5, bootstrapSample = 100,
                saveCSV = TRUE, saveLog = FALSE, blackAndWhite = FALSE){
  
  # Inits
  n       = length(myData)
  h       = ygrid
  myData  = sort(myData)/mean(myData)
  
  
  # Timing stuff
  previosTime    = Sys.time()
  currentTime    = Sys.time()
  previosTime2   = Sys.time()
  currentTime2   = Sys.time()
  righNow        = gsub("-","",gsub(" ","_",gsub(":","",as_datetime(Sys.time()))))
  varianceTime   = 0
  covarianceTime = 0
  thetaTime      = 0
  phiTime        = 0
  logTime        = 0
  poliCubeTime   = 0
  kernelTime     = 0
  
  # Set folders for saving plots and plot titles
  # -----------------------------------------------------------------------------
  
  runFolder   = file.path(paste(getwd(),"/results/", righNow, "/", sep = ""))
  
  dataDensityPlotPath         <- file.path(runFolder, "dataPlot.png")
  zeroDerivativePlotPath      <- file.path(runFolder, "zeroPlot.png")
  firstDerivativePlotPath     <- file.path(runFolder, "firstPlot.png")
  secondDerivativePlotPath    <- file.path(runFolder, "secondPlot.png")
  thirdDerivativePlotPath     <- file.path(runFolder, "thirdPlot.png")
  siZerZeroPlotPath           <- file.path(runFolder, "SiZerPlotZeroDerivative.png")
  siZerOnePlotPath            <- file.path(runFolder, "SiZerPlotFirstDerivative.png")
  siZerTwoPlotPath            <- file.path(runFolder, "SiZerPlotSecondDerivative.png")
  distanceSiZerZeroPlotPath   <- file.path(runFolder, "distanceSiZerZeroPlot.png")
  distanceSiZerFirstPlotPath  <- file.path(runFolder, "distanceSiZerFirstPlot.png")
  distanceSiZerSecondPlotPath <- file.path(runFolder, "distanceSiZerSecondPlot.png")
  phiNplotPath                <- file.path(runFolder, "phiNPlot.png")
  cumPlotPath                 <- file.path(runFolder, "cumPlot.png")
  allFamiliesPlotPath         <- file.path(runFolder, "allFamiliesPlot.png")
  allSizersPlotPath           <- file.path(runFolder, "allSizersPlot.png")
  sampleSpacePlotPath         <- file.path(runFolder, "samplespacePlot.png")
  
  csvPath                     <- file.path(runFolder, "resultData.csv")
  logPath                     <- file.path(runFolder, "log.txt")
  
  # Prepare folders and files basic structure
  # -----------------------------------------------------------------------------
  dir.create(runFolder)
  
  # Init some variables
  # -----------------------------------------------------------------------------
  
  # Init the kernels
  kernelZeroBarCube      = 0
  kernelZeroBarBarCube   = 0
  kernelFirstBarCube     = 0
  kernelFirstBarBarCube  = 0
  kernelSecondBarCube    = 0
  kernelSecondBarBarCube = 0
  
  # Init the matrices where we are going to store the A's and a's values
  A0Matrix = matrix(0, nrow = xgrid, ncol = ygrid)
  A1Matrix = matrix(0, nrow = xgrid, ncol = ygrid)
  A2Matrix = matrix(0, nrow = xgrid, ncol = ygrid)
  A3Matrix = matrix(0, nrow = xgrid, ncol = ygrid)
  a0Matrix = matrix(0, nrow = xgrid, ncol = ygrid)
  a1Matrix = matrix(0, nrow = xgrid, ncol = ygrid)
  a2Matrix = matrix(0, nrow = xgrid, ncol = ygrid)
  a3Matrix = matrix(0, nrow = xgrid, ncol = ygrid)
  a4Matrix = matrix(0, nrow = xgrid, ncol = ygrid)
  a5Matrix = matrix(0, nrow = xgrid, ncol = ygrid)
  a6Matrix = matrix(0, nrow = xgrid, ncol = ygrid)
  
  # Init the plots
  sizerZero   = 0
  sizerFirst  = 0
  sizerSecond = 0
  
  # Find out all Pis and P0s
  pisPoints = seq(1,n)/n
  p0sPoints = seq(1,xgrid)/xgrid
  delta = p0sPoints[2]-p0sPoints[1]
  
  # Set out the data that we are going to evaluate and use to make the plots
  minimumXi = min(myData)        # Minimum point to represent in the X axys
  maximumXi = max(myData)        # Maximum point to represent in the X axys
  minimumPi = min(pisPoints)     # Minimum point to represent in the X axys (0 to 1)
  maximumPi = max(pisPoints)     # Maximum point to represent in the X axys (0 to 1)
  minimumH  = 2/(xgrid-1)        # Minimum point to represent in the Y axys 
  maximumH  = 1                  # Maximum point to represent in the Y axys 
  
  if(hMin != 0 && hMin < hMax){ # Correct the defaults hs if needed
    
    minimumH = hMin
    
  }
  
  if(hMax != 0 && hMax > hMin){ 
    
    maximumH = hMax
    
  }
  
  # Set out the vectors with p0 frecuencies and h's to evaluate
  #xiValues  = seq(minimumXi,maximumXi,length = xgrid)                # All the valid xi values to evaluate
  #p0Values  = seq(minimumPi,maximumPi,length = xgrid)                # All the valid p0 values to evaluate (frecuencies)
  
  hValues   = 10^seq(log10(minimumH),log10(maximumH),length = ygrid) # All the valid h values to evaluate
  
  maxHeight = hValues[h] - hValues[h-1] # ggplot will draw only lines in the h's, this serves
                                        # to fill the gap in between if you don't do the log(h) trick
  
  # Generate constant data and time the algorithms
  
  print("Please wait, I'm doing hardcore math")
  
  print("Generating weights matrix...")
  wMatrix          = generateWeightMatrix(n) # your weight matrix, constant for the whole problem
  phiVector        = as.vector(wMatrix %*% myData)
  phiVector[1]     = 0
  phiVector[n]     = 1
  phiVector[phiVector<0]=0
  phiVector[phiVector>1]=1
  
  print("Finding polinomial...")
  previosTime      = Sys.time()
  poliCube         = generatePolinomialCube(pisPoints, p0sPoints)
  currentTime      = Sys.time()
  poliCubeTime     = currentTime - previosTime
  
  print("Calculating kernels...")
  previosTime      = Sys.time()
  kernelCube       = generateKernelHCube(pisPoints,p0sPoints,hValues,kernel)
  currentTime      = Sys.time()
  kernelTime       = currentTime - previosTime
  
  print("Doing variances...")
  previosTime      = Sys.time()
  varianceVector   = generateVarianceXVectorBOOTSRAP(myData, bootstrapSample)
  currentTime      = Sys.time()
  varianceTime     = currentTime - previosTime
  
  print("Making covariances...")
  previosTime      = Sys.time()
  covarianceMatrix = generateCovarianceXMatrixBOOTSTRAP(myData, bootstrapSample)
  currentTime      = Sys.time()
  covarianceTime   = currentTime - previosTime
  
  print("Making sigma matrix")
  sigmaMatrix      = generateSigmaMatrix(varianceVector, covarianceMatrix)

  # Prepare the variable where we are going to collect all data for the plot
  x                     = p0sPoints
  y                     = hValues
  CIData                = expand.grid(p0=x, h=y)
  CIData$phiZero        = 0
  CIData$phiOne         = 0
  CIData$phiTwo         = 0
  CIData$zeroVariance   = 0
  CIData$firstVariance  = 0
  CIData$secondVariance = 0
  CIData$ESS            = 0
  
  # Init other variables
  currentZeroVariance   = 0
  currentFirstVariance  = 0
  currentSecondVariance = 0
  currentThetas         = 0
  
  # Prepare a dataframe for the base data, for later ggplot
  BaseData = data.frame(xis=myData, pis=pisPoints, phi=phiVector, var=varianceVector)
  
  
  # First, find out the thetas
  # -----------------------------------------------------------------------------
  previosTime2 = Sys.time()
  
  print("Calculating phi_h(p0)...")
  value = 0
  for (i in 1:xgrid) { # For each of the p0 values, 401 by default
    
    # Feedback to the user about how much is done and how much time is left
    print(paste(round(i/xgrid*100,2)," %"))
    previosTime = Sys.time()
    
    for (j in 1:ygrid) { # For each of the hs, 11 by default
      
      # Find the right row
      myRowIndex = (j-1)*(xgrid)+i # This do first h's, and then the p0
      
      # Find the right p0 and h
      # currentXi = xiValues[i]
      currentP0 = p0sPoints[i]
      currentH  = hValues[j]
      
      # Find the effective sample space
      CIData$ESS[myRowIndex] = sum( kernelHFunction( currentP0 - pisPoints , currentH , kernel ) )/ kernelHFunction(0,currentH,kernel)
      
      # Find the A's and a's to generate the matrices for the theta function
      currentBigAs    = findBigAvalues(i,j,phiVector,kernelCube,poliCube)
      currentLittleAs = findLittleAvalues(i,j,kernelCube,poliCube)
      
      # Find all thetas
      if(myMethod == "quadratic"){
        
        currentThetas = findThetas(currentBigAs, currentLittleAs)
        
      }
      else{
        
        # currentThetas = findThetasCubic(currentLittleAs, kernelCube, poliCube, j, i, phiVector)
        currentThetas = findThetascubic_Otro(currentBigAs, currentLittleAs)
      }
      
      
      
      # Put them into the data to plot
      CIData$phiZero[myRowIndex]   = currentThetas[1]
      CIData$phiOne[myRowIndex]    = currentThetas[2]
      CIData$phiTwo[myRowIndex]    = 2 * currentThetas[3]
      
      # Save the A's for later
      A0Matrix[i,j] = currentBigAs[1]
      A1Matrix[i,j] = currentBigAs[2]
      A2Matrix[i,j] = currentBigAs[3]
      A3Matrix[i,j] = currentBigAs[4]
      
      # Save the a's for later
      a0Matrix[i,j] = currentLittleAs[1]
      a1Matrix[i,j] = currentLittleAs[2]
      a2Matrix[i,j] = currentLittleAs[3]
      a3Matrix[i,j] = currentLittleAs[4]
      a4Matrix[i,j] = currentLittleAs[5]
      a5Matrix[i,j] = currentLittleAs[6]
      a6Matrix[i,j] = currentLittleAs[7]
    }
    
    # Find out timing statistics
    currentTime = Sys.time()
    deltaTime = currentTime - previosTime
    deltaTime = round(deltaTime,2)
    
    # Show time to go
    print((xgrid-i) * deltaTime)
    
  }
  
  currentTime2 = Sys.time()
  thetaTime    = currentTime2 - previosTime2
  
  # Generate the kernelBar and kernelBarBarCubes
  # -----------------------------------------------------------
  # You will do this using all the little a's information that you calculated before
  # (so you don't have to find evertything again)
  
  
  # Note that we call kernel Bar to everything, whether is cubic or quadratic.
  # In the paper, the kernel tilde is the kernel cube to make notation more understandable.
  
  
  print("Calculating new kernels...")
  if(myMethod == "quadratic"){
  
    kernelZeroBarCube      = generateKernelBarHCube(kernelCube, poliCube, a0Matrix, a1Matrix, a2Matrix, a3Matrix, a4Matrix, 0)
    
    kernelFirstBarCube     = generateKernelBarHCube(kernelCube, poliCube, a0Matrix, a1Matrix, a2Matrix, a3Matrix, a4Matrix, 1)
    
    kernelSecondBarCube    = generateKernelBarHCube(kernelCube, poliCube, a0Matrix, a1Matrix, a2Matrix, a3Matrix, a4Matrix, 2)
    
    
  }
  else{
    
    kernelZeroBarCube     = generateKernelBarHCubeCubic(kernelCube, poliCube, a0Matrix, a1Matrix, a2Matrix, a3Matrix, a4Matrix, a5Matrix, a6Matrix, 0)
    kernelFirstBarCube    = generateKernelBarHCubeCubic(kernelCube, poliCube, a0Matrix, a1Matrix, a2Matrix, a3Matrix, a4Matrix, a5Matrix, a6Matrix, 1)
    kernelSecondBarCube   = generateKernelBarHCubeCubic(kernelCube, poliCube, a0Matrix, a1Matrix, a2Matrix, a3Matrix, a4Matrix, a5Matrix, a6Matrix, 2)

  }
  
 
  
  # Now, find out the variance
  # -----------------------------------------------------------
  
  previosTime2 = Sys.time()
  
  print("Calculating variance of phi_h(p0)")
  value = 0
  for (i in 1:xgrid) { # For each of the p0 values, 401 by default
    
    # Feedback to the user about how much is done and how much time is left
    print(paste(i/xgrid*100," %"))
    previosTime = Sys.time()
    
    for (j in 1:ygrid) { # For each of the hs, 11 by default
      
      # Find the right row
      myRowIndex = (j-1)*(xgrid)+i # This do first h's, and then the p0

      if(myMethod == "quadratic"){

        currentZeroVariance   = findVariancePhiWithKernel(j,i, kernelZeroBarCube,   sigmaMatrix, wMatrix, 0)
        currentFirstVariance  = findVariancePhiWithKernel(j,i, kernelFirstBarCube,  sigmaMatrix, wMatrix, 1)
        currentSecondVariance = findVariancePhiWithKernel(j,i, kernelSecondBarCube, sigmaMatrix, wMatrix, 2)

      }
      
      else{
        
        currentZeroVariance   = findVariancePhiWithKernelCubic(j,i, kernelZeroBarCube,   sigmaMatrix, wMatrix, 0)
        currentFirstVariance  = findVariancePhiWithKernelCubic(j,i, kernelFirstBarCube,  sigmaMatrix, wMatrix, 1)
        currentSecondVariance = findVariancePhiWithKernelCubic(j,i, kernelSecondBarCube, sigmaMatrix, wMatrix, 2)
          
          
      }
      
      
      # Put the results into the dataframe
      CIData$zeroVariance[myRowIndex]   = currentZeroVariance
      CIData$firstVariance[myRowIndex]  = currentFirstVariance
      CIData$secondVariance[myRowIndex] = currentSecondVariance
      
    }
    
    # Find out timing statistics
    currentTime = Sys.time()
    deltaTime = currentTime - previosTime
    
    # Show time to go
    print((xgrid-i) * deltaTime)
    
  }
  
  currentTime2 = Sys.time()
  thetaTime    = currentTime2 - previosTime2
  
  
  
  # Find all the information about confident intervals
  # and set the colors according
  # -----------------------------------------------------------
  #
  
  if(quantileMethod == "simultaneous"){
    
    previosTime2 = Sys.time()
    
    print("Calculating quantiles and confident intervals")
    value = 0
    for (i in 1:xgrid) { # For each of the p0 values, 401 by default
      
      # Feedback to the user about how much is done and how much time is left
      print(paste(i/xgrid*100," %"))
      previosTime = Sys.time()
      
      for (j in 1:ygrid) { # For each of the hs, 11 by default
        
        # Find the right row
        myRowIndex = (j-1)*(xgrid)+i # This do first h's, and then the p0
        
        # Find out the z-score for that particular combination of p0 and h
        thetaDelta = 2*pnorm(delta*sqrt(3*log10(xgrid)/(2*hValues[j]))) - 1
        
        q = qnorm((1-alpha/2)^(1/(thetaDelta*xgrid)))
        
        myZScore = q
       
        # Fill the confident interval data one by one
        CIData$LeftIntervalZero[myRowIndex]   = CIData$phiZero[myRowIndex] - myZScore*sqrt(CIData$zeroVariance[myRowIndex])
        CIData$RightIntervalZero[myRowIndex]  = CIData$phiZero[myRowIndex] + myZScore*sqrt(CIData$zeroVariance[myRowIndex])
        CIData$LeftIntervalZero[myRowIndex]   = CIData$LeftIntervalZero[myRowIndex]/CIData$p0[myRowIndex] - p0sPoints
        CIData$RightIntervalZero[myRowIndex]  = CIData$RightIntervalZero[myRowIndex]/CIData$p0[myRowIndex] - p0sPoints
        
        CIData$LeftIntervalFirst[myRowIndex]  = CIData$phiOne[myRowIndex] - myZScore*sqrt(CIData$firstVariance[myRowIndex])
        CIData$RightIntervalFirst[myRowIndex] = CIData$phiOne[myRowIndex] + myZScore*sqrt(CIData$firstVariance[myRowIndex])
        
        CIData$LeftIntervalSecond[myRowIndex]  = CIData$phiTwo[myRowIndex] - myZScore*sqrt(CIData$secondVariance[myRowIndex])
        CIData$RightIntervalSecond[myRowIndex] = CIData$phiTwo[myRowIndex] + myZScore*sqrt(CIData$secondVariance[myRowIndex])

         
      }
      
      # Find out timing statistics
      currentTime = Sys.time()
      deltaTime = currentTime - previosTime
      
      # Show time to go
      print((xgrid-i) * deltaTime)
      
    }
    
    currentTime2 = Sys.time()
    thetaTime    = currentTime2 - previosTime2

    
  }
  
  else{ # Default for normal, or if you don't know how to write simultaneous
  
    myZScore = qnorm(1 - alpha/2)
      
    CIData$LeftIntervalZero   = CIData$phiZero - myZScore*sqrt(CIData$zeroVariance) 
    CIData$RightIntervalZero  = CIData$phiZero + myZScore*sqrt(CIData$zeroVariance)
    #CIData$LeftIntervalZero   = CIData$LeftIntervalZero/CIData$p0 
    #CIData$RightIntervalZero  = CIData$RightIntervalZero/CIData$p0 
    
    CIData$LeftIntervalZero   = CIData$LeftIntervalZero - CIData$p0 
    CIData$RightIntervalZero  = CIData$RightIntervalZero - CIData$p0
    

    CIData$LeftIntervalFirst  = CIData$phiOne - myZScore*sqrt(CIData$firstVariance)
    CIData$RightIntervalFirst = CIData$phiOne + myZScore*sqrt(CIData$firstVariance)

    CIData$LeftIntervalSecond  = CIData$phiTwo - myZScore*sqrt(CIData$secondVariance)
    CIData$RightIntervalSecond = CIData$phiTwo + myZScore*sqrt(CIData$secondVariance)
    

    
  }
  
  # Now you have all the CIData that you need, you can find out wheter the 0 is inside or not
  # for all of them at the same time
  
  # Check if 0 is above, under, or inside the interval of the zero derivative
  CIData$ZeroInZero         = 0 > CIData$LeftIntervalZero & 0 < CIData$RightIntervalZero
  CIData$ZeroSmallerZero    = 0 < CIData$LeftIntervalZero
  CIData$ZeroBiggerZero     = 0 > CIData$RightIntervalZero
  
  # Check if 0 is above, under, or inside the interval of the first derivative
  CIData$ZeroInFirst        = 0 > CIData$LeftIntervalFirst & 0 < CIData$RightIntervalFirst
  CIData$ZeroSmallerFirst   = 0 < CIData$LeftIntervalFirst
  CIData$ZeroBiggerFirst    = 0 > CIData$RightIntervalFirst
  
  # Check if 0 is above, under, or inside the interval of the second derivative
  CIData$ZeroInSecond        = 0 > CIData$LeftIntervalSecond & 0 < CIData$RightIntervalSecond
  CIData$ZeroSmallerSecond   = 0 < CIData$LeftIntervalSecond
  CIData$ZeroBiggerSecond    = 0 > CIData$RightIntervalSecond

  
  
  # Change the NAs values to FALSE for the intervals
  CIData[c("ZeroInZero",   "ZeroSmallerZero",   "ZeroBiggerZero" )][is.na(CIData[c("ZeroInZero",    "ZeroSmallerZero",   "ZeroBiggerZero")])]   <- FALSE
  CIData[c("ZeroInFirst",  "ZeroSmallerFirst",  "ZeroBiggerFirst" )][is.na(CIData[c("ZeroInFirst",  "ZeroSmallerFirst",  "ZeroBiggerFirst")])]  <- FALSE
  CIData[c("ZeroInSecond", "ZeroSmallerSecond", "ZeroBiggerSecond")][is.na(CIData[c("ZeroInSecond", "ZeroSmallerSecond", "ZeroBiggerSecond")])] <- FALSE
  
  # Get the color coding for the SiZer
  # Init to non valid color
  CIData$ColorCodeZero   = "grey" 
  CIData$ColorCodeFirst  = "grey" 
  CIData$ColorCodeSecond = "grey"
  
  # -- For the zero derivative
  
  if(sum(CIData$ZeroSmallerZero)>0){ # Check that we have stuff to change before changing
    CIData[CIData$ZeroSmallerZero,]$ColorCodeZero       = "yellow"  
  }
  
  if(sum(CIData$ZeroInZero)>0){ 
    CIData[CIData$ZeroInZero,]$ColorCodeZero            = "camel"  
  }
  
  if(sum(CIData$ZeroBiggerZero)>0){ 
    CIData[CIData$ZeroBiggerZero,]$ColorCodeZero        = "olivedrab"  
  }
  
  
  # -- For the first derivative
  
  if(sum(CIData$ZeroSmallerFirst)>0){ # Check that we have stuff to change before changing
    CIData[CIData$ZeroSmallerFirst,]$ColorCodeFirst   = "blue"  
  }
  
  if(sum(CIData$ZeroInFirst)>0){ 
    CIData[CIData$ZeroInFirst,]$ColorCodeFirst        = "purple"  
  }
  
  if(sum(CIData$ZeroBiggerFirst)>0){ 
    CIData[CIData$ZeroBiggerFirst,]$ColorCodeFirst    = "red"  
  }
  
  # -- For the second derivative
  
  if(sum(CIData$ZeroSmallerSecond)>0){ 
    CIData[CIData$ZeroSmallerSecond,]$ColorCodeSecond = "orange"  
  }
  
  if(sum(CIData$ZeroInSecond)>0){ 
    CIData[CIData$ZeroInSecond,]$ColorCodeSecond      = "green"  
  }
  
  if(sum(CIData$ZeroBiggerSecond)>0){ 
    CIData[CIData$ZeroBiggerSecond,]$ColorCodeSecond  = "cyan"   #TODO: Warm colors are not consistance between plots
  }
  
  # Find the variance distance for each derivative
  # How many sigmas away is the 0 from the average
  CIData$distanceZero   <- (0-CIData$phiZero)/sqrt(CIData$zeroVariance) 
  CIData$distanceFirst  <- (0-CIData$phiOne)/sqrt(CIData$firstVariance) 
  CIData$distanceSecond <- (0-CIData$phiTwo)/sqrt(CIData$secondVariance) 
  
  
  # Addjust +3 and -3 distance so you get a softer plot
  
  # -- For the zero derivative
  
  validNumberCondition = !is.na(CIData$distanceZero) & !is.nan(CIData$distanceZero)
  
  if(sum( (validNumberCondition & CIData$distanceZero>3)  > 0 )){
    CIData[validNumberCondition & CIData$distanceZero >  3,]$distanceZero = 3
  }
  
  if(sum( (validNumberCondition & CIData$distanceZero< -3) > 0 )){
    CIData[validNumberCondition & CIData$distanceZero <  -3,]$distanceZero = -3
  }
  
  # Change NA and NAN distances to something really big
  CIData$distanceZero[is.na(CIData$distanceZero)]  = 9999
  CIData$distanceZero[is.nan(CIData$distanceZero)] = 9999
  
  
  # -- For the first derivative
  
  validNumberCondition = !is.na(CIData$distanceFirst) & !is.nan(CIData$distanceFirst)
  
  if(sum( (validNumberCondition & CIData$distanceFirst>3)  > 0 )){
    CIData[validNumberCondition & CIData$distanceFirst >  3,]$distanceFirst =   3
  }
  
  if(sum( (validNumberCondition & CIData$distanceFirst< -3) > 0 )){
    CIData[validNumberCondition & CIData$distanceFirst <  -3,]$distanceFirst = -3
  }
  
  # Change NA and NAN distances to something really big
  CIData$distanceFirst[is.na(CIData$distanceFirst)]  = 9999
  CIData$distanceFirst[is.nan(CIData$distanceFirst)] = 9999
  
  
  # -- For the second derivative
  
  validNumberCondition = !is.na(CIData$distanceSecond) & !is.nan(CIData$distanceSecond)
  
  if(sum( (validNumberCondition & CIData$distanceSecond>3)  > 0 )){
    CIData[validNumberCondition & CIData$distanceSecond >  3,]$distanceSecond =   3
  }
  
  if(sum( (validNumberCondition & CIData$distanceSecond< -3) > 0 )){
    CIData[validNumberCondition & CIData$distanceSecond <  -3,]$distanceSecond = -3
  }
  
  # Change NA and NAN distances to something really big
  
  CIData$distanceSecond[is.na(CIData$distanceSecond)]  = 9999
  CIData$distanceSecond[is.nan(CIData$distanceSecond)] = 9999
  
  
  # Apply the ESS layer
  # Everything bellow the ESS limit is set to grey as there is not enought data
  
  CIData$ColorCodeZero[CIData$ESS < ESSLimit]   = "grey"
  CIData$ColorCodeFirst[CIData$ESS < ESSLimit]  = "grey"
  CIData$ColorCodeSecond[CIData$ESS < ESSLimit] = "grey"
  
  
  # Do the plots
  # -----------------------------------------------------------------------------
  # -- Fills vector
  fillsZero  <-  c("yellow"   = "yellow", 
                   "olivedrab"= "olivedrab",
                   "camel"    = "#C19A6B",
                   "grey"     = "grey")
  
  
  fillsFirst <-  c("red"      = "red", 
                   "blue"     = "blue",
                   "purple"   = "purple",
                   "grey"     = "grey")
  
  fillsSecond <- c("orange"   = "orange", 
                   "cyan"     = "cyan",
                   "green"    = "green",
                   "grey"     = "grey")
  
  
  if(blackAndWhite == TRUE){
    
    fillsZero  <-  c("yellow"       = "black", 
                     "olivedrab"    = "white",
                     "camel"        = "grey50",
                     "grey"         = "grey")
    
    
    fillsFirst <-  c("red"          = "white",
                     "blue"         = "black",
                     "purple"       = "grey50",
                     "grey"         = "grey")
    
    
    fillsSecond <- c("orange"   = "white", 
                     "cyan"     = "black",
                     "green"    = "grey50",
                     "grey"     = "grey")
    
  }
  
  
  # Input data plot
  
  # -- Transform data array into data frame
  myDataFrame          = as.data.frame(as.table(myData))
  myDataFrame$Var1     = NULL
  myDataFrame$Relative = myDataFrame$Freq/maximumXi
  
  
  # Plot with the raw data, using density
  dataPlot <- ggplot(myDataFrame, aes(x=Freq)) +
    # Do a density plot
    geom_density(alpha = 0.7, fill="red")+
    # Set the limits 
    xlim(minimumXi, maximumXi) + 
    labs(title = "Density plot for given data",
         x = "x", y = "Density") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_blank(),
          panel.background = element_blank()) 
  
  
  ggsave(dataDensityPlotPath , width = 8)
  
  
  # Plot with the phi vector.
  # The red line is Y = X. Purple line is Y = 0
  phiPlot <- ggplot(BaseData, aes(x=pis, y=phi)) +
    geom_point()+
    geom_line(linetype = 3, color="darkgrey") +
    geom_abline(intercept = 0, slope = 1, linetype = 1, color="red") + 
    geom_abline(intercept = 0, slope = 0, linetype = 1, color="purple") + 
    geom_abline(intercept = 1, slope = 0, linetype = 1, color = "purple" ) +
    labs(title = "TTT-Plot",
         x = "p", y = "") +
    
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_blank(),
          panel.background = element_blank()) 
  
  ggsave(phiNplotPath , width = 8)
  
  
  # The Theta / Phi_h vector / Family plots
  
  # -- Phi_H Zero

  familyZeroPlot <- ggplot(CIData, aes(x=p0, y=phiZero, group = factor(round(h,2)), colour = factor(round(h,2)))) +
    geom_line(linetype = 1) +
    xlim(0, 1) +
    # geom_abline(intercept = 0, slope = 1, linetype = 2, color="red") + 
    geom_abline(intercept = 0, slope = 0, linetype = 1, color="purple") + 
    geom_abline(intercept = 1, slope = 0, linetype = 1, color="purple") + 
    labs(title = TeX("Family plot: $\\widehat{\\phi_h } = \\widehat{\\theta _0}"),
         
         x = "p", y = TeX("$\\widehat{\\phi_h }"),
         
         color='h') +
    
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_blank(),
          panel.background = element_blank()) 
  
  ggsave(zeroDerivativePlotPath , width = 8)
  
  
  # -- Phi_H One
  
  familyFirstPlot <- ggplot(CIData, aes(x=p0, y=phiOne, group = factor(round(h,2)), colour = factor(round(h,2)))) +
    geom_line(linetype = 1) +
    xlim(0, 1) +
    # geom_abline(intercept = 0, slope = 1, linetype = 2, color="red") + 
    geom_abline(intercept = 0, slope = 0, linetype = 1, color="purple") + 
    geom_abline(intercept = 1, slope = 0, linetype = 1, color="purple") + 
    labs(title = TeX("Family plot: $\\widehat{\\phi'_h } = \\widehat{\\theta _1}"),
         
         x = "p", y = TeX("$\\widehat{\\phi'_h }"),
         
         color='h') +
    
    
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_blank(),
          panel.background = element_blank()) 
  
  
  ggsave(firstDerivativePlotPath , width = 8)
  
  
  # -- Phi_H Two
  familySecondPlot <- ggplot(CIData, aes(x=p0, y=phiTwo, group = factor(round(h,2)), colour = factor(round(h,2)))) +
    geom_line(linetype = 1) +
    xlim(0, 1) +
    # geom_abline(intercept = 0, slope = 1, linetype = 2, color="red") + 
    geom_abline(intercept = 0, slope = 0, linetype = 1, color="purple") + 
    geom_abline(intercept = 1, slope = 0, linetype = 1, color="purple") + 
    labs(title = TeX("Family plot: $\\widehat{\\phi''_h } = 2 \\widehat{\\theta _2}"),
         
         x = "p", y = TeX("$\\widehat{\\phi''_h }"),
         
         color='h') +
    
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_blank(),
          panel.background = element_blank()) 
  
  
  ggsave(secondDerivativePlotPath , width = 8)
  
  # -- All together 
  phisandfamilyPlot <- ggarrange(phiPlot,
                                 familyZeroPlot,
                                 familyFirstPlot,
                                 familySecondPlot,
                                 # labels = c("A", "B", "C","D"),
                                 
                                 ncol = 1, nrow = 4)
  
  ggsave(allFamiliesPlotPath , width = 8)
  
  
  
  # FOR THE ZERO DERIVATIVE:
  
  # -- SiZer 
  sizerZero <- ggplot(CIData, aes(p0, log(h), z= ColorCodeZero)) +
    geom_tile(aes(fill = ColorCodeZero)) +
    scale_fill_manual(values = fillsZero) + 
    theme(legend.position = "none") + 
    labs(title = "SiZer-0",
        #  subtitle = "olive => 1 < [CI] ; yellow => [CI] < 1",
         x = "p", y = "")
  
  
  ggsave(siZerZeroPlotPath , width = 8)
  
  # -- SiZer with distance
  ggplot(CIData, aes(p0, log(h), z= distanceZero)) +
    geom_tile(aes(fill = distanceZero)) +
    
    scale_fill_gradient2(midpoint=0,
                         low="green", mid="steelblue2",high="brown",
                         limits = c(-3,3),
                         space ="Lab" ) + #space is gradient algorithm 
    
    labs(title = bquote("Distance:" ~ y[i] == (0 - phi[0]) / sigma(phi[0])),
         x = "p", y = "")
  ggsave(distanceSiZerZeroPlotPath , width = 8)
  
  
  
  # FOR THE FIRST DERIVATIVE:
  
  # -- SiZer 
  sizerFirst <- ggplot(CIData, aes(p0, log(h), z= ColorCodeFirst)) +
    geom_tile(aes(fill = ColorCodeFirst)) +
    scale_fill_manual(values = fillsFirst) + 
    theme(legend.position = "none") + 
    labs(title = "SiZer-1",
         # subtitle = "red => 0 < [CI] ; blue => [CI] < 0",
         x = "p", y = "")
  
  ggsave(siZerOnePlotPath , width = 8)
  
  # -- SiZer with distance
  ggplot(CIData, aes(p0, log(h), z= distanceFirst)) +
    geom_tile(aes(fill = distanceFirst)) +
    
    scale_fill_gradient2(midpoint=0,
                         low="blue", mid="yellow",high="red",
                         limits = c(-3,3),
                         space ="Lab" ) + #space is gradient algorithm 
    
    labs(title = bquote("Distance:" ~ y[i] == (0 - phi[1]) / sigma(phi[1])),
         x = "p", y = "")
  
  ggsave(distanceSiZerFirstPlotPath , width = 8)
  
  
  # FOR THE SECOND DERIVATIVE:
  
  # -- SiZer 
  sizerSecond <- ggplot(CIData, aes(p0, log(h), z= ColorCodeSecond)) +
    geom_tile(aes(fill = ColorCodeSecond)) +
    scale_fill_manual(values = fillsSecond) + 
    theme(legend.position = "none") + 
    labs(title = "SiZer-2",
         #subtitle = "cyan => 0 < [CI] ; orange => [CI] < 0",
         x = "p", y = "")
  
  ggsave(siZerTwoPlotPath , width = 8)
  
  # -- SiZer with distance
  ggplot(CIData, aes(p0, log(h), z= distanceSecond)) +
    
    geom_tile(aes(fill = distanceSecond)) +
    
    scale_fill_gradient2(midpoint=0,
                         low="turquoise", mid="purple",high="orange",
                         limits = c(-3,3),
                         space ="Lab" ) + #space is gradient algorithm 
    
    labs(title = expression("Distance (0-phi_TWO)/sd"),
         x = "p", y = "")
  
  ggsave(distanceSiZerSecondPlotPath , width = 8)
  
  # -- All sizers with the family plot
  phisandfamilyPlot <- ggarrange(phiPlot, sizerZero,
                                 sizerFirst, sizerSecond,
                                 # labels = c("A", "B", "C","D"),
                                 ncol = 2, nrow = 2)
  annotate_figure(phisandfamilyPlot, top = text_grob("TTT-SiZer"))
  ggsave(allSizersPlotPath , width = 8)
  
  
  # Effective Sample Space
  ggplot(CIData, aes(p0, log(h), z= ESS)) +
    geom_tile(aes(fill = ESS)) +
    
    scale_fill_gradient2(midpoint=5,
                         low="black", mid="grey",high="green",
                         space ="Lab" ) + #space is gradient algorithm 
    
    labs(title = "ESS",
        # subtitle = "Black to grey < 5 = not enought points ; grey to green > 5 = enought points",
         x = "p", y = "")
  ggsave(sampleSpacePlotPath , width = 8)
  
  
  
  
  
  # Write the result of everything into the log file
  # -----------------------------------------------------------------------------
  
  if( saveCSV == TRUE){
    
    write.csv2(CIData, file = csvPath)
    
    
  }
  
  if(saveLog == TRUE){
    
    writeLog(logPath, n, h, kernel, xgrid, ygrid, myData, pisPoints, CIData,
             wMatrix, phiVector, poliCube, kernelCube, varianceVector, covarianceMatrix,
             a0Matrix, a1Matrix, a2Matrix, a4Matrix, kernelBarCube, kernelBarBarCube,
             A0Matrix, A1Matrix, A2Matrix)   
    
  }
  
  
  return (CIData)
  
  
}
