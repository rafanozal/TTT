



# Add the needed libraries
# ---- Basics
library(ggplot2)       # Basic ggplot2 library
library(ggnewscale)    # Allows for multiple color scales in one plot
library(RColorBrewer)  # Color settins and palettes
library(shadowtext)    # Drop shadows in text for better viewing
library(latex2exp)     # Allows latex expressions in ggplot
# ---- Our own
source("toolsTTT.R",   encoding="utf-8")


# Return the proper Vd and Vg values from a given bias
# For example, for "vd05vg05", the return values would be
# Vd = 0.5 , Vg = 0.5
getBiasValues <- function(biasString){
  
  returnVd = 1
  returnVg = 1
  
  if(biasString == "vd05vg05"){
    returnVd = 0.5
    returnVg = 0.5
  }
  else{
    if(biasString == "vd1vg05"){
      returnVg = 0.5
    }
    else{
      if(biasString == "vd05vg1"){
        returnVd = 0.5
      } 
    }
  }
  
  myReturn = vector("list", length = 2)
  myReturn[[1]] = returnVd
  myReturn[[2]] = returnVg
  
  return(myReturn)
  
  
}

# Return the proper X and Y values from a given bias
# For example, for 1_m2 the return values will be:
# X = 1 , Y = -2
getXYValues <- function(coordinatesString){
  
  x = 0
  y = 0
  
  negativeX = FALSE
  negativeY = FALSE
  
  # Split by "_"
  values  = strsplit(coordinatesString, "_")[[1]]
  xString = values[1]
  yString = values[2]
  
  # For x and y
  # If has m,
  if(grepl("m", xString)){
    #     flag as negative
    negativeX = TRUE
    #     delete the m
    xString   = substring(xString, 2)
    
  }
  if(grepl("m", yString)){
    negativeY = TRUE
    yString   = substring(yString, 2)
    
  }
  
  # Cast the number to integers and multiply by negative if needed
  x = as.integer(xString)
  y = as.integer(yString)
  if(negativeX) x = x * (-1)
  if(negativeY) y = y * (-1)
  
  myReturn = vector("list", length = 2)
  myReturn[[1]] = x
  myReturn[[2]] = y
  
  return(myReturn)
  
}

# From a table of values, get the element I,J where,
# I is the individual (each row) from 1 to N
# J is the Frequency.
#     J = 1 , F = 100
#     J = 2 , F = 1000
#     J = 3 , F = 10000
#
#     In any case J is NOT the index of the table, J only means which one of the
#     frequencies do you want.
getZIJ <- function(tableBase, I, J){
  
  return(tableBase[I,J+1]) # +1 Because all the tables have the same order of columns, Vth, F102, F103, F104, Noise
  
}

# From a table of values, get the Z.. value
getZDotDot <- function(tableBase){
  
  totalRows = nrow(tableBase)
  totalJ    = 3
  
  sumVariable = 0
  
  for (i in 1:totalRows) {
    for(j in 1:totalJ){
      
      sumVariable = sumVariable + getZIJ(tableBase, i, j)
      
    }
    
  }
  
  return(1/totalRows * 1/totalJ * sumVariable)
  
}

# From a table of values, generate all possibles yij
getYIJTable <- function(tableBase){
  
  # Get the dimensions
  totalRows = nrow(tableBase)
  totalJ    = 3
  
  # Get the Zdotdot constant
  Zdotdot = getZDotDot(tableBase)

  # Prepare the dataframe where everything goes
  yDF           =  data.frame(matrix(NA, nrow = totalRows, ncol = totalJ))
  colnames(yDF) = c("F10e2", "F10e3", "F10e4")
  
  # Find all the cells
  for (i in 1:totalRows){
    for(j in 1:totalJ){
      
      yDF[i,j] = getZIJ(tableBase, i, j) - Zdotdot
      
    }
    
  }
  
  return(yDF)
  
}


# Backfitting algorithm
#
# totalSteps is how many times do you want to run the loop to test for
# convergence.
#
# tableBase is a table of size N and, at least 5 columns, in this order:
#     Vth, F10e2, F10e3, F10e4, Noise.
#
# h is the value for the kernel that is use in Kh(xi-xr)
#
# kernel  - gaussian
#         - biweight
#         - triweight
#         - epanechnikov
#
# Return the following:
#   -- A vector of size totalSteps with all the calculated M
#   -- A dataframe of size totalSteps x 3 with the B vectors. Each row
#      correspond to each of the steps calculated
backfitting <- function(totalSteps, tableBase, h, kernel = "gaussian",
                        stopThreshold = 0.000001, minSteps = 5){
  
  # Init the basic variables
  mDF = 0
  bDF = 0
  
  # If you have at least one step do something
  if(totalSteps > 0){

    # Get the dimensions and Yij matrix
    totalRows = nrow(tableBase)
    totalJ    = 3
    yDF       = getYIJTable(tableBase)
    
    # Find the averages for each row in the yDF matrix
    # We need this for later in the backfitting to find the M
    yRowAverages = rowMeans(yDF, na.rm=TRUE)
    
    # Find the averages for each column in the yDF matrix
    # We need this for later in the backfitting to find the B
    yColumnAverages = as.numeric(colMeans(yDF, na.rm=TRUE))
    
    # alias for the variables, so this is easy to read based on the formulas
    n = totalRows
    
    # Get the X vector (Vth values)
    X = tableBase[,1]
    
    # Prepare the variables where to write the results
    mDF = data.frame(matrix(NA, nrow = totalSteps, ncol = n))
    bDF = data.frame(matrix(NA, nrow = totalSteps, ncol = totalJ))

    # Do the first step manually
    # -- M(xi)
    for(i in 1:n){
      mDF[1,i] = 0
    }
    # -- B(j)
    for(j in 1:totalJ){
      
      bDF[1,j] = yColumnAverages[j]
      
    }
    
    # In order to find the the a's values later we need these two variables.
    #
    # -- First we need the kernel cube. In this case, the cube is only a slide,
    #    but the function is called kernelcube anyway
    myKernelSlide    = generateKernelHCube(X,X,h,kernel)
    # -- Now we need the polynomial cube. The polinomial goes from 0 to 6, the
    #    find a's function calculate all the a's from 0 to 6 also. Although we
    #    will doing a bit of overwork, we don't care.
    myPolinomialCube = generatePolinomialCube(X, X)
    
    # Do the rest of the steps
    # We go from step 2 to whatever is the final step
    # We also might stop if the difference from previous is too low
    stopThresholdReached = FALSE
    for(s in 2:totalSteps){
      
      # If we haven't find the stop condition yet, keep doing stuff
      if(stopThresholdReached == FALSE){
        
        # mAverage is the average of all the m of the previous step
        mAverage = mean(as.numeric(mDF[s-1,]))
        
        # First calculate the M(Xi)
        # -- M
        #     For each of the values in the table base with the Vth
        for(i in 1:n){
          
          # A's vector, this is common for all the Rs
          # hIndex = 1 because we only have 1 h, which is the one given in the function
          # kernelCube is only kernelSlide because we only have 1 h
          aVector = findLittleAvalues(i,1,myKernelSlide,myPolinomialCube)
          
          # Accumulate the sum here
          sumVariable = 0
          
          # For each of the r, also from 1 to n
          for(r in 1:n){
            
            # Complete a's fraction including (xi-xr)
            # R is a horrible language and the indexes start at 1 not 0.
            # So in real math notation, we get a +1 in a bunch of indexes here:
            # a0 = aVector[1]
            # a1 = aVector[2], and so on
            # myPolinomialCube[[2]] means that we are using the (xi-xr)^1
            # myPolinomialCube[[2]] means that we are using the (xi-xr)^1
            aFraction = (aVector[3] - aVector[2] * myPolinomialCube[[2]][i,r])/(aVector[3]*aVector[1] - aVector[2]^2) # (xi-xr)^1 always, the exponent doesn't change.
            
            # Kh (xi-xr)
            kernelValue = kernelHFunction(tableBase[i,1],h,kernel)
            
            # Finally, all together 
            sumVariable = sumVariable + ( aFraction * kernelValue * (yRowAverages[r] + mAverage) )
            
          }
          
          # Write the result in the Xi 
          mDF[s,i] = sumVariable
        }
        
        # -- B
        #    For each of the j variables which are the frequencies
        #    In this case is only 1,2,3
        for(j in 1:totalJ){
          
          bDF[s,j] = yColumnAverages[j] - mean(as.numeric(mDF[s,]))
          
        }
        
        # We might have reach the given stop threshold, check that out
        stepsDifference   = mDF[s,] - mDF[(s-1),]
        averageDifference = mean(as.numeric(stepsDifference))
        print(averageDifference)
        if(averageDifference <= stopThreshold){
          # If we also performed enough steps, then stop
          if(minSteps < s){
            stopThresholdReached = TRUE  
          }
          
        } 
        
      }
      
      # If we find the stop condition, the rest of the DF remains NA
      # First NA indicates which step was reached

    }
    
     
  }
  # Otherwise, return an error to the user and complain that you need more steps
  else{
    
    print("I need to do at least one step!")
    
  }
  
  # Return
  myReturn = vector("list", length = 2)
  myReturn[[1]] = mDF
  myReturn[[2]] = bDF
  return (myReturn)
  
  
}
