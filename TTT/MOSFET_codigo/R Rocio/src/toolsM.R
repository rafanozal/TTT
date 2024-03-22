# ---- Our own
source("toolsT.R",   encoding="utf-8")

# Backfitting algorithm
#
# totalSteps is how many times do you want to run the loop to test for
# convergence.
#
# tableBase is a table of size N and, at least 5 columns, in this order:
#     Vth, F10e4, F10e3, F10e2, Noise.
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
                        stopThreshold = 0.0001, minSteps = 2){
  
  # Init the basic variables
  mDF = 0
  bDF = 0
  pDF = 0
  y.ij = 0
  ym = 0
  yt = 0
  
  # If you have at least one step do something
  if(totalSteps > 0){
    
    # Get the dimensions and Yij matrix
    n = nrow(tableBase)
    totalJ    = 3
    y.ij = tableBase[,3:5]
    ym  = c(tableBase[,3],tableBase[,4],tableBase[,5])
    
    
    ni.factor = c(rep(1e4,n),rep(1e3,n),rep(1e2,n))
    x.i = rep(tableBase[,2],3)
    indiv.i = rep(tableBase[,1],3)
    
    fit = lme(ym~factor(ni.factor),random=~1|indiv.i)
    pred = predict(fit,asList=T,level=0:1)
    pi.i = pred$predict.indiv.i - pred$predict.fixed
    bj = pred$predict.fixed[(0:2)*n+1]
    
    yt = ym - pred$predict.indiv.i
    
     # Prepare the variables where to write the results
    mDF = data.frame(matrix(NA, ncol = totalSteps, nrow = n))
    bDF = data.frame(matrix(NA, nrow = totalSteps, ncol = totalJ))
    pDF = data.frame(matrix(NA, nrow = n, ncol = totalSteps))
    
    # Do the first step manually
    # -- M(xl)
    for(l in 1:n){
      mDF[l,1] = 0
      pDF[l,1] = pi.i[l]
    }
    # -- B(j)
    for(j in 1:totalJ){
      
      bDF[1,j] = bj[j]
      
    }
    
    
    # In order to find the the a's values later we need these two variables.
    #
    # --- First we need the kernel cube. In this case, the cube is only a slide,
    #    but the function is called kernelcube anyway
    myKernelSlide    = generateKernelH0Cube(tableBase[,2],tableBase[,2],h,kernel)
    # -- Now we need the polynomial cube. The polinomial goes from 0 to 6, the
    #    find a's function calculate all the a's from 0 to 6 also. Although we
    #    will doing a bit of overwork, we don't care.
    myPolinomialCube = generatePolinomialCube1(tableBase[,2], tableBase[,2])
    
    # Do the rest of the steps
    # We go from step 2 to whatever is the final step
    # We also might stop if the difference from previous is too low
    stopThresholdReached = FALSE
     # lo que vamos a suavizar
    for(s in 2:totalSteps){
      
      # If we haven't find the stop condition yet, keep doing stuff
      if(stopThresholdReached == FALSE){
        
        ys = rowMeans(matrix(yt,n,3))
                # First calculate the M(Xl)
        # -- M
        #     For each of the values in the table base with the Vth
        for(l in 1:n){
          
          # A's vector, this is common for all the Rs
          # hIndex = 1 because we only have 1 h, which is the one given in the function
          # kernelCube is only kernelSlide because we only have 1 h
          aVector = findLittleAvaluesb(l,1,myKernelSlide,myPolinomialCube)
          
          # Accumulate the sum here
          sumVariable = 0
          
          # For each of the i, also from 1 to n
          for(i in 1:n){
            
            # Complete a's fraction including (xi-xr)
            # R is a horrible language and the indexes start at 1 not 0.
            # So in real math notation, we get a +1 in a bunch of indexes here:
            # a0 = aVector[1]
            # a1 = aVector[2], and so on
            # myPolinomialCube[[2]] means that we are using the (xi-xr)^1
            # myPolinomialCube[[2]] means that we are using the (xi-xr)^1
            aFraction = (aVector[3] - aVector[2] * myPolinomialCube[[2]][i,l])/(aVector[3]*aVector[1] - aVector[2]^2) # (xi-xr)^1 always, the exponent doesn't change.
            
            # Kh (xi-xr)
           # kernelValue = kernelHFunction((tableBase[r,1]-rpoints[i]),h,kernel)
            kernelValue = myKernelSlide[i,l]
            # Finally, all together 
            sumVariable = sumVariable + ( aFraction * kernelValue * (ys[i]) )
            
          }
          
          # Write the result in the Xi 
          mDF[l,s] = sumVariable
        }
        ml = rep(mDF[,s],3)
        y = ym-ml
        # -- B y p
        #    For each of the j variables which are the frequencies
        #    In this case is only 1,2,3
        fit = lme(y~factor(ni.factor),random=~1|indiv.i)
        pred = predict(fit,asList=T,level=0:1)
        pi.i = pred$predict.indiv.i - pred$predict.fixed
        bj = pred$predict.fixed[(0:2)*n+1]
        
        yt = ym - pred$predict.indiv.i
        for(l in 1:n){
          
          pDF[l,s] = pi.i[l]
        }
        # -- B(j)
        for(j in 1:totalJ){
          
          bDF[s,j] = bj[j]
          
        }
        
        
        # We might have reach the given stop threshold, check that out
        stepsDifference   = sum((mDF[,s] - mDF[,(s-1)])^2)/(sum(mDF[,s-1]^2)+0.0001)
       
        print(stepsDifference)
        if(stepsDifference <= stopThreshold){
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
  myReturn = vector("list", length = 3)
  myReturn[[1]] = mDF
  myReturn[[2]] = bDF
  myReturn[[3]] = pDF
  return (myReturn)
  
  
}
backfit.sizer <- function(totalSteps, tableBase, h, kernel = "gaussian",
                        stopThreshold = 0.0001, minSteps = 2){
  
  # Init the basic variables
  mDF = 0
  bDF = 0
  pDF = 0
  y.ij = 0
  ym = 0
  yt = 0
  dmDF=0
  n = nrow(tableBase)
  
  # If you have at least one step do something
  if(totalSteps > 0){
    
    # Get the dimensions and Yij matrix
    
    totalJ    = 3
    y.ij = tableBase[,3:5]
    ym  = c(tableBase[,3],tableBase[,4],tableBase[,5])
    
    
    ni.factor = c(rep(1e4,n),rep(1e3,n),rep(1e2,n))
    x.i = rep(tableBase[,2],3)
    indiv.i = rep(tableBase[,1],3)
    
    fit = lme(ym~factor(ni.factor),random=~1|indiv.i)
    pred = predict(fit,asList=T,level=0:1)
    pi.i = pred$predict.indiv.i - pred$predict.fixed
    bj = pred$predict.fixed[(0:2)*n+1]
    
    yt = ym - pred$predict.indiv.i
    
    # Prepare the variables where to write the results
    mDF = data.frame(matrix(NA, ncol = totalSteps, nrow = n))
    bDF = data.frame(matrix(NA, nrow = totalSteps, ncol = totalJ))
    pDF = data.frame(matrix(NA, nrow = n, ncol = totalSteps))
    
    # Do the first step manually
    # -- M(xl)
    for(l in 1:n){
      mDF[l,1] = 0
      pDF[l,1] = pi.i[l]
    }
    # -- B(j)
    for(j in 1:totalJ){
      
      bDF[1,j] = bj[j]
      
    }
    
    
    # In order to find the the a's values later we need these two variables.
    #
    # --- First we need the kernel cube. In this case, the cube is only a slide,
    #    but the function is called kernelcube anyway
    #myKernelSlide    = generateKernelH0Cube(tableBase[,2],xgrid,h,kernel)
    myKernelSlide    = generateKernelH0Cube(tableBase[,2],tableBase[,2],h,kernel)
    # -- Now we need the polynomial cube. The polinomial goes from 0 to 6, the
    #    find a's function calculate all the a's from 0 to 6 also. Although we
    #    will doing a bit of overwork, we don't care.
    #myPolinomialCube = generatePolinomialCube1(tableBase[,2], xgrid)
    myPolinomialCube = generatePolinomialCube1(tableBase[,2], tableBase[,2])
    
    # Do the rest of the steps
    # We go from step 2 to whatever is the final step
    # We also might stop if the difference from previous is too low
    stopThresholdReached = FALSE
    # lo que vamos a suavizar
    for(s in 2:totalSteps){
      
      # If we haven't find the stop condition yet, keep doing stuff
      if(stopThresholdReached == FALSE){
        
        ys = rowMeans(matrix(yt,n,3))
        # First calculate the M(Xl)
        # -- M
        #     For each of the values in the table base with the Vth
        for(l in 1:n){
          
          # A's vector, this is common for all the Rs
          # hIndex = 1 because we only have 1 h, which is the one given in the function
          # kernelCube is only kernelSlide because we only have 1 h
          aVector = findLittleAvaluesb(l,1,myKernelSlide,myPolinomialCube)
          
          # Accumulate the sum here
          sumVariable = 0
          
          # For each of the i, also from 1 to n
          for(i in 1:n){
            
            # Complete a's fraction including (xi-xr)
            # R is a horrible language and the indexes start at 1 not 0.
            # So in real math notation, we get a +1 in a bunch of indexes here:
            # a0 = aVector[1]
            # a1 = aVector[2], and so on
            # myPolinomialCube[[2]] means that we are using the (xi-xr)^1
            # myPolinomialCube[[2]] means that we are using the (xi-xr)^1
            aFraction = (aVector[3] - aVector[2] * myPolinomialCube[[2]][i,l])/(aVector[3]*aVector[1] - aVector[2]^2) # (xi-xr)^1 always, the exponent doesn't change.
            
            # Kh (xi-xr)
            # kernelValue = kernelHFunction((tableBase[r,1]-rpoints[i]),h,kernel)
            kernelValue = myKernelSlide[i,l]
            # Finally, all together 
            sumVariable = sumVariable + ( aFraction * kernelValue * (ys[i]) )
            
          }
          
          # Write the result in the Xi 
          mDF[l,s] = sumVariable
        }
        ml = rep(mDF[,s],3)
        y = ym-ml
        # -- B y p
        #    For each of the j variables which are the frequencies
        #    In this case is only 1,2,3
        fit = lme(y~factor(ni.factor),random=~1|indiv.i)
        pred = predict(fit,asList=T,level=0:1)
        pi.i = pred$predict.indiv.i - pred$predict.fixed
        bj = pred$predict.fixed[(0:2)*n+1]
        
        yt = ym - pred$predict.indiv.i
        for(l in 1:n){
          
          pDF[l,s] = pi.i[l]
        }
        # -- B(j)
        for(j in 1:totalJ){
          
          bDF[s,j] = bj[j]
          
        }
        
        
        # We might have reach the given stop threshold, check that out
        stepsDifference   = sum((mDF[,s] - mDF[,(s-1)])^2)/(sum(mDF[,s-1]^2)+0.0001)
        
        print(stepsDifference)
        if(stepsDifference <= stopThreshold){
          # If we also performed enough steps, then stop
          if(minSteps < s){
            stopThresholdReached = TRUE  
            smax = s
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
  
  for(l in 1:n){
    
    # A's vector, this is common for all the Rs
    # hIndex = 1 because we only have 1 h, which is the one given in the function
    # kernelCube is only kernelSlide because we only have 1 h
    aVector = findLittleAvaluesb(l,1,myKernelSlide,myPolinomialCube)
    
    # Accumulate the sum here
    sumder = 0
    
    # For each of the i, also from 1 to n
    for(i in 1:n){
      
      # Complete a's fraction including (xi-xr)
      # R is a horrible language and the indexes start at 1 not 0.
      # So in real math notation, we get a +1 in a bunch of indexes here:
      # a0 = aVector[1]
      # a1 = aVector[2], and so on
      # myPolinomialCube[[2]] means that we are using the (xi-xr)^1
      # myPolinomialCube[[2]] means that we are using the (xi-xr)^1
      aFraction = (aVector[2] - aVector[1] * myPolinomialCube[[2]][i,l])/(aVector[3]*aVector[1] - aVector[2]^2) # (xi-xr)^1 always, the exponent doesn't change.
      
      # Kh (xi-xr)
      # kernelValue = kernelHFunction((tableBase[r,1]-rpoints[i]),h,kernel)
      kernelValue = myKernelSlide[i,l]
      # Finally, all together 
      sumder = sumder + ( aFraction * kernelValue * (ys[i]) )
      
    }
    
    # Write the result in the Xi 
    dmDF[l] = sumder
  }
  

  myReturn = vector("list", length = 4)
  myReturn[[1]] = dmDF
  myReturn[[2]] = mDF[,smax]
  myReturn[[3]] = bDF[smax,]
  myReturn[[4]] = pDF[,smax]
  return (myReturn)
  
  
}
