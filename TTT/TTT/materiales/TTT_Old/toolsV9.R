#' Generate random data
#'
#' @param n How many numbers do you want.
#' 
#' @param dist "normal"  ------- Normal distribution
#'                  p1 = mean
#'                  p2 = sd
#'
#'             "poisson" ------- Poisson distribution
#'                  p1 = lamda
#'                  p2 = nothing
#'
#'             "weibull" ------- Weibull distribution
#'                  p1 = shape
#'                  p2 = scale
#'                  
#'             "lognormal" ----- Log Normal distribution 
#'                  p1 = meanlog
#'                  p2 = sdlog
#'                  
#'              "logmix"   ----- Mix of two log distribution
#'                               with the same average and sd
#'                  p1 = meanlog
#'                  p2 = sdlog   
#'                  
#'              "exponential" -- Exponential distribution
#'                  p1 = rate
#'                  p2 = nothing
#'                  
#'              "gamma" -------- Gamma distribution
#'                  p1 = shape
#'                  p2 = scale
#'               
#'                  
#' @param p1 first parameter for your selected distribution
#' 
#' @param p2 second parameter for your selected distribution
#' 
#' @param forceOneMean if this is TRUE, the function will only take into
#'                     account the p1 parameter, and based on that value
#'                     and whatever distribution you selected, addjust
#'                     the p2 parameter and return some random data which
#'                     have an average of 1.
#'                     
#'                     TODO: Finish this
#'
#' @return
#' @export
#'
#' @examples
#' 
#' getRandomData(10,"normal", 100, 5)
#' getRandomData(10,"weibull", 3, 1)
#' 
#' 
getRandomData <- function(n, dist, p1, p2, forceOneMean = FALSE){
  
  myFinalData = 0
  
  if(dist=="normal"){
    myFinalData = rnorm(n = n, mean = p1, sd = p2)  
    
  }
  if(dist=="poisson"){
    myFinalData =  rpois(n, lambda=p1)
    
  }
  
  if(dist=="weibull"){
    
    if(forceOneMean == TRUE){
        p2 = gamma(1+1/p1)^{-1}
        myFinalData = rweibull(n, shape = p1, scale = p2)
    }
    else{
        myFinalData = rweibull(n, shape = p1, scale = p2)
        
    }
  }
    
  if(dist=="lognormal"){
    
    if(forceOneMean == TRUE){
      
      p1 = -(p2^2)/2  
      p1 = (-1) * abs(p1) # Force mu to be negative
            
      myFinalData = rlnorm(n,p1,p2)
      
    }
    else{
      
      myFinalData = rlnorm(n,p1,p2)
      
    }
      
    
  }

  if(dist == "exponential"){
      
        myFinalData = rexp(n,rate=p1)
          
  }
      
  if(dist == "gamma"){
       
       if(forceOneMean == TRUE){
      
          p1 = 1/p2
      
       }
    
        myFinalData = rgamma(n,shape=p1,scale=p2) 
        
  }

  return(myFinalData)
  
}


# Generate random data based on an aditivive weibull model  ###
### Esta funci?n hay que reemplazarla
getRandomDataMixedWeibull <- function(n, shapes, scales, mixing){
  
  totalWeibull = length(shapes)
  returnArray  = array(0,n)
  
  # Generate one weibull for each parameter we have
  weibullArray = list()
  
  for(i in 1:totalWeibull){
    
    weibullArray[[i]] <- rweibull(n, shape=shapes[i],scale=scales[i])
    
  }
  
  # Mix them all toguether:
  
  # -- If you want an exponential mix: 
  if(mixing == "exponential"){

    for (i in totalWeibull) {
      
      returnArray = returnArray + shapes[i]*weibullArray[[i]]^shapes[i]
      
    }
    
  }
  
  # -- If you want a minimum mix:
  if(mixing == "minimum"){
    
    for (i in 1:n){
      
      myMinimum = 99999999
      
      for (j in 1:totalWeibull) {
        
        myMinimum = min(myMinimum, weibullArray[[j]][i])
        
      }
      
      returnArray[i] = myMinimum
      
    }

  }
  
  return (returnArray)
  
}



# Read a file with numbers and return them as vector
getDataFromFile <- function(filePath){
  
  df <- read.table(filePath, header = FALSE)

  return (df[[1]])
  
}

# For a given number, return the value of the kernel function, with the specified kernel
#
# SEE:
#
#     TODO: For the cumulative kernel function that finds the integral, check integralKernelFunction()
#
# INPUT:
#
# x       - the data to find the kernel, can be a vector
# kernel  - gaussian
#         - biweight
#         - triweight
#         - epanechnikov
#
# OUTPUT:
#
# a float vector of size length(x) with the result of the kernel

kernelFunction <- function(x, kernel="gaussian"){
  
  myResult = x*0
  
  if(missing(kernel)){
    kernel = "gaussian"  
  }
  
  if(kernel=="gaussian") {
    myResult = exp(-(x^2) / 2) / sqrt(2 * pi)
  }
  
  if(kernel=="biweight"){
    myResult = (abs(x)<=1)*15/16*(1-x^2)^2
  }
  
  if(kernel=="triweight"){
    myResult = (abs(x)<=1)*35/32*(1-x^2)^3
  }
  
  if(kernel=="epanechnikov"){
    myResult = (abs(x)<=1)*((1-x^2)*3/4)
  }
  
  return(myResult)
  
} 


# Correct the kernel of a given data with a given badwith h and k
#
#
# INPUT:
#
# x       - the data to find the kernel, can be a vector
# h       - bandwith
# kernel  - gaussian
#         - biweight
#         - triweight
#         - epanechnikov
#
# OUTPUT:
#
# a float vector of size length(x) with the result of the kernel
kernelHFunction <- function(x,h,kernel){
  
  myCorrection = (1/h) * kernelFunction(x/h,kernel)
  
  return(myCorrection)
  
}


# For a given vector of pi's and h's, generate all possible K_h(pi-p0) combinations.
#
# INPUT:
#
# piVector - The vector with all the frecuencies (can be a vector of size one)
# p0Vector - The vector with only the p0 values (can be a vector of size one)
# hVector  - The vector with only the h's values (can be a vector of size one)
# kernel   - gaussian
#          - biweight
#          - triweight
#          - epanechnikov
#
# OUTPUT:
#
# A list of matrixs with all the K_h(pi-p0) combinations pre-calculated.
# Where the index of the list represent the h value, the row is pi and the column p0
#
# The row is the h, the column is the sum(K_h(pi-p0))

generateKernelHCube <- function(piVector,p0Vector,hVector,kernel){
  
  totalPs  = length(piVector)
  totalP0s = length(p0Vector) 
  totalHs  = length(hVector)
  
  matrixByHIndex <- list()
  
  # The index of the list refers to the h that you want.
  # For example, hCube[1][2,3] gives you the kernelH value of h=[1], p-p0=[2][3]

  for(i in 1:totalHs){
    
    matrixByHIndex[[i]] <- array(0,dim = c(totalPs, totalP0s))
    
    for(j in 1:totalP0s){
      
      for(k in 1:totalPs){
        
        matrixByHIndex[[i]][k,j] <- kernelHFunction(piVector[k]-p0Vector[j],hVector[i],kernel)
        
      }
      
    }  
    
  }
  
  # Create an empty matrix
  return(matrixByHIndex)
  
}


# For a given vector of pi's and generate all possible (pi-p0)^n combinations.
# Where n=0,1,2,3,4,5,6
#
#
# INPUT:
#
# piVector - The vector with all the frecuencies (can be a vector of size one)
# p0Vector - The vector with only the p0 values (can be a vector of size one)
#
# OUTPUT:
#
# A list of matrix with all the(pi-p0)^n combinations pre-calculated.
#
# The index of the list is n, the row of the matrix is pi and the column is p0

generatePolinomialCube <- function(piVector, p0Vector){
  
  totalPs  = length(piVector)
  totalP0s = length(p0Vector) 
  
  matrixByNIndex <- list()
  
  # First, we generate the cube of n x p x p0
  
  for (k in 1:7) { # We only need 0,1,2,3,4,5,6 n's
    
    matrixByNIndex[[k]] <- matrix(0,totalPs,totalP0s)
    
    for(i in 1:totalPs){
      
      for (j in 1:totalP0s) {
        
        matrixByNIndex[[k]][i,j] = (piVector[i] - p0Vector[j])^(k-1)
        
      }
      
    } 
    
  }
  
  return(matrixByNIndex)
  
}


# Finds a single A_r(p0). To use in the next function that finds all A's.
#
# PRE:
# - a kernel cube with all the K_h(pi-p0) values
# - a polinomial cube with all the (pi-p0)^n values
# - a vector with all the phi(xis) values
#
# INPUT:
#
# p0           - An index number refering to the p0 value
# h            - An index number refering to the bandwith
# r            - TODO: Name? index with the ^r value. 1=0, 2=1, ... 5=4
# phiVector    - The vector with the phi values
# kernelCube   - The cube with all the Kh(pi-p0) values
# poliCube     - The cube with all the (pi-p0)^n values
#
# OUTPUT:
#
# A single float with the Ar(p0) 
findBigA <- function(p0,h,r,phiVector,kernelCube,poliCube){
  
  myAr = sum(phiVector*poliCube[[r]][,p0]*kernelCube[[h]][,p0])
  
  return(myAr)
  
}



# Find As (uppercase)
#
# For a given set of data, find the A matrix of dimession i x 1
#
# PRE:
# - a kernel cube with all the K_h(pi-p0) values
# - a polinomial cube with all the (pi-p0)^n values
# - a vector with all the phi(xis) values
#
# INPUT:
#
# p0Index      - An index number refering to the p0 value
# hIndex       - An index number refering to the bandwith
# phiVector    - The vector with the phi values
# kernelCube   - The cube with all the Kh(pi-p0) values
# poliCube     - The cube with all the (pi-p0)^n values

# OUTPUT:
#
# Return a vector with the (A0,A1,A2, A3) values

findBigAvalues <- function(p0Index,hIndex,phiVector,kernelCube,poliCube){
  
  A0 = findBigA(p0Index,hIndex,1,phiVector,kernelCube,poliCube)
  A1 = findBigA(p0Index,hIndex,2,phiVector,kernelCube,poliCube)
  A2 = findBigA(p0Index,hIndex,3,phiVector,kernelCube,poliCube)
  A3 = findBigA(p0Index,hIndex,4,phiVector,kernelCube,poliCube)
  
  return(c(A0,A1,A2,A3))
  
}


# Finds a single a_r(p0). To use in the next function that finds all a's. (lowercase)
#
# PRE:
# - a kernel cube with all the K_h(pi-p0) values
# - a polinomial cube with all the (pi-p0)^n values
# - a vector with all the phi(xis) values
#
# INPUT:
#
# p0           - An index number refering to the p0 value
# h            - An index number refering to the bandwith
# r            - TODO: Name? index with the ^r value. 1=0, 2=1, ... 5=4
# kernelCube   - The cube with all the Kh(pi-p0) values
# poliCube     - The cube with all the (pi-p0)^n values
#
# OUTPUT:
#
# A single float with the Ar(p0) 
findLittleA <- function(p0,h,r,kernelCube,poliCube){
  
  myAr = sum(poliCube[[r]][,p0]*kernelCube[[h]][,p0])
  
  return(myAr)
  
}



# Find as (lowercase)
#
# For a given set of data, find the a matrix of dimession n x n
#
# PRE:
# - a kernel cube with all the K_h(pi-p0) values
# - a polinomial cube with all the (pi-p0)^n values
# - a vector with all the phi(xis) values
#
# INPUT:
#
# p0Index      - An index number refering to the p0 value
# hIndex       - An index number refering to the bandwith
# kernelCube   - The cube with all the Kh(pi-p0) values
# poliCube     - The cube with all the (pi-p0)^n values
#
# OUTPUT:
#
# Return a vector with the (a0,a1,a2,a3,a4) values

findLittleAvalues <- function(p0Index,hIndex,kernelCube,poliCube){
  
  a0 = findLittleA(p0Index,hIndex,1,kernelCube,poliCube)
  a1 = findLittleA(p0Index,hIndex,2,kernelCube,poliCube)
  a2 = findLittleA(p0Index,hIndex,3,kernelCube,poliCube)
  a3 = findLittleA(p0Index,hIndex,4,kernelCube,poliCube)
  a4 = findLittleA(p0Index,hIndex,5,kernelCube,poliCube)
  a5 = findLittleA(p0Index,hIndex,6,kernelCube,poliCube)
  a6 = findLittleA(p0Index,hIndex,7,kernelCube,poliCube)
  
  return(c(a0,a1,a2,a3,a4,a5,a6))
  
}


# Find thetas and quadratic method
#
# For a given set of data, find the variable matrix solution of theta0, theta1, and theta2
#
# INPUT:
#
# bigAs     - vector with the big As associated with that p0,x,p,h and kernel
# littleAs  - vector with the little As associated with that p0, p,h and kernel
#
# OUTPUT:
#
# c(NaN, NaN, NaN)             - If is not a Cramer system
# c(theta0, theta1, theta2)    - Otherwise
#
# Note that selecting the quadratic algorithm, will not return a valid theta3 value.

findThetas <- function(bigAs, littleAs){
    
    # Init stuff
    myThetas      = c(NaN,NaN,NaN)
    matrixLittleA = 0
    matrixTheta0  = 0
    matrixTheta1  = 0
    matrixTheta2  = 0
    
    # First, generate the matrix with the little a values
    
    matrixLittleA = matrix(c(littleAs[1],littleAs[2], littleAs[3],
                             littleAs[2],littleAs[3], littleAs[4],
                             littleAs[3],littleAs[4], littleAs[5]),3,3, byrow = TRUE)

    # Check out that we have indeed a Cramer system
    cramerValue = det(matrixLittleA)
    
    # Correct for very tiny numbers
    # if(abs(cramerValue)<0.001){ #TODO: 0.001 as it says in ??????
    #     
    #     cramerValue = 0
    # }
    
    if(cramerValue !=0){
        
        # Make the matrix for eath theta
        matrixTheta0 = matrix(c(bigAs[1],    littleAs[2], littleAs[3],
                                bigAs[2],    littleAs[3], littleAs[4],
                                bigAs[3],    littleAs[4], littleAs[5]),3,3, byrow = TRUE)
            
        matrixTheta1 = matrix(c(littleAs[1], bigAs[1],    littleAs[3],
                                littleAs[2], bigAs[2],    littleAs[4],
                                littleAs[3], bigAs[3],    littleAs[5]),3,3, byrow = TRUE)
            
        matrixTheta2 = matrix(c(littleAs[1], littleAs[2], bigAs[1],
                                littleAs[2], littleAs[3], bigAs[2],
                                littleAs[3], littleAs[4], bigAs[3]),3,3, byrow = TRUE)

        # Find the final thetas
        theta0 = det(matrixTheta0)/cramerValue
        theta1 = det(matrixTheta1)/cramerValue
        theta2 = det(matrixTheta2)/cramerValue
        
        # Cut the values that are smaller than 0 to 0
        if(theta0<0) theta0 = 0
        if (theta0>1) theta0 = 1
        # if(theta1<0) theta1 = 0
        # if(theta2<0) theta2 = 0
        # 
        # 
        myThetas = c(theta0, theta1, theta2)
        
    }
    
    return(myThetas)
    
}










# Find thetas for a cubic method
#
# For a given set of data, find the variable matrix solution of theta0, theta1, theta2, and theta3
#
# INPUT:
#
# bigAs     - vector with the big As associated with that p0,x,p,h and kernel
# littleAs  - vector with the little As associated with that p0, p,h and kernel
#
# OUTPUT:
#
# c(NaN, NaN, NaN)          - If is not a Cramer system
# c(theta0, theta1, theta2) - Otherwise
#
# Note that selecting the quadratic algorithm, will not return a valid theta3 value.

findThetasCubic <- function(littleAs, kernelCube, poliCube, hIndex, p0Index,phiVector){
    
    # Init stuff
    myThetas      = c(NaN,NaN,NaN)
    
    theta0 = 0
    theta1 = 0
    theta2 = 0
    theta3 = 0
    
    numerator   = 0
    denominator = 0
    
    delta00 = 0
    delta01 = 0
    delta02 = 0
    delta03 = 0
    
    delta10 = 0
    delta11 = 0
    delta12 = 0
    delta13 = 0
    
    delta20 = 0
    delta21 = 0
    delta22 = 0
    delta23 = 0
    
    deltaComplete = 0
    
    matrix00 = 0
    matrix01 = 0
    matrix02 = 0
    matrix03 = 0
    
    matrix10 = 0
    matrix11 = 0
    matrix12 = 0
    matrix13 = 0
    
    matrix20 = 0
    matrix21 = 0
    matrix22 = 0
    matrix23 = 0
    
    matrixComplete = 0
 
    # Grab all a's values
    a0 = littleAs[1]
    a1 = littleAs[2]
    a2 = littleAs[3]
    a3 = littleAs[4]
    a4 = littleAs[5]
    a5 = littleAs[6]
    a6 = littleAs[7]
    
    
    # First of all, check that we have a cramer system
    matrixComplete = matrix(c( a0, a1, a2, a3,
                               a1, a2, a3, a4,
                               a2, a3, a4, a5,
                               a3, a4, a5, a6 ),4,4, byrow = TRUE)
    
    deltaComplete = det(matrixComplete)
    
    denominator = deltaComplete
    
    
    # Correct for very tiny numbers
    # if(abs(denominator)<0.001){ #TODO: 0.001 as it says in ??????
    #     
    #     denominator = 0
    # }
    # 
    if(denominator !=0){

        # Find all maxtrices and deltas determinants
        
        matrix00 = matrix(c( a2, a3, a4,
                             a3, a4, a5,
                             a4, a5, a6 ),3,3, byrow = TRUE)
        
        matrix01 = matrix(c( a1, a2, a3,
                             a3, a4, a5,
                             a4, a5, a6 ),3,3, byrow = TRUE)
        
        matrix02 = matrix(c( a1, a2, a3,
                             a2, a3, a4,
                             a4, a5, a6 ),3,3, byrow = TRUE)
        
        matrix03 = matrix(c( a1, a2, a3,
                             a2, a3, a4,
                             a3, a4, a5 ),3,3, byrow = TRUE)
        
        
        matrix10 = matrix(c( a1, a3, a4,
                             a2, a4, a5,
                             a3, a5, a6 ),3,3, byrow = TRUE)
        
        matrix11 = matrix(c( a0, a2, a3,
                             a2, a4, a5,
                             a3, a5, a6 ),3,3, byrow = TRUE)
        
        matrix12 = matrix(c( a0, a2, a3,
                             a1, a3, a4,
                             a3, a5, a6 ),3,3, byrow = TRUE)
        
        matrix13 = matrix(c( a0, a2, a3,
                             a1, a3, a4,
                             a2, a4, a5 ),3,3, byrow = TRUE)
        
        
        matrix20 = matrix(c( a1, a2, a4,
                             a2, a3, a5,
                             a3, a4, a6 ),3,3, byrow = TRUE)
        
        matrix21 = matrix(c( a0, a1, a3,
                             a2, a3, a5,
                             a3, a4, a6 ),3,3, byrow = TRUE)
        
        matrix22 = matrix(c( a0, a1, a3,
                             a1, a2, a4,
                             a3, a4, a6 ),3,3, byrow = TRUE)
        
        matrix23 = matrix(c( a0, a1, a3,
                             a1, a2, a4,
                             a2, a3, a5 ),3,3, byrow = TRUE)

        delta00 = det(matrix00)
        delta01 = det(matrix01)
        delta02 = det(matrix02)
        delta03 = det(matrix03)
        
        delta10 = det(matrix10)
        delta11 = det(matrix11)
        delta12 = det(matrix12)
        delta13 = det(matrix13)
        
        delta20 = det(matrix20)
        delta21 = det(matrix21)
        delta22 = det(matrix22)
        delta23 = det(matrix23)
        
        # Now we need to cumsum all values for all pis
        totalPis = ncol(poliCube[[1]])

        theta0 = 0
        theta1 = 0
        theta2 = 0

        for(i in 1:totalPis){
            
            
            # Grab all the precalculated values
            kernelValue = kernelCube[[hIndex]][i,p0Index]
            p1Value     = poliCube[[2]][i,p0Index]
            p2Value     = poliCube[[3]][i,p0Index]
            p3Value     = poliCube[[4]][i,p0Index]
            phiValue    = phiVector[i]


            numerator = delta00 - delta01*p1Value + delta02*p2Value - delta03*p3Value
            
            # print(numerator)
            
            theta0    = theta0 + (numerator * kernelValue * phiValue / denominator)
            
            numerator = - delta10 + delta11*p1Value - delta12*p2Value + delta13*p3Value
            theta1    = theta1 + (numerator * kernelValue * phiValue / denominator)
            
            numerator = delta20 - delta21*p1Value + delta22*p2Value - delta23*p3Value
            theta2    = theta2 + (numerator * kernelValue * phiValue / denominator)
            
        }

        # Cut the values that are smaller than 0 to 0
        # TODO: Do the same with greater than 1??????
        
        if(theta0 < 0) theta0 = 0
        if(theta1 < 0) theta1 = 0
        if(theta2 < 0) theta2 = 0
        
        # if(theta0 > 1) theta0 = 1
        # if(theta1 > 1) theta1 = 1
        # if(theta2 > 1) theta2 = 1
        
        
        myThetas = c(theta0, theta1, theta2)
           
    }
 
    return(myThetas)
   
}


# For a given kernelCube and polinomialCube, and a's matrix, generate all possible Kbar_h(pi-p0) combinations.
# This is only for the quadratic version, thus you only need 5 matrices
#
# PRE:
# - a kernel cube with all the K_h(pi-p0) values
# - a polinomial cube with all the (pi-p0)^n values
#
# INPUT:
#
#     kernelCube   - The cube with all the Kh(pi-p0) values
#     poliCube     - The cube with all the (pi-p0)^n values
#     anMatrix     - A matrix with the pi,h values of the a_n's
#     order        - TODO: The order of the kernel cube you want to generate
#
# OUTPUT:
#
# A list of matrixs with all the Kbar_h(pi-p0) combinations pre-calculated.
# Where the index of the list represent the h value, the row is pi and the column p0
#
generateKernelBarHCube <- function(kernelCube, poliCube,
                                   a0Matrix, a1Matrix, a2Matrix, a3Matrix, a4Matrix,
                                   order){
  
  totalPs  = length(kernelCube[[1]][,1])
  totalP0s = length(kernelCube[[1]][1,]) 
  totalHs  = length(kernelCube)

  numerator   = 0
  denominator = 0
  
  matrixByHIndex <- list()
  
  # The index of the list refers to the h that you want.
  # For example, hCube[1][2,3] gives you the kernelH value of h=[1], p-p0=[2][3]

  for(i in 1:totalHs){
    
    # Feedback to the user about how much is done and how much time is left
    print(paste(round(i/xgrid*100,2)," %"))
    previosTime = Sys.time()
    
    # We have a list of matrices, create a new one in position i of the list
    matrixByHIndex[[i]] <- array(0,dim = c(totalPs, totalP0s))
    
    for(j in 1:totalP0s){
      
      for(k in 1:totalPs){

        # Grab all a's values
        a0 = a0Matrix[j,i]
        a1 = a1Matrix[j,i]
        a2 = a2Matrix[j,i]
        a3 = a3Matrix[j,i]
        a4 = a4Matrix[j,i]
        
        # Grab all the precalculated values
        kernelValue = kernelCube[[i]][k,j]
        p1Value     = poliCube[[2]][k,j]  # Remember poliCube[[n]][i,j] is the (pi-pj)^(n-1)
        p2Value     = poliCube[[3]][k,j]
        
        
        # Put everything toguether based on order and save
        if(order == 0){
            
            numerator   = (a1*a3-a2^2)*p2Value + (a2*a3 - a1*a4)*p1Value + (a2*a4-a3^2)
            denominator = a0*a2*a4 + 2*a1*a2*a3 - a2^3 - a0*a3^2 - a1^2*a4
            
        }
        
        if(order == 1){
            
            numerator   = (a1*a2-a1*a3)*p2Value + (a0*a4 - a2^2)*p1Value + (a2*a3-a1*a4)
            denominator = a0*a2*a4 + 2*a1*a2*a3 - a2^3 - a0*a3^2 - a1^2*a4
            
        }
        
        if(order == 2){
        
            numerator   = (a0*a2-a1^2)*p2Value + (a1*a2 - a0*a3)*p1Value + (a1*a3-a2^2)
            denominator = a0*a2*a4 + 2*a1*a2*a3 - a2^3 - a0*a3^2 - a1^2*a4
                
        }
        # denominator is the same for the three orders, maybe not for the quadratic, leave there for now
        
        # Check that we have "enought value"
        if(abs(denominator)<0.001){ #TODO: 0.001 as it says in ??????
          
          denominator = 0
          
        }
        
        # Save it into the matrix
        matrixByHIndex[[i]][k,j] <- kernelValue * numerator/denominator
        
      }
      
    }  
    
    # Find out timing statistics
    currentTime = Sys.time()
    deltaTime = currentTime - previosTime
    deltaTime = round(deltaTime,2)
    
    # Show time to go
    print((xgrid-i) * deltaTime)
    
  }
  
  # Create an empty matrix
  return(matrixByHIndex)
  
}



# For a given kernelCube and polinomialCube, and a's matrix, generate all possible Kbar_h(pi-p0) combinations.
# This is for the cubic version, and thus you need the 7 matrices.
#
# In the paper, this is notated as K tilde
#
# PRE:
# - a kernel cube with all the K_h(pi-p0) values
# - a polinomial cube with all the (pi-p0)^n values
#
# INPUT:
#
#     kernelCube   - The cube with all the Kh(pi-p0) values
#     poliCube     - The cube with all the (pi-p0)^n values
#     phiVector    - The vector that contains all phi(pi) values TODO: Delete this? Not in use now
#     anMatrix     - A matrix with the pi,h values of the a_n's
#     order        - TODO: The order of the kernel cube you want to generate
#
# OUTPUT:
#
# A list of matrixs with all the Kbar_h(pi-p0) combinations pre-calculated.
# Where the index of the list represent the h value, the row is pi and the column p0
#
generateKernelBarHCubeCubic <- function(kernelCube, poliCube, 
                                        a0Matrix, a1Matrix, a2Matrix, a3Matrix, a4Matrix, a5Matrix, a6Matrix,
                                        order){
    
    totalPs  = length(kernelCube[[1]][,1])
    totalP0s = length(kernelCube[[1]][1,]) 
    totalHs  = length(kernelCube)
    
    # Init all the matrices
    numerator   = 0
    denominator = 0
    
    delta00 = 0
    delta01 = 0
    delta02 = 0
    delta03 = 0
    
    delta10 = 0
    delta11 = 0
    delta12 = 0
    delta13 = 0
    
    delta20 = 0
    delta21 = 0
    delta22 = 0
    delta23 = 0
    
    deltaComplete = 0
    
    matrix00 = 0
    matrix01 = 0
    matrix02 = 0
    matrix03 = 0
    
    matrix10 = 0
    matrix11 = 0
    matrix12 = 0
    matrix13 = 0
    
    matrix20 = 0
    matrix21 = 0
    matrix22 = 0
    matrix23 = 0
    
    matrixComplete = 0
    
    matrixByHIndex <- list()
    
    # The index of the list refers to the h that you want.
    # For example, hCube[1][2,3] gives you the kernelH value of h=[1], p-p0=[2][3]
    
    for(i in 1:totalHs){
        
        matrixByHIndex[[i]] <- array(0,dim = c(totalPs, totalP0s))
        
        for(j in 1:totalP0s){
            
            for(k in 1:totalPs){
                
                # Grab all a's values
                a0 = a0Matrix[j,i]
                a1 = a1Matrix[j,i]
                a2 = a2Matrix[j,i]
                a3 = a3Matrix[j,i]
                a4 = a4Matrix[j,i]
                a5 = a5Matrix[j,i]
                a6 = a6Matrix[j,i]
                
                # Grab all the precalculated values
                kernelValue = kernelCube[[i]][k,j]
                p1Value     = poliCube[[2]][k,j]
                p2Value     = poliCube[[3]][k,j]
                p3Value     = poliCube[[4]][k,j]
                
                # Find all maxtrices and deltas determinants
                
                matrix00 = matrix(c( a2, a3, a4,
                                     a3, a4, a5,
                                     a4, a5, a6 ),3,3, byrow = TRUE)
            
                matrix01 = matrix(c( a1, a2, a3,
                                     a3, a4, a5,
                                     a4, a5, a6 ),3,3, byrow = TRUE)
                
                matrix02 = matrix(c( a1, a2, a3,
                                     a2, a3, a4,
                                     a4, a5, a6 ),3,3, byrow = TRUE)
                
                matrix03 = matrix(c( a1, a2, a3,
                                     a2, a3, a4,
                                     a3, a4, a5 ),3,3, byrow = TRUE)
                
                
                matrix10 = matrix(c( a1, a3, a4,
                                     a2, a4, a5,
                                     a3, a5, a6 ),3,3, byrow = TRUE)
                
                matrix11 = matrix(c( a0, a2, a3,
                                     a2, a4, a5,
                                     a3, a5, a6 ),3,3, byrow = TRUE)
                
                matrix12 = matrix(c( a0, a2, a3,
                                     a1, a3, a4,
                                     a3, a5, a6 ),3,3, byrow = TRUE)
                
                matrix13 = matrix(c( a0, a2, a3,
                                     a1, a3, a4,
                                     a2, a4, a5 ),3,3, byrow = TRUE)
                
                
                matrix20 = matrix(c( a1, a2, a4,
                                     a2, a3, a5,
                                     a3, a4, a6 ),3,3, byrow = TRUE)
                
                matrix21 = matrix(c( a0, a1, a3,
                                     a2, a3, a5,
                                     a3, a4, a6 ),3,3, byrow = TRUE)
                
                matrix22 = matrix(c( a0, a1, a3,
                                     a1, a2, a4,
                                     a3, a4, a6 ),3,3, byrow = TRUE)
                
                matrix23 = matrix(c( a0, a1, a3,
                                     a1, a2, a4,
                                     a2, a3, a5 ),3,3, byrow = TRUE)
                
                
                
                matrixComplete = matrix(c( a0, a1, a2, a3,
                                           a1, a2, a3, a4,
                                           a2, a3, a4, a5,
                                           a3, a4, a5, a6 ),4,4, byrow = TRUE)
                
                
                delta00 = det(matrix00)
                delta01 = det(matrix01)
                delta02 = det(matrix02)
                delta03 = det(matrix03)
                
                delta10 = det(matrix10)
                delta11 = det(matrix11)
                delta12 = det(matrix12)
                delta13 = det(matrix13)
                
                delta20 = det(matrix20)
                delta21 = det(matrix21)
                delta22 = det(matrix22)
                delta23 = det(matrix23)
                
                deltaComplete = det(matrixComplete)
                
                # Put everything toguether based on order and save
                if(order == 0){
                    
                    numerator   = delta00 - delta01*p1Value + delta02*p2Value - delta03*p3Value 
                    denominator = deltaComplete
                    
                }
                
                if(order == 1){
                    
                    numerator   = -delta10 + delta11*p1Value + delta12*p2Value + delta13*p3Value
                    denominator =  deltaComplete
                    
                }
                
                if(order == 2){
                    
                    numerator   = 2 * delta20 - delta21*p1Value + delta22*p2Value - delta23*p3Value
                    denominator = deltaComplete
                    
                }
                
                # Check that we have "enought value"
                if(abs(denominator)<0.001){ #TODO: 0.001 as it says in ??????
                    
                    denominator = 0
                    
                }
                
                # Save it into the matrix
                matrixByHIndex[[i]][k,j] <- kernelValue * numerator/denominator
                
            }
            
        }  
        
    }
    
    # Create an empty matrix
    return(matrixByHIndex)
    
}


# For a given kernelCube and polinomialCube, and a's matrix, generate all possible Kbar_h(pi-p0) combinations.
#
# This is deprecated, as the paper doesn't use the Kernel Bar Bar notation anymore.
#
# PRE:
# - a kernel cube with all the K_h(pi-p0) values
# - a polinomial cube with all the (pi-p0)^n values
#
# INPUT:
#
# kernelCube   - The cube with all the Kh(pi-p0) values
# poliCube     - The cube with all the (pi-p0)^n values
# anMatrix     - A matrix with the pi,h values of the a_n's
#
# OUTPUT:
#
# A list of matrixs with all the Kbar_h(pi-p0) combinations pre-calculated.
# Where the index of the list represent the h value, the row is pi and the column p0
#
generateKernelBarBarHCube <- function(kernelBarCube, weightsMatrix){
  
  totalPs  = length(kernelBarCube[[1]][,1])
  totalP0s = length(kernelBarCube[[1]][1,]) 
  totalHs  = length(kernelBarCube)
  
  matrixByHIndex <- list()
  
  # The index of the list refers to the h that you want.
  # For example, hCube[1][2,3] gives you the kernelH value of h=[1], p-p0=[2][3]
  
  for(i in 1:totalHs){
    
    matrixByHIndex[[i]] <- array(0,dim = c(totalPs, totalP0s))
    
    for(j in 1:totalP0s){
      
      for(k in 1:totalPs){
        
          auxiliar = 0
      
          for(t in k:totalPs){
            
            currentWeight      = weightsMatrix[t,j]
            currentKernelValue = kernelBarCube[[i]][t,j]
            
            auxiliar = auxiliar + currentWeight*currentKernelValue
            
          }

          matrixByHIndex[[i]][k,j] = auxiliar
          
      }
      
    }
    
  }
  
  return(matrixByHIndex)
  
}

# Generate a weight matrix object of size n x n with this format
#
# 1    0       0      ... 0
# 1/n (n-1)/n  0      ... 0
# 1/n  1/n    (n-2)/n ... 0
#          ...
# 1/n  1/n     1/n    ... 1/n
#
# INPUT:
# 
# n - Size of the matrix
#
# OUTPUT:
#
# As specified in the description
generateWeightMatrix <- function(n){
  
  
  myWeightMatrix = matrix(0, nrow = n, ncol = n)
  
  # This loop create a triangular 0 above, the main diagonal to 1/n, and triangular 1/n under
  for (i in 1:n) {
    
    for (j in 1:n) {
      
      myWeightMatrix[i,j] = 1/n
      
      if(i==j){
        i = n
      }
      
    }
    
  }
  
  # This loop correct the main diagonal to the actual values
  for (i in 1:n) {
    
    myWeightMatrix[i,i] = myWeightMatrix[i,i] * (n-(i-1))
    
  }
  
  return(myWeightMatrix)
  
}



# Generate a weight matrix object of size n x n 
#
# Where the element i,j has this value:
#
# w_ij = F(x =   i/n   ; alpha = j ; beta = n-j+1) -
#        F(x = (i-1)/n ; alpha = j ; beta = n-j+1) -
#
# Where F is the cumulative distribution given by a Beta distribution
#
#
# INPUT:
# 
# n - Size of the matrix
#
# OUTPUT:
generateBetaWeightMatrix <- function(n){
    
    
    myWeightMatrix = matrix(0, nrow = n, ncol = n)
    
    # This loop create a triangular 0 above, the main diagonal to 1/n, and triangular 1/n under
    for (i in 1:n) {
        
        for (j in 1:n) {
            
            myWeightMatrix[i,j] = pbeta(i/n, j, n-j+1)- pbeta((i-1)/n, j, n-j+1)
                
        }
        
    }

    return(myWeightMatrix)
    
}



#Generate all possible u_i^m moments, where m is 1 and 2
#
# INPUT:
#
# x       - The vector with the values
#
# OUTPUT:
#
# a matrix of size n x 2, where the columns are the moments order, and the row are the i of u_i^m
generateCentralMomentMatrix <- function(x, printFeedback = FALSE){
  
  n = length(x)
  
  myMomentsMatrix = matrix(0, nrow = n, ncol = 2)
  
  # Feedback for the user
  # ---- Timing stuff
  previosTime = Sys.time()
  currentTime = Sys.time()
  
  if(printFeedback == TRUE){
    
    print("Generating central moment matrix")
    print(n)
    print("--------------------------------")    

  }

  
  for (i in 1:n) {
    
    if(printFeedback == TRUE){
    
        # Feedback to the user about how much is done and how much time is left
        print(paste(i/n*100," %"))
        previosTime = Sys.time()
        # --------------------------------------------------------------------
    
    }
    
    combinationNumber = choose(n-1, i-1)
    
    for (m in 1:2) {
      
      auxiliar = 0
        
      for (r in 1:n) {
        
        px  = x[r]^m
        pr  = (r/n)^(i-1)
        p1r = (1-r/n)^(n-i)
        
        subTotal = px * pr * p1r

        auxiliar = auxiliar + subTotal
        
      }
      
      myMomentsMatrix[i,m] = combinationNumber * auxiliar
      
    }
    
    if(printFeedback == TRUE){
    
        # Find out timing statistics
        currentTime = Sys.time()
        deltaTime   = currentTime - previosTime
    
        # Show time to go
        print((n-i) * deltaTime)
    }
    
  }
  
  return(myMomentsMatrix)
  
}


# Generate all possible u_i^m moments, where m is 1 and 2
# Using the Beta method
#
# INPUT:
#
# x       - The vector with the values
#
# OUTPUT:
#
# a matrix of size n x 2, where the columns are the moments order, and the row are the i of u_i^m
generateCentralMomentMatrixBETA <- function(x, printFeedback = FALSE){
    
    n = length(x)
    myWeights = generateBetaWeightMatrix(n)
    myMomentsMatrix = matrix(0, nrow = n, ncol = 2)
    
    # Feedback for the user
    # ---- Timing stuff
    previosTime = Sys.time()
    currentTime = Sys.time()
    
    if(printFeedback == TRUE){
        
        print("Generating central moment matrix")
        print(n)
        print("--------------------------------")    
        
    }
    
    
    for (i in 1:n) {
        
        if(printFeedback == TRUE){
            
            # Feedback to the user about how much is done and how much time is left
            print(paste(i/n*100," %"))
            previosTime = Sys.time()
            # --------------------------------------------------------------------
            
        }
        
        for (m in 1:2) {
            
            auxiliar = 0
            
            for (r in 1:n) {

                auxiliar = auxiliar + myWeights[i,r] * x[r]^m
                
            }
            
            myMomentsMatrix[i,m] = auxiliar
            
        }
        
        if(printFeedback == TRUE){
            
            # Find out timing statistics
            currentTime = Sys.time()
            deltaTime   = currentTime - previosTime
            
            # Show time to go
            print((n-i) * deltaTime)
        }
        
    }
    
    return(myMomentsMatrix)
    
}


# Generate all possible variance of Xi values
#
# INPUT:
#
# x       - The vector with the values
#
# OUTPUT:
#
# A vector of size length(x) with the variance for each Xi
generateVarianceXVector <- function(x){
  
  # Init stuff
  n = length(x)
  myMomentsMatrix  = generateCentralMomentMatrixBETA(x)
  myVarianceVector = 1:n

  for (i in 1:n) {
    
    myVarianceVector[i] = myMomentsMatrix[i,2] - myMomentsMatrix[i,1]^2
    
  }
  
  return(myVarianceVector)

}

# Generate all possible variance of Xi values
#
# INPUT:
#
# x       - The vector with the values
#
# OUTPUT:
#
# A vector of size length(x) with the variance for each Xi
generateVarianceXVectorBOOTSRAP <- function(x, bootFactor = 100){
    
    # Init stuff
    n                = length(x)
    myVarianceVector = 1:n
    
    # print("Pre function")
    # print(myVarianceVector)
    
    # This first loop goes throught all the variance vector that needs to be returned
    for (i in 1:n) {
        
        # xcopy     = x    
        auxVector = 1:bootFactor
    
        # print("Pre auxVector")
        # print(auxVector)
        
        # This second loops does the bootstraps
        
        # xcopy = sample( xcopy , n , replace=T )
        # xcopySorted = sort(xcopy)
        
        for (j in 1:bootFactor) {
            
            xcopy       = x    
            xcopy       = sample( xcopy , n , replace=TRUE )
            xcopySorted = sort(xcopy)
          
            auxVector[j] = xcopySorted[i]
        
        }
        # 
        # print("post bootstrap")
        # print(auxVector)
        
        myVarianceVector[i] = var( auxVector )
        
    }
    
    # print("Inside function")
    # print(myVarianceVector)
    
    return(myVarianceVector)
    
}


#Generate all possible u_i_j^1 moments
#
# INPUT:
#
# x       - The vector with the values
#
# OUTPUT:
#
# a matrix of size n x n, where the main diagonal and the lower triangle are 0
# and the upper triangle is the u_i_j product moment

generateProductMomentMatrix <- function(x, printFeedback = FALSE){
  
  n = length(x)
  
  myMomentsMatrix = matrix(0, nrow = n, ncol = n)
  
  # Feedback for the user
  # ---- Timing stuff
  previosTime = Sys.time()
  currentTime = Sys.time()
  
  if(printFeedback == TRUE){
    
    print("Generating product moment matrix")
    print(n)
    print("--------------------------------")    
    
  }
  
  
  for (i in 1:n) {
  
    if(printFeedback == TRUE){
      
        # Feedback to the user about how much is done and how much time is left
        print(paste(i/n*100," %"))
        previosTime = Sys.time()
    
    }
    
    for (j in 1:n) {
      
      if(i<j){
        
        # TODO: Factorial for big numbers doesn't work
        # TODO: n > 170
        
        factorialNumber = factorial(n)/(factorial(i-1)*factorial(j-i-1)*factorial(n-j))
        
        auxiliarR = 0
        auxiliarS = 0
        
        for (r in 1:n) { #TODO: Check limits
          
          auxiliarS = 0

          for (s in 1:r) {
                
                if(s<=r){

                    Xs   = x[r]*x[s]/n^2
                    ps   = (s/n)^(i-1)
                    prps = ((r-s)/n)^(j-i-1)
                    pr   = (1-r/n)^(n-j)
                    
                    subTotal = Xs * ps * prps * pr
                    
                    auxiliarS = auxiliarS + subTotal

                }
            }

          auxiliarR = auxiliarR + auxiliarS
          
        }

        myMomentsMatrix[i,j] = factorialNumber * auxiliarR
        
      }
      
    }
    
    if(printFeedback == TRUE){
    
      # Find out timing statistics
      currentTime = Sys.time()
      deltaTime   = currentTime - previosTime
      
      # Show time to go
      print((n-i) * deltaTime)
        
    }
    
  }
  
  return(myMomentsMatrix)
  
}


# Generate all possible u_i_j^1 moments
# Using the beta method
#
# INPUT:
#
# x       - The vector with the values
#
# OUTPUT:
#
# a matrix of size n x n, where the main diagonal and the lower triangle are 0
# and the upper triangle is the u_i_j product moment

generateProductMomentMatrixBETA <- function(x, printFeedback = FALSE){
    
    n = length(x)
    
    myMomentsMatrix = matrix(0, nrow = n, ncol = n)
    
    # Feedback for the user
    # ---- Timing stuff
    previosTime = Sys.time()
    currentTime = Sys.time()
    
    if(printFeedback == TRUE){
        
        print("Generating product moment matrix")
        print(n)
        print("--------------------------------")    
        
    }
    
    
    for (i in 1:n) {
        
        if(printFeedback == TRUE){
            
            # Feedback to the user about how much is done and how much time is left
            print(paste(i/n*100," %"))
            previosTime = Sys.time()
            
        }
        
        for (j in 1:n) {
            
            if(i<j){

                auxiliarR = 0
                auxiliarS = 0
                
                for (r in 1:n) { #TODO: Check limits
                    
                    auxiliarS = 0
                    
                    for (s in 1:r) {
                        
                        if(s<=r){

                            auxiliarS = auxiliarS + 0

                        }
                    }

                    auxiliarR = auxiliarR + auxiliarS
                    
                }
                
                myMomentsMatrix[i,j] = auxiliarR
                
            }
            
        }
        
        if(printFeedback == TRUE){
            
            # Find out timing statistics
            currentTime = Sys.time()
            deltaTime   = currentTime - previosTime
            
            
            # Show time to go
            print((n-i) * deltaTime)
            
        }
        
    }
    
    return(myMomentsMatrix)
    
}


# TODO: Description
#
# INPUT:
# x       - The vector with the values
#
# OUTPUT:
#
# A matrix with the main diagonal and the lower triangle equal to 0
# and the upper triangle with the covariance of XiXj where i is the row and j the column
#
generateCovarianceXMatrix <- function(x, printFeedback = FALSE){
  
  # Init stuff
  n = length(x)
  myCentralMomentsMatrix  = generateCentralMomentMatrixBETA(x, printFeedback)
  myProdcutMomentsMatrix  = generateProductMomentMatrixBETA(x, printFeedback)
  
  myCovarianceMatrix = matrix(0, nrow = n, ncol = n)
  
  for (i in 1:n) {
    
    for (j in 1:n) {
      
      if(i<j){
    
        myCovarianceMatrix[i,j] = myProdcutMomentsMatrix[i,j] - myCentralMomentsMatrix[i,1]*myCentralMomentsMatrix[j,1]
            
      }
    }
    
  }
  
  return(myCovarianceMatrix)
}


# TODO: Description
#
# INPUT:
# x       - The vector with the values
#
# OUTPUT:
#
# A matrix with the main diagonal and the lower triangle equal to 0
# and the upper triangle with the covariance of XiXj where i is the row and j the column
#
generateCovarianceXMatrixBOOTSTRAP <- function(x, bootFactor = 100, printFeedback = TRUE){
    
    # Init stuff
    n                  = length(x)
    myCovarianceMatrix = matrix(0, nrow = n, ncol = n)
 
    # Feedback for the user
    # ---- Timing stuff
    previosTime = Sys.time()
    currentTime = Sys.time()
    
    
    if(printFeedback == TRUE){
        
        print("Generating covariance matrix")
        print(n)
        print("--------------------------------")    
        
    }
    
    
    # The first two loops goes throught the covariance matrix
    for (i in 1:n) {
        
        if(printFeedback == TRUE){
            
            # Feedback to the user about how much is done and how much time is left
            print(paste(round(i/n*100,2)," %"))
            previosTime = Sys.time()
            
        }
        
        for (j in 1:n) {
            
            if(i<j){
                
                # xcopy        = x
                # 
                # firstVector  = x
                # secondVector = x
                    
                auxVectorX   = 1:bootFactor
                auxVectorY   = 1:bootFactor

                # firstVector   = sample( xcopy , n , replace=T ) #TODO: Also check this is ok too
                # secondVector  = sample( xcopy , n , replace=T )
                # 
                # firstVectorSorted  = sort(firstVector)
                # secondVectorSorted = sort(secondVector)
                                
                # The third loops, does the bootstrap
                for (k in 1:bootFactor) {
                    
                    xcopy        = x
                  
                    firstVector  = x
                    secondVector = x
                    
                    firstVector   = sample( xcopy , n , replace=TRUE ) #TODO: Also check this is ok too
                    secondVector  = sample( xcopy , n , replace=TRUE )
                    
                    firstVectorSorted  = sort(firstVector)
                    secondVectorSorted = sort(secondVector)
                  
                    auxVectorX[k] = firstVectorSorted[i]
                    auxVectorY[k] = secondVectorSorted[j]
                    
                }
                
                myCovarianceMatrix[i,j] = cov(auxVectorX,auxVectorY)

            }
        }
    
        
        
        if(printFeedback == TRUE){
            
            # Find out timing statistics
            currentTime = Sys.time()
            deltaTime   = currentTime - previosTime
            
            # Show time to go
            print((n-i) * deltaTime)
            
        }
        
    }
    
    return(myCovarianceMatrix)
}

# Generate the symetric matrix of covariances as specify in (14)
generateSigmaMatrix <- function(varianceVector, covarianceMatrix){
  
  mySigmaMatrix = matrix(0, nrow = n, ncol = n)
  
  
  for (i in 1:n) {
    
    for (j in 1:n) {
      
      # For the main diagonal cases
      if(i == j){
      
        mySigmaMatrix[i,j] = varianceVector[i]
          
      }
      
      # For the upper and lower triangles
      else{
        
        if(i<j){
          
          mySigmaMatrix[i,j] = covarianceMatrix[i,j]
          mySigmaMatrix[j,i] = covarianceMatrix[i,j]
          
        }

      }

    }
    
  }
  
  # # This loop correct the main diagonal to the actual values
  # for (i in 1:n) {
  #   
  #   mySigmaMatrix[i,i] = mySigmaMatrix[i,i] * (n-(i-1))
  #   
  # }
  # 
   return(mySigmaMatrix)
  
}





# Find the variance for a given derivate with a given kernelCube
#
# TODO: Fulldoc
# h is an index
# p0 is an index

findVariancePhiWithKernel <- function(h,p0, kernelBarBarCube, varianceVector, covarianceMatrix, order){
  
  totalPs  = length(kernelBarBarCube[[1]][,1])
  totalP0s = length(kernelBarBarCube[[1]][1,]) 
  totalHs  = length(kernelBarBarCube)
  
  totalValue = 0
  
  auxiliarRight = 0
  auxiliarLeft  = 0
  
  # print(kernelBarBarCube)
  # print(varianceVector)
  # print()
  # print(totalPs)
  # print(totalP0s)
  # print(totalHs)
  # 
  # stop()
  
  # Do first, the left side
  for (j in 1:totalPs) {

      auxiliarLeft = auxiliarLeft + (kernelBarBarCube[[h]][j,p0]^2)*varianceVector[j]
    
  }

  # Do now, the right side
  for (i in 1:totalPs) {
    
    for (j in 1:totalPs) {
      
      if(i<j){ 

        kernelValueHI   = kernelBarBarCube[[h]][i,p0]
        kernelValueHJ   = kernelBarBarCube[[h]][j,p0]
        covarianceValue = covarianceMatrix[i,j]

        auxiliarRight = auxiliarRight + kernelValueHI * kernelValueHJ * covarianceValue

      }
      
    }
    
  }
  
  # Finally, put it toghether based on the order
  if(order == 0){
      
      totalValue = auxiliarLeft + 2*auxiliarRight
      
  }
  
  if(order == 1){ #TODO: order 1 and order 2 are equals????
      
      totalValue = auxiliarLeft + 2*auxiliarRight
      
  }
  
  if(order == 2){
  
      totalValue = 4*auxiliarLeft + 8*auxiliarRight
          
  }
  
  return(totalValue)

}






findVariancePhiWithKernelNEW <- function(h,p0, kernelBarCube, sigmaMatrix, weightMatrix, order){
  
  
  returnValue  = 0
  
  transposeWeightMatrix = t(weightMatrix)

  # This is a Pi x 1 Matrix
  kernelColumn = t(t(kernelBarCube[[h]][,p0])) # kernelBarCube[[h]][,p0] is a column, but R return it as an array.
                                               # In order to get a matrix, the simplest thing is to transpose it twice.
                                               # If you transpose it once, you get a 1 x Pi Matrix
                                               # If you transpose it again, you get your column as Pi x 1
  # This is a 1 x Pi Matrix
  kernelRow    = t(kernelBarCube[[h]][,p0])

  if(order != 2){
    
    returnValue = kernelRow %*% weightMatrix %*% sigmaMatrix %*% transposeWeightMatrix %*% kernelColumn  
    
  }
  else{
    
    returnValue = 4 * kernelRow %*% weightMatrix %*% sigmaMatrix %*% transposeWeightMatrix %*% kernelColumn  
    
  }
  
  return (returnValue)
  
}



findVariancePhiWithKernelCubicNEW <- function(h,p0, kernelBarCube, sigmaMatrix, weightMatrix, order){
  
  
  returnValue  = 0
  
  transposeWeightMatrix = t(weightMatrix)
  
  # This is a Pi x 1 Matrix
  kernelColumn = t(t(kernelBarCube[[h]][,p0])) # kernelBarCube[[h]][,p0] is a column, but R return it as an array.
                                               # In order to get a matrix, the simplest thing is to transpose it twice.
                                               # If you transpose it once, you get a 1 x Pi Matrix
                                               # If you transpose it again, you get your column as Pi x 1
                                              
  # This is a 1 x Pi Matrix
  kernelRow    = t(kernelBarCube[[h]][,p0])
  
  returnValue = kernelRow %*% weightMatrix %*% sigmaMatrix %*% transposeWeightMatrix %*% kernelColumn
  
  return (returnValue)
  
}


####################################
# From here, only functions to log
# system stuff and debugging
####################################


# Gets a cube-like data structure and return a human readable string
#
kernelToString <- function(cube, mainHeader, subHeader, columnsString, rowsNameVector){
  
  # -- Write the main Header
  cubeString       = ""
  cubeHeaderString = mainHeader
  
  # For each matrix of the cube
  #
  for (i in 1:length(cube)) {
    
    # Write the sub-header and the columns names
    cubeString = paste(cubeString, subHeader,i, sep = '')
    cubeString = paste(cubeString, columnsString, sep = '\n')
   
    # For each row of the matrix 
    for (j in 1:length(rowsNameVector)) {
      
      rowString  = paste(j, " | ", round(rowsNameVector[j],4), " | ", sep = '')
      
      for (k in 1:length(cube[[i]][j,])) {
        
        cellString = paste(round(cube[[i]][j,k],4), " ", sep = '')  
        rowString  = paste(rowString, cellString, sep = '') 
        
      }
      
      cubeString = paste(cubeString, rowString, sep = '\n') 
      
    }
    
    cubeString = paste(cubeString, "    ----------    ", sep = '\n\n') 
    
  }
  
  return( cubeString )
  
}

# Gets a matrix-like data structure and return a human readable string
matrixToString <- function(myMatrix, myMatrixHeaderString){

  totalRows    = nrow(myMatrix)
  totalColumns = ncol(myMatrix)
  
  myMatrixString = paste(myMatrixHeaderString,"\n",sep = '')
  
  
  for (i in 1:totalRows) {
    
    rowString  = paste(i, "  | ", sep = '')
    
    for (j in 1:totalColumns) {

      cellString = paste(round(myMatrix[i,j],4), " ", sep = '')  
      
      rowString = paste(rowString, cellString, sep = '') 
      
    }
    
    myMatrixString = paste(myMatrixString, rowString, sep = '\n') 
    
  }
  
  return(myMatrixString)
  
}

# Writes everything into human reading format
#
writeLog <- function(logPath, n, h, kernel, xgrid, ygrid, data, frequencies, CIData,
                     weights, phis, poliCube, kernelCube, varianceVector, covarianceMatrix,
                     a0Matrix, a1Matrix, a2Matrix, a4Matrix, kernelBarCube, kernelBarBarCube,
                     A0Matrix, A1Matrix, A2Matrix){
  
  # Open the FD
  #textFile = paste("log.txt", sep = "")
  fileConn<-file(logPath)
  
  
  # Formating strings
  blankString        = ""
  separatorString    = "----------------------------------------------"

  # Get all p0s into one string that will be the row several times later
  allP0s    = round(CIData$p0[1:xgrid],4)
  p0sString = ""
  for (i in 1:length(allP0s)) {
    
    p0sString = paste(p0sString, " | ", allP0s[i])  
    
  }
  
  #Write  Basic input data
  totalSamplesString = paste("Total samples: ",n,sep = '')
  totalHsString      = paste("Total hs:      ",h,sep = '')
  kernelTypeString   = paste("Kernel in use: ",kernel,sep = '')
  dimensionsString   = paste("SiZer plot dimensions: ", xgrid, " x " , ygrid, sep = '')
  
  # Write samples and frecuencies
  # --------------------------------------------------
  
  dataString = ""
  
  # -- Sort the data
  sortedXData = sort(data)
  
  # -- Write the header
  headerString = "row | -- xi -- | -- pi --"

  # -- For each data and frequency, write it next to each other
  for (i in 1:n) {
    
    separatorString = ""
    separatorStringSize = 4 - nchar(i)
    for (j in 1:separatorStringSize) { #TODO: Overkill for , do better
      separatorString = paste(separatorString, " ", sep = '')
    }
    
    rowString  = paste(i, separatorString, "|  ", round(sortedXData[i],4), "  | ", round(frequencies[i],4) , sep = '')
    dataString = paste(dataString, rowString, sep = '\n') 
 
    
  }
  
  # Write the selections we takes from samples and frequencies
  # The h's selected, average, variance, and so on:
  # ----------------------------------------------------------
  
  selectedDataString   = ""
  
  selectedHeaderString = "row | -- p0 -- | -- h -- | -- Phi2 -- | -- Variance(p0) -- "
  
  # -- For each data and frequency, write it next to each other
  for (i in 1:nrow(CIData)) {
    
    separatorString = ""
    separatorStringSize = 4 - nchar(i)
    for (j in 1:separatorStringSize) { #TODO: Overkill for, do better
      separatorString = paste(separatorString, " ", sep = '')
    }
    
    currentP0       = "NA/NaN"
    currentH        = "NA/NaN"
    currentPhiTwo   = "NA/NaN"
    currentVariance = "NA/NaN"
    
    if( !is.na(CIData$p0[i]) & !is.nan(CIData$p0[i]) )            currentP0       = round(CIData$p0[i],4)  
    if( !is.na(CIData$h[i]) & !is.nan(CIData$h[i]) )              currentH        = CIData$h[i]  
    if( !is.na(CIData$phiTwo[i]) & !is.nan(CIData$phiTwo[i]) )    currentPhiTwo   = round(CIData$phiTwo[i],4)  
    if( !is.na(CIData$variance[i]) & !is.nan(CIData$variance[i])) currentVariance = round(CIData$variance[i],4)  
    
    rowString  = paste(i, separatorString, "| ", currentP0 , " | ", currentH , " | ", currentPhiTwo ," | ", currentVariance , sep = '')
    
    selectedDataString = paste(selectedDataString, rowString, sep = '\n') 
    
    
  }
  
  # Write all the constants data
  # --------------------------------------------------
  
  # -- Write the weight matrix
  weightMatrixString  = ""
  weigthsHeaderString = " Generated Weights"
  
  for (i in 1:n) {
    
    rowString  = paste(i, "  | ", sep = '')
    
    for (j in 1:n) {
    
      
      cellString = paste(weights[i,j], " ", sep = '')  
    
      rowString = paste(rowString, cellString, sep = '') 
        
    }
    
    weightMatrixString = paste(weightMatrixString, rowString, sep = '\n') 

  }

  # -- Write the phi vector  
  phisString       = ""
  phisHeaderString = "Phi vector"
  
  for (i in 1:n) {
    
    rowString  = paste(i, "  | ", round(phis[i],4), sep = '')
    
    phisString = paste(phisString, rowString, sep = '\n') 
    
    
  }
  
  # -- Write the polinomial cube
  poliString          = ""
  poliHeaderString    = "Polinomial Cube"
  poliSubHeaderString = " -- Matrix for polinomials (pi-p0)^(i-1) for i = "
  poliString          = kernelToString(poliCube, poliHeaderString, poliSubHeaderString, p0sString, frequencies)
  
  # Write the kernel h cube
  kernelString          = ""
  kernelHeaderString    = "Kernel H Cube"
  kernelSubHeaderString = " -- Matrix Kh = ((pi-p0)/h)/h for h = "
  kernelString          = kernelToString(kernelCube, kernelHeaderString, kernelSubHeaderString, p0sString, frequencies)
  
  # Write the kernel bar h cube
  kernelBarString          = ""
  kernelBarHeaderString    = "Kernel Bar H Cube"
  kernelBarSubHeaderString = " -- Matrix KBarh = (pi-p0) for h = "
  kernelBarString          = kernelToString(kernelBarCube, kernelBarHeaderString, kernelBarSubHeaderString, p0sString, frequencies)
  
  # Write the kernel h cube
  kernelBarBarString          = ""
  kernelBarBarHeaderString    = "Kernel BarBar H Cube"
  kernelBarBarSubHeaderString = " -- Matrix KBar_hj = (p0) for h = "
  kernelBarBarString          = kernelToString(kernelBarBarCube, kernelBarBarHeaderString, kernelBarBarSubHeaderString, p0sString, frequencies)
  
  
  # Write variance vector
  # --------------------------------------------------
  
  varianceString       = ""
  varianceHeaderString = "Variance vector"
  
  for (i in 1:n) {
    
    rowString  = paste(i, " = ", round(sortedXData[i],4),"  | ", round(varianceVector[i],4), sep = '')
    
    varianceString = paste(varianceString, rowString, sep = '\n') 
    
    
  }
  
  # Write covariance matrix
  # --------------------------------------------------
  
  # -- Write the weight matrix
  coovarianceMatrixString = ""
  covarianceHeaderString  = " Generated Coovariances"
  
  for (i in 1:n) {
    
    rowString  = paste(i, "  | ", sep = '')
    
    for (j in 1:n) {
      
      
      cellString = paste(round(covarianceMatrix[i,j],4), " ", sep = '')  
      
      rowString = paste(rowString, cellString, sep = '') 
      
    }
    
    coovarianceMatrixString = paste(coovarianceMatrixString, rowString, sep = '\n') 
    
  }
  
  
  # Write A's and a's matrices
  # --------------------------------------------------
  
  A0MatrixString = matrixToString(A0Matrix, "A0 Matrix p0s - hs")
  A1MatrixString = matrixToString(A1Matrix, "A1 Matrix p0s - hs")
  A2MatrixString = matrixToString(A2Matrix, "A2 Matrix p0s - hs")
  a0MatrixString = matrixToString(a0Matrix, "a0 Matrix p0s - hs")
  a1MatrixString = matrixToString(a1Matrix, "a1 Matrix p0s - hs")
  a2MatrixString = matrixToString(a2Matrix, "a2 Matrix p0s - hs")
  a3MatrixString = matrixToString(a3Matrix, "a3 Matrix p0s - hs")
  a4MatrixString = matrixToString(a4Matrix, "a4 Matrix p0s - hs")
  
  
  
  
  # Write everything in the log
  # --------------------------------------------------
  writeLines(c(totalSamplesString,
               totalHsString,
               kernelTypeString,
               dimensionsString,
               blankString,
               headerString,
               separatorString,
               dataString,
               blankString,
               selectedHeaderString,
               separatorString,
               selectedDataString,
               blankString,
               weigthsHeaderString,
               separatorString,
               weightMatrixString,
               blankString,
               phisHeaderString,
               separatorString,
               phisString,
               blankString,
               varianceHeaderString,
               separatorString,
               varianceString,
               blankString,
               covarianceHeaderString,
               separatorString,
               coovarianceMatrixString,
               blankString,
               poliHeaderString,
               separatorString,
               poliString,
               blankString,
               kernelHeaderString,
               separatorString,
               kernelString,
               blankString,
               kernelBarHeaderString,
               separatorString,
               kernelBarString,
               blankString,
               kernelBarBarHeaderString,
               separatorString,
               kernelBarBarString,
               blankString,
               separatorString,
               A0MatrixString,
               blankString,
               separatorString,
               A1MatrixString,
               blankString,
               separatorString,
               A2MatrixString,
               blankString,
               separatorString,
               a0MatrixString,
               blankString,
               separatorString,
               a1MatrixString,
               blankString,
               separatorString,
               a2MatrixString,
               blankString,
               separatorString,
               a3MatrixString,
               blankString,
               separatorString,
               a4MatrixString,
               blankString
               
  ), fileConn)
  
  # Close the FD
  close(fileConn)
  
  
  return(0)
}

############ cubico como cuadratico
findThetascubic_Otro <- function(bigAs, littleAs){
    
    # Init stuff
    myThetas      = c(NaN,NaN,NaN)
    matrixLittleA = 0
    matrixTheta0  = 0
    matrixTheta1  = 0
    matrixTheta2  = 0
    
    # First, generate the matrix with the little a values
    
    matrixLittleA = matrix(c(littleAs[1],littleAs[2], littleAs[3],littleAs[4],
                             littleAs[2],littleAs[3], littleAs[4],littleAs[5],
                             littleAs[3],littleAs[4], littleAs[5],littleAs[6],
                             littleAs[4],littleAs[5], littleAs[6],littleAs[7]),4,4, byrow = TRUE)

    # Check out that we have indeed a Cramer system
    cramerValue = det(matrixLittleA)
    
    # Correct for very tiny numbers
    # if(abs(cramerValue)<0.001){ #TODO: 0.001 as it says in ??????
    #     
    #     cramerValue = 0
    # }
    
    if(cramerValue !=0){
        
        # Make the matrix for eath theta
        matrixTheta0 = matrix(c(bigAs[1],littleAs[2], littleAs[3],littleAs[4],
                                bigAs[2],littleAs[3], littleAs[4],littleAs[5],
                                bigAs[3],littleAs[4], littleAs[5],littleAs[6],
                                bigAs[4],littleAs[5], littleAs[6],littleAs[7]),4,4, byrow = TRUE)

            
        matrixTheta1 = matrix(c(littleAs[1],bigAs[1], littleAs[3],littleAs[4],
                                littleAs[2],bigAs[2], littleAs[4],littleAs[5],
                                littleAs[3],bigAs[3], littleAs[5],littleAs[6],
                                littleAs[4],bigAs[4], littleAs[6],littleAs[7]),4,4, byrow = TRUE)

            
        matrixTheta2 = matrix(c(littleAs[1],littleAs[2],bigAs[1],littleAs[4],
                               littleAs[2],littleAs[3], bigAs[2],littleAs[5],
                               littleAs[3],littleAs[4], bigAs[3],littleAs[6],
                               littleAs[4],littleAs[5], bigAs[4],littleAs[7]),4,4, byrow = TRUE)



        # Find the final thetas
        theta0 = det(matrixTheta0)/cramerValue
        theta1 = det(matrixTheta1)/cramerValue
        theta2 = det(matrixTheta2)/cramerValue
        
        # Cut the values that are smaller than 0 to 0
        if(theta0<0) theta0 = 0
        if (theta0>1) theta0 = 1
        # if(theta1<0) theta1 = 0
        # if(theta2<0) theta2 = 0
        # 
        # 
        myThetas = c(theta0, theta1, theta2)
        
    }
    
    return(myThetas)
    
}

