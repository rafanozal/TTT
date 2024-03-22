#' Generate random data
#'
#' @param n How many numbers do you want.
#' 
#' @param dist "weibull" ------- Weibull distribution
#'                  p1 = shape
#'                  p2 = scale
#'                  
#'             "lognormal" ----- Log Normal distribution 
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
#' getRandomData(10,"weibull", 3, 1)
#' 
#' 
getRandomData <- function(n, dist, p1, p2, forceOneMean = FALSE){
  
  myFinalData = 0

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

    # This first loop goes throught all the variance vector that needs to be returned
    for (i in 1:n) {
      
        auxVector = 1:bootFactor

        # This second loops does the bootstraps

        for (j in 1:bootFactor) {
            
            xcopy       = x    
            xcopy       = sample( xcopy , n , replace=TRUE )
            xcopySorted = sort(xcopy)
          
            auxVector[j] = xcopySorted[i]
        
        }

        myVarianceVector[i] = var( auxVector )
        
    }

    return(myVarianceVector)
    
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
# TODO: Fulldoc
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

findVariancePhiWithKernel <- function(h,p0, kernelBarCube, sigmaMatrix, weightMatrix, order){
  
  
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

# Find the variance for a given derivate with a given kernelCube for the cubic method
#
# TODO: Fulldoc
# h is an index
# p0 is an index

findVariancePhiWithKernelCubic <- function(h,p0, kernelBarCube, sigmaMatrix, weightMatrix, order){
  
  
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