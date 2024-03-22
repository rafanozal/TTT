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
  
  # Return the matrix
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


