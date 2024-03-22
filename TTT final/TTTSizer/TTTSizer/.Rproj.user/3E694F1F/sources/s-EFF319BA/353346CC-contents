#' @title TTTSizer
#'
#' @description Delete this. Function to test Roxygen2
#'
#' Very complicated equations runs this function:
#'
#' \deqn{a + b}
#'\deqn{ \sigma = \sqrt{ \frac{Z}{n} \sum
#'  \left[ \textstyle\frac{1}{2}\displaystyle
#'    \left( \log \frac{H_i}{L_i} \right)^2  - (2\log 2-1)
#'    \left( \log \frac{C_i}{O_i} \right)^2 \right] }
#'}{sqrt(N/n * runSum(0.5 * log(OHLC[,2]/OHLC[,3])^2 -
#'           (2*log(2)-1) * log(OHLC[,4]/OHLC[,1])^2, n))}
#'
#' @return A "Hello, world!" string
#'
#' @examples hello()
#'
#' @export

hello <- function() {
  print("Hello, world!")
}

#' @title TTTSizer
#'
#' @description Read a file with numbers and return them as vector
#'
#' @param filePath The file path of the file
#'
#' @return A vector with the numbers
#'
#' @examples getDataFromFile("/home/me/myNumbers.txt")
#'
#' @export

getDataFromFile <- function(filePath){

  df <- read.table(filePath, header = FALSE)

  return (df[[1]])

}

#' @title TTTSizer
#'
#' @description For a given number, return the value of the kernel function, with the specified kernel
#'
#' @param x (float) A value for the symetric kernel function
#' @param kernel The selected kernel function
#'
#'               \describe{
#'                   \item{gaussian}{exp(-(x^2) / 2) / sqrt(2 * pi)}
#'                   \item{biweight}{(abs(x)<=1)*15/16*(1-x^2)^2}
#'                   \item{triweight}{(abs(x)<=1)*35/32*(1-x^2)^3}
#'                   \item{epanechnikov}{(abs(x)<=1)*((1-x^2)*3/4)}
#'                }
#'
#'
#' @return (float) The kernel function result.
#'
#' @examples kernelFunction(1.5)
#'           kernelFunction(0, kernel="biweight")
#'
#' @export

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

#' @title TTTSizer
#'
#' @description Correct the kernel of a given data with a given badwith h and k
#'
#' @param x the data to find the kernel, can be a vector
#' @param h Any given bandwith
#' @param kernel The selected kernel function ("gaussian", "biweight", "triweight", "epanechnikov")
#'
#' @return a float vector of size length(x) with the result of the kernel
#'
#' @examples kernelHFunction(1,2,"gaussian")
#'
#' @export
kernelHFunction <- function(x,h,kernel){

  myCorrection = (1/h) * kernelFunction(x/h,kernel)

  return(myCorrection)

}

#' @title TTTSizer
#'
#' @description For a given vector of pi's and h's, generate all possible K_h(pi-p0) combinations.
#'
#' @param piVector The vector with all the frecuencies (can be a vector of size one)
#' @param p0Vector The vector with only the p0 values (can be a vector of size one)
#' @param hVector  The vector with only the h's values (can be a vector of size one)
#' @param kernel The selected kernel function ("gaussian", "biweight", "triweight", "epanechnikov")
#'
#' @return A list of matrixs with all the K_h(pi-p0) combinations pre-calculated.
#'     Where the index of the list represent the h value, the row is pi and the column p0
#'
#'     The row is the h, the column is the sum(K_h(pi-p0))
#'
#' @export

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

#' @title TTTSizer
#'
#' @description For a given vector of pi's and generate all possible (pi-p0)^n combinations.
#'     Where n=0,1,2,3,4,5,6
#'
#' @param piVector The vector with all the frecuencies (can be a vector of size one)
#' @param p0Vector The vector with only the p0 values (can be a vector of size one)
#'
#' @return A list of matrix with all the(pi-p0)^n combinations pre-calculated.
#'     The index of the list is n, the row of the matrix is pi and the column is p0
#'
#' @export
#'

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

#' @title TTTSizer
#'
#' @description Finds a single A_r(p0). This is use the function that finds all A's. (uppercase)
#' @seealso findBigAvalues
#'
#' @note The function needs you to calculate first the following constants
#'
#'     \enumerate{
#'         \item - a kernel cube with all the K_h(pi-p0) values
#'         \item - a polinomial cube with all the (pi-p0)^n values
#'         \item - a vector with all the phi(xis) values
#'     }
#'
#' @param p0           - An index number refering to the p0 value
#' @param h            - An index number refering to the bandwith
#' @param r            - index with the ^r value. (Note, R indexing start at 1, 1 = -> ^0, 2 -> ^1, and so on )
#' @param phiVector    - The vector with the phi values
#' @param kernelCube   - The cube with all the Kh(pi-p0) values
#' @param poliCube     - The cube with all the (pi-p0)^n values
#'
#' @return A single float with the Ar(p0)
#'
#' @export

findBigA <- function(p0,h,r,phiVector,kernelCube,poliCube){

  myAr = sum(phiVector*poliCube[[r]][,p0]*kernelCube[[h]][,p0])

  return(myAr)

}

#' @title TTTSizer
#'
#' @description For a given set of data, find the A matrix of dimession i x 1
#' @seealso findThetas
#'
#' @note The function needs you to calculate first the following constants
#'
#'     \enumerate{
#'         \item - a kernel cube with all the K_h(pi-p0) values
#'         \item - a polinomial cube with all the (pi-p0)^n values
#'         \item - a vector with all the phi(xis) values
#'     }
#'
#' @param p0Index      - An index number refering to the p0 value
#' @param hIndex       - An index number refering to the bandwith
#' @param phiVector    - The vector with the phi values
#' @param kernelCube   - The cube with all the Kh(pi-p0) values
#' @param poliCube     - The cube with all the (pi-p0)^n values
#'
#' @return Return a vector with the (A0,A1,A2, A3) values
#'
#' @export

findBigAvalues <- function(p0Index,hIndex,phiVector,kernelCube,poliCube){

  A0 = findBigA(p0Index,hIndex,1,phiVector,kernelCube,poliCube)
  A1 = findBigA(p0Index,hIndex,2,phiVector,kernelCube,poliCube)
  A2 = findBigA(p0Index,hIndex,3,phiVector,kernelCube,poliCube)
  A3 = findBigA(p0Index,hIndex,4,phiVector,kernelCube,poliCube)

  return(c(A0,A1,A2,A3))

}

#' @title TTTSizer
#'
#' @description Finds a single a_r(p0). This is use the function that finds all a's. (lowercase)
#' @seealso findThetas
#'
#' @note The function needs you to calculate first the following constants
#'
#'     \enumerate{
#'         \item - a kernel cube with all the K_h(pi-p0) values
#'         \item - a polinomial cube with all the (pi-p0)^n values
#'         \item - a vector with all the phi(xis) values
#'     }
#'
#' @param p0           - An index number refering to the p0 value
#' @param h            - An index number refering to the bandwith
#' @param r            - index with the ^r value. (Note, R indexing start at 1, 1 = -> ^0, 2 -> ^1, and so on )
#' @param kernelCube   - The cube with all the Kh(pi-p0) values
#' @param poliCube     - The cube with all the (pi-p0)^n values
#'
#' @return A single float with the Ar(p0)
#'
#' @export

findLittleA <- function(p0,h,r,kernelCube,poliCube){

  myAr = sum(poliCube[[r]][,p0]*kernelCube[[h]][,p0])

  return(myAr)

}


#' @title TTTSizer
#'
#' @description Find thetas for a quadratic method.
#'     For a given set of data, find the variable matrix solution of theta0, theta1, and theta2
#'
#' @param bigAs     vector with the big As associated with that p0,x,p,h and kernel
#' @param littleAs  vector with the little As associated with that p0, p,h and kernel
#'
#' @return c(NaN, NaN, NaN) - If is not a Cramer system
#'         c(theta0, theta1, theta2)    - Otherwise
#'
#' @export

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

    myThetas = c(theta0, theta1, theta2)

  }

  return(myThetas)

}

#' @title TTTSizer
#'
#' @description Find thetas for a cubic method.
#'     For a given set of data, find the variable matrix solution of theta0, theta1, and theta2
#'
#' @param bigAs     vector with the big As associated with that p0,x,p,h and kernel
#' @param littleAs  vector with the little As associated with that p0, p,h and kernel
#'
#' @return c(NaN, NaN, NaN) - If is not a Cramer system
#'         c(theta0, theta1, theta2)    - Otherwise
#'
#' @export

findThetasCubic <- function(bigAs, littleAs){

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

    myThetas = c(theta0, theta1, theta2)

  }

  return(myThetas)

}

#' @title TTTSizer
#'
#' @description For a given kernelCube and polinomialCube, and a's matrix, generate all possible Kbar_h(pi-p0) combinations.
#'              This is only for the quadratic version, thus you only need 5 matrices
#'
#' @note The function needs you to calculate first the following constants
#'
#'     \enumerate{
#'         \item - a kernel cube with all the K_h(pi-p0) values
#'         \item - a polinomial cube with all the (pi-p0)^n values
#'      }
#'
#' @param kernelCube   - The cube with all the Kh(pi-p0) values
#' @param poliCube     - The cube with all the (pi-p0)^n values
#' @param anMatrix     - A matrix with the pi,h values of the a_n's
#' @param order        - The order of the kernel cube you want to generate
#'
#' @return  A list of matrixs with all the Kbar_h(pi-p0) combinations pre-calculated.
#'          Where the index of the list represent the h value, the row is pi and the column p0
#'
#' @export

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

        denominator = a0*a2*a4 + 2*a1*a2*a3 - a2^3 - a0*a3^2 - a1^2*a4

        if(order == 0) numerator   = (a1*a3-a2^2)*p2Value + (a2*a3 - a1*a4)*p1Value + (a2*a4-a3^2)
        if(order == 1) numerator   = (a1*a2-a0*a3)*p2Value + (a0*a4 - a2^2)*p1Value + (a2*a3-a1*a4)
        if(order == 2) numerator   = (a0*a2-a1^2)*p2Value + (a1*a2 - a0*a3)*p1Value + (a1*a3-a2^2)

        # Check that we have "enought value"

        # Save it into the matrix
        matrixByHIndex[[i]][k,j] <- kernelValue * numerator/denominator

      }

    }

  }

  return(matrixByHIndex)

}

#' @title TTTSizer
#'
#' @description For a given kernelCube and polinomialCube, and a's matrix, generate all possible Kbar_h(pi-p0) combinations.
#'              This is for the cubic version, and thus you need the 7 matrices. In the paper, this is notated as K tilde
#'
#' @note The function needs you to calculate first the following constants
#'
#'     \enumerate{
#'         \item - a kernel cube with all the K_h(pi-p0) values
#'         \item - a polinomial cube with all the (pi-p0)^n values
#'      }
#'
#' @param kernelCube   - The cube with all the Kh(pi-p0) values
#' @param poliCube     - The cube with all the (pi-p0)^n values
#' @param anMatrix     - A matrix with the pi,h values of the a_n's
#' @param order        - The order of the kernel cube you want to generate
#'
#' @return A list of matrixs with all the Kbar_h(pi-p0) combinations pre-calculated.
#'         Where the index of the list represent the h value, the row is pi and the column p0
#'
#' @export

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
        denominator   = deltaComplete

        # Put everything toguether based on order and save
        if(order == 0) numerator   = delta00 - delta01*p1Value + delta02*p2Value - delta03*p3Value
        if(order == 1) numerator   = -delta10 + delta11*p1Value - delta12*p2Value + delta13*p3Value
        if(order == 2) numerator   = delta20 - delta21*p1Value + delta22*p2Value - delta23*p3Value

        # Save it into the matrix
        matrixByHIndex[[i]][k,j] <- kernelValue * numerator/denominator

      }

    }

  }

  return(matrixByHIndex)

}

#' @title TTTSizer
#'
#' @description Generate a weight matrix object of size n x n with this format
#'
#' 1    0       0      ... 0
#' 1/n (n-1)/n  0      ... 0
#' 1/n  1/n    (n-2)/n ... 0
#'          ...
#' 1/n  1/n     1/n    ... 1/n
#'
#' @param n Size of the matrix
#'
#' @return A matrix of floats as specified in the description
#'
#' @export

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

#' @title TTTSizer
#'
#' @description Generate the Variance/Coovariance matrix with a bootstrap algorithm
#'
#' @param x The original vector with your data points.
#' @param bootFactor The number of bootstrap iterations, default is 100
#'
#' @return A matrix with all the coovariances and variances. The main diagonal contains the variance vector
#'
#' @export

generateVarianceXVectorBOOTSTRAP <- function(x, bootFactor = 100){

  # Init stuff
  n                = length(x);
  mySamplingMatrix = matrix(NA,n,bootFactor) # Each column represent each of the bootstrap iterations
  myResultMatrix   = matrix(0,n,n)

  # For each of the bootstrap iteration, select a sample from the original datavector and put it into the sampling matrix
  for(j in 1:bootFactor) {
    mySamplingMatrix[,j]<-sort(sample(x=x,size=n,replace=T))
  }

  # Create the diagonal covariance matrix
  # The main diagonal is still 0 after this loop.
  for(i in 1:(n-1)){
    for(j in (i+1):n) {

      myResultMatrix[i,j] = cov(mySamplingMatrix[i,],mySamplingMatrix[j,])

    }
  }

  # Generate the main diagonal in the matrix
  # The main diagonal now contains the variance
  for(i in 1:n){

    myResultMatrix[i,i] = var(mySamplingMatrix[i,])


  }

  return(myResultMatrix)

}

#' @title TTTSizer
#'
#' @description Find the variance for a given derivate with a given kernelCube
#'
#' @param h             - An index number refering to the bandwith.
#' @param p0            - An index number refering to the p0 value.
#' @param kernelBarCube - A matrix with all the Kbar_h(pi-p0) combinations.
#' @param sigmaMatrix   - A Variance/Coovariance matrix.
#' @param weightMatrix  - A weight matrix object.
#' @param order         - The order of the kernel cube you want to generate
#'
#' @return The variance value that correspond to the given indexes
#'
#' @export

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

#' @title TTTSizer
#'
#' @description Find the variance for a given derivate with a given kernelCube (cubic version)
#'
#' @param h             - An index number refering to the bandwith.
#' @param p0            - An index number refering to the p0 value.
#' @param kernelBarCube - A matrix with all the Kbar_h(pi-p0) combinations.
#' @param sigmaMatrix   - A Variance/Coovariance matrix.
#' @param weightMatrix  - A weight matrix object.
#' @param order         - The order of the kernel cube you want to generate
#'
#' @return The variance value that correspond to the given indexes
#'
#' @export

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
  if(order != 2){

    returnValue = kernelRow %*% weightMatrix %*% sigmaMatrix %*% transposeWeightMatrix %*% kernelColumn

  }
  else{

    returnValue = 4 * kernelRow %*% weightMatrix %*% sigmaMatrix %*% transposeWeightMatrix %*% kernelColumn

  }


  return (returnValue)

}

#' @title TTTSizer
#'
#' @description After running a TTT function, summarize the results of each SiZer map into a dataframe.
#'              - The first column of the dataframe is the color name
#'              - The second column of the dataframe is proportion of pixels for that SiZer map.
#'              - Notice that rows are group 4 by 4. The first 4 is the SiZer-0, the next 4 SiZer-1, and the last 4 the SiZer-2
#'
#' @param sizerData The summary returned by the TTT function
#'
#' @return A 12 x 2 dataframe with the color info summarized.
#'
#' @export

summaryColors <- function(sizerData){


  # Find the total of pixel for that map
  totalPixels = nrow(sizerData)

  # Create the dataframe for each color and each percentage.
  # There are 4 colors in each of the three maps

  colorsZero   = c("yellow", "olivedrab", "camel",  "grey")
  colorsFirst  = c("red",    "blue",      "purple", "grey")
  colorsSecond = c("orange", "cyan",      "green",  "grey")

  totalColors  = c(colorsZero, colorsFirst, colorsSecond)

  mySummary = expand.grid(color = totalColors, percentage = 0)

  # Calculate the percentage of each color

  # -- Zero SiZer
  percentageGreen  = sum(sizerData$ColorCodeZero == "olivedrab")/totalPixels    # < 0 white
  percentageLemon  = sum(sizerData$ColorCodeZero == "yellow")/totalPixels       # > 0 black
  percentageBrown  = sum(sizerData$ColorCodeZero == "camel")/totalPixels        # = 0 dark-gray
  percentageGrey0  = sum(sizerData$ColorCodeZero == "grey")/totalPixels

  # -- First SiZer
  percentageRed    = sum(sizerData$ColorCodeFirst == "red")/totalPixels        # <0 white
  percentageBlue   = sum(sizerData$ColorCodeFirst == "blue")/totalPixels       # > 0 black
  percentagePurple = sum(sizerData$ColorCodeFirst == "purple")/totalPixels     # = 0 dark-gray
  percentageGrey1  = sum(sizerData$ColorCodeFirst == "grey")/totalPixels

  # -- Second SiZer
  percentageOrange = sum(sizerData$ColorCodeSecond == "orange")/totalPixels     # <0 white
  percentageCyan   = sum(sizerData$ColorCodeSecond == "cyan")/totalPixels       # > 0 black
  percentageVerde  = sum(sizerData$ColorCodeSecond == "green")/totalPixels      # = 0 dark-gray
  percentageGrey2  = sum(sizerData$ColorCodeSecond == "grey")/totalPixels

  # Write it into the dataframe and return it

  mySummary[1,2]  = percentageGreen
  mySummary[2,2]  = percentageLemon
  mySummary[3,2]  = percentageBrown
  mySummary[4,2]  = percentageGrey0
  mySummary[5,2]  = percentageRed
  mySummary[6,2]  = percentageBlue
  mySummary[7,2]  = percentagePurple
  mySummary[8,2]  = percentageGrey1
  mySummary[9,2]  = percentageOrange
  mySummary[10,2] = percentageCyan
  mySummary[11,2] = percentageVerde
  mySummary[12,2] = percentageGrey2

  return(mySummary)

}

