set.seed(123)

n        = 100
mySample = rweibull(n,2,1)
mySample = sort(mySample)



# Return an array of betas
myFunction <- function(n,j){
  
  betaValue1 = pbeta(q=(1:n)/n,     shape1 = j, shape2 = n-j+1)
  betaValue2 = pbeta(q=(0:(n-1))/n, shape1 = j, shape2 = n-j+1)
 
  return(betaValue1 - betaValue2) 
  
}


# Return a single Beta
myFunction <- function(index,n){
  
  currentinverseValue  = index/n
  previousInverseValue = (index-1)/n
  
  betaValue1 = pbeta(q= inverseValue,         shape1 = inverseValue, shape2 = n - inverseValue + 1)
  betaValue2 = pbeta(q= previousInverseValue, shape1 = inverseValue, shape2 = n - inverseValue + 1)
  
  return(betaValue1 - betaValue2) 
  
}



inverseList = (1:n)/n # 1/n , 2/n ... n/n


myMatrix = matrix(NA,n,n)

for (i in 1:n) {
  
  myValues = myFunction(n,i)
  
  for (j in 1:n) {
    
    myMatrix[i,j]  = myValues[j]
    
  }
  
}


w.r<-matrix(NA,n,n)
w.r<-sapply(1:n, myFunction())
w.r<-t(w.r)

##para la covarianza

nCrs<-function(n,r,s)
{
  return(factorial(n)/(factorial(r-1)*factorial(s-r-1)*factorial(n-s)))
}


f.rs<-function(n,r,s,ur,us) 
{
  return(nCrs(n,r,s)*ur^(r-1)*(us-ur)^(s-r-1)*(1-us)^(n-s))
}