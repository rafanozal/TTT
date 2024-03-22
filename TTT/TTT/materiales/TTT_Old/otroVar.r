generateVarianceXVectorBOOTSRAP_otro <- function(x, bootFactor = 100)
  {
  n<-length(x);M<-bootFactor
  X.boot<-matrix(NA,n,M)
  Sigma.mat<-matrix(0,n,n)
  for(j in 1:M) {
    X.boot[,j]<-sort(sample(x=x,size=n,replace=T))
  }
 
  for(i in 1:(n-1)){for(j in (i+1):n) {Sigma.mat[i,j]<-cov(X.boot[i,],X.boot[j,])}}
#  Sigma.aux<-t(Sigma.mat);
#  Sigma.mat<-Sigma.aux+Sigma.mat
  
  diag(Sigma.mat)<-sapply(1:n,function(i){var(X.boot[i,])})
   return(Sigma.mat)                       
}