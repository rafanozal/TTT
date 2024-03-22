# Stirling<- function (n)((2*pi*n)^.5)*(n/exp(1))^n
# 
# f<-function(n,x){return(Stirling(n)/(Stirling(x)*Stirling(n-x)))}
# 
# n<-1e2
# f.n<-function(x)f(n,x)

n<-100
ss<-1;set.seed(ss)

r<-80;s<-90
i<-1;x.r<-c();x.s<-c()
while(i<=10000)
{
xi<-rweibull(n,shape=2,scale=1)
oxi<-sort(xi); 
x.r<-c(x.r,oxi[r]);x.s<-c(x.s,oxi[s])
i<-i+1}

n<-100
ss<-1;set.seed(ss)
xi<-rweibull(n,shape=2,scale=1)
oxi<-sort(xi); 


x.r.boot<-c();x.s.boot<-c()
i<-1; R<-10000
while(i<=R)
{
  xi.boot<-sample(xi,n,replace=T)
  oxi.boot<-sort(xi.boot); 
  x.r.boot<-c(x.r.boot,oxi.boot[r])
  x.s.boot<-c(x.s.boot,oxi.boot[s])
  i<-i+1
}

cov(x.r,x.s);cov(x.r.boot,x.s.boot)
