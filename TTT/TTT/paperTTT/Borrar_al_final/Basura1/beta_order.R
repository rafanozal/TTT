
rm(list=ls(all=TRUE))
#mu_r^(m)=sum_{i=1}^n w_(i,r) * X_r^m  #momento de orden m

#w_(i,r): es una matriz cuya entrada (i,r) es
#        w_(i,r)=F(x=i/n,a=r,b=n-r+1)-F(x=(i-1)/n,a=r,b=n-r+1)
#donde F(x,a,b) es la función de distribución de una v.a. Beta con parámetros a y b
#la clave está en partir la integral de 0 a 1 en subintegrales con límites en los puntos de la red ((i-1)/n, i/n)
#pero no se aproximan estas subintegrales sino que junto con el coeficiente Ci (del documento anterior)
#dan lugar al área de la densidad beta entre dos puntos consecutivos de la red de {i/n}
#Así nosotros no calculamos los números factoriales, lo hace R 

#Esto es una prueba sencilla con una Weibull
n<-100
ss<-1;set.seed(ss)
xi<-rweibull(n,2,1)
xxi<-sort(xi)



#############################################################
###Hutson and Ernst (2000)

###para la media y la varianza

j<-(1:n)/n
w.r<-matrix(NA,n,n)

w.r<-sapply(1:n,function(j)return(pbeta(q=(1:n)/n,shape1=j,shape2=n-j+1)-
                                     pbeta(q=(0:(n-1))/n,shape1=j,shape2=n-j+1)))
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

# ###### ###esta funcion "f.rs" es "ddirichlet" del paquete MCMCpack,
# da la densidad, no el área que es lo que necesitamos

# n<-100;r<-20;s<-30
# us<-0.3;ur<-0.1
# f.rs(n,r,s,ur,us)
# library(MCMCpack)
# ddirichlet(c(ur,us-ur,1-us),c(r,s-r,n-s+1))
  #####


###calculo los pesos del artículo de Hutson y Ernst (2000)

w.rs<-matrix(0,n,n) #r<s

v.rs<-double(n) #r<s

library(pracma) #para integral doble numérica
#calculo como caso particular: cov(X.20,X.30)
#hay que hacerlo para todo r<s

r=80
s=90

#defino la función del integrando

fun<-function(x,y){return(f.rs(n=n,r=r,s=s,ur=y,us=x))}
#x<-us;y<-ur

#hay que hacer una doble integral
for(j in 2:n) #estos ciclos for hay que sustituirlos por algo mejor
 {
  for(i in 1:(j-1))
  {
    w.rs[i,j]<-integral2(fun, xmin=(j-1)/n,xmax=j/n,ymin=(i-1)/n,ymax=i/n, reltol = 1e-10)$Q
  }
   ymax<-function(x){x}
   v.rs[j]<-integral2(fun, xmin=(j-1)/n,xmax=j/n,ymin=(j-1)/n,ymax=ymax, reltol = 1e-10)$Q
   rm(ymax)
}

#calculo media y varianza de X.r y X.s

 r<-80
 m.r<-sum(w.r[r,]*xxi)
 V.r<-sum(w.r[r,]*xxi^2)-(sum(w.r[r,]*xxi))^2

 s<-90
 m.s<-sum(w.r[s,]*xxi)
 V.s<-sum(w.r[s,]*xxi^2)-(sum(w.r[s,]*xxi))^2
 

XX.r<-matrix(xxi-m.r,n,1);
XX.s<-matrix(xxi-m.s,n,1)


#la covarianza se calcula como la fórmula (3.2) de la pg 92 de Hutson y Ernst

cov.rs1<-t(XX.r)%*%w.rs%*%XX.s
cov.rs2<-sum(v.rs*XX.r*XX.s)
cov.rs<-cov.rs1+cov.rs2

cov.rs






