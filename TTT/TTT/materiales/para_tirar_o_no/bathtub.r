



#1. Modelo Weibull aditivo (3 componentes)
#1.1. Función de azar
hw<-function(t,param){a<-param[1];b<-param[2];c<-param[3]
return(a*t^(a-1)+b*t^(b-1)+c*t^(c-1))}

#1.2. Generador aleatorio

r.w<-function(n,param)
{a<-param[1];b<-param[2];c<-param[3]
t1<-rweibull(n,shape=a,scale=1);t2<-rweibull(n,shape=b,scale=1);
t3<-rweibull(n,shape=c,scale=1)

t<-sapply(1:n,function(i){return(min(t1[i],t2[i],t3[i]))})
return(t)

}


#1.3 Prueba este caso: param<-c(3,2,0.5)

param<-c(3,2,0.5)
hw.0<-function(t)return(hw(t,param))
curve(hw.0(x),0,1.5,ylim=c(0,10))

#TTT
rw.0<-function(n){return(r.w(n,param))}

t<-rw.0(100)
mean(t); #alrededor de 0.4
library(AdequacyModel)
TTT(t)


