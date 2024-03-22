library(AdequacyModel)

sigma<-1;mu<--(sigma)^2/2
t<-rlnorm(100,meanlog=mu,sdlog=sigma)
mean(t)
par(mfrow=c(1,2))
TTT(t)

h.ln<-function(x)
{return(dlnorm(x,meanlog=mu,sdlog=sigma)/(1-plnorm(x,meanlog=mu,sdlog=sigma)))}
qlnorm(0.999,meanlog=mu,sdlog=sigma)
curve(h.ln(x),0,12)



###bathtub
r.bt<-function(n,param)
{
sc<-2.5;a<-param[1];b<-param[2];c<-param[3]
t<-double(n)
for(i in 1:n)
{
t1<-rweibull(1,shape=a,scale=sc)
t2<-rweibull(1,shape=b,scale=sc)
t3<-rweibull(1,shape=c,scale=sc)
t[i]<-min(t1,t2,t3)
}
return(t)
}
param<-c(3,2,0.5)

rbt.0<-function(n){return(r.bt(n,param=param))}
times<-rbt.0(100)
mean(times)

TTT(times)
