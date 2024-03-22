sigma<-1;mu<--(sigma)^2/2
t<-rlnorm(100,meanlog=mu,sdlog=sigma)
mean(t)
TTT(t)

h.ln<-function(x)
{return(dlnorm(x,meanlog=mu,sdlog=sigma)/(1-plnorm(x,meanlog=mu,sdlog=sigma)))}
qlnorm(0.999,meanlog=mu,sdlog=sigma)
curve(h.ln(x),0,12)
