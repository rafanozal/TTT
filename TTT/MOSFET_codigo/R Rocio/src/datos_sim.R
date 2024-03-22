# calculos sobre las simulaciones
mean(b[,1])
# Sesgo

mean(b[,1]-(-1.91))
# Varianza
sd(b[,1])

mean(b[,2])
# Sesgo
mean(b[,2]-(-0.32))
# Varianza
sd(b[,2])

mean(b[,3])
# Sesgo
mean(b[,3]-(2.23))
# Varianza
sd(b[,3])
mean(steps)
sd(steps)
min(steps)
max(steps)
mean(sigp)
mean(sigp-(0.05))
sd(sigp)

Vth = seq(0,1,length=n)
mxi = (cos(Vth*3*pi)+1)/10
mxi = dbeta(Vth,2,2)
mg=0
for (j in 1:n){
  mg[j]=mean(m[j,])
}

error = 0
for (r in 1:B){
  for (l in 1:n){
    error = error + (m[l,r]-mxi[l])^2
  }
}
error_t = error/(n*B)


# Gráfico de la simulacion 1
x11();
par(mar=c(3,3.5,3.5,1.5),oma=c(2,0.5,0.5,0.2),mgp=c(2,0.5,0),cex.axis=1,cex.lab=1.1,cex.main=1.2)

plot(Vth,mxi,type="l", col="black",xlab="",ylab="", ylim=c(0,0.22), main="Nonparametric component estimation: Model 1",)
lines(Vth,mg,type="l", lty="dashed",col = "black",xlab="",ylab="")
legend(x = "topright", bty="n",  legend = c("true curve", "estimated curve"),lty = c(1, 2), col = c("black", "black"),
       lwd = 2)
# Gráfico de la simulación 2

#n=500
x11();
par(mar=c(3,3.5,3.5,1.5),oma=c(2,0.5,0.5,0.2),mgp=c(2,0.5,0),cex.axis=1,cex.lab=1.1,cex.main=1.2)

plot(Vth,mxi,type="l", col="black",xlab="",ylab="", main="Nonparametric component estimation: Model 2",ylim=c(0,1.6))
lines(Vth,mg,type="l", lty="dashed", col = "black",xlab="",ylab="")
legend(x = "topright", bty="n",  legend = c("true curve", "estimated curve"),lty = c(1, 2), col = c("black", "black"),
       lwd = 2) 
x11();
par(mar=c(3,3.5,3.5,1.5),oma=c(2,0.5,0.5,0.2),mgp=c(2,0.5,0),cex.axis=1,cex.lab=1.1,cex.main=1.2)

plot(Vth,mxi+2.23,type='l',xlab="x",ylab="",col="black",main="Model 2, m + b1",ylim=c(2,3.85))
lines(Vth,mg+mean(b[,5]),col="black",type="l", lty="dashed",xlab="",ylab="")
legend(x = "topright", bty="n",  legend = c("true curve", "estimated curve"),lty = c(1, 2), col = c("black", "black"),
       lwd = 2) 
