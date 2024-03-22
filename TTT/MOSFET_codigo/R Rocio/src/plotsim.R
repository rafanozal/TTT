#gráficos



dev.off()
win.graph() # Abrimos el primer dispositivo

plot(Vth,mxi-1.91,type='l',xlab="x",ylab="")
lines(Vth,mg+mean(b[,1]),col="red")


plot(Vth,mxi-0.32,type='l',ylim = c(-2.1,2.5),xlab="x",ylab="")
points(Vth,mg+mean(b[,2]),col="red")

plot(Vth,mxi+2.23,type='l',ylim = c(-2.1,2.5),xlab="x",ylab="")
points(Vth,mg+mean(b[,3]),col="red")

