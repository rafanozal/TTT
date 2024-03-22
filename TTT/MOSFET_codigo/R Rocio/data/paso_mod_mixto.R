newdata<-read.table('newdata_backfit.txt',h=T)
head(newdata)
# DIE_Location   Vth   F10e4Hz    F10e3Hz  F10e2Hz
# 1          0_0 0.422 -2.472500 -1.2291675 1.817817
# 2          0_1 0.414 -1.915262 -0.1220609 2.012245
# 3          0_2 0.390 -2.399097 -0.3145988 1.773732
# 4          0_6 0.398 -1.715417 -0.4498027 2.404191
# 5         0_m1 0.438 -2.243464 -0.4109625 1.726670
# 6         0_m2 0.394 -2.091202 -0.3848370 1.891361

### tenemos que hacer esto en cada iteración del algoritmo en el que ajustamos
### el modelo mixto
### y.ij[r]=y.ij-b.j[r]-pi.i[r]
### Paso 0: m0(x)=0
### empezamos con los datos originales:

y.ij<-c(newdata[,3],newdata[,4],newdata[,5]) ### lognoise en formato largo
x.i<-rep(newdata[,2],3) ### la covariable  (Vth)   ### este no se usa en el paso del modelo mixto!!!!
freq.j<-c(rep(1e2,89),rep(1e3,89),rep(1e4,89))  ### el factor 
subject.i<-rep(newdata[,1],3)  ### la localizacion del individuo son los grupos:89 en total
dat.mix<-data.frame(subject.i,freq.j,y.ij,x.i)

## para ajustar el modelo mixto y predecir los pi.i, usamos la funcion lme del paquete nlme
library(nlme)
fit.0<-lme(y.ij~factor(freq.j),random=~1|subject.i)

pred<-predict(fit.0,asList=T,level=0:1)
pi.i.1<-pred$predict.subject-pred$predict.fixed  ### esta es la prediccion de pi.i con el modelo mixto
ngrupos<-89 #numero de sujetos
b.j.1<-pred$predict.fixed[(0:2)*ngrupos+1]  ## la suma vale 0 (efectos fijos)

## Paso 1: para el suavizado:
y.ij.1<-y.ij-pred$predict.fixed ## esto es lo mismo que y.ij.1<-y.ij-b.j.1-pi.i.1


## ahora habría que suavizar: y.ij.1=m(xi)
### SUM_i{SUM_j{(y.ij.1-beta0-beta1(x-xi))^2 K_h(x-xi)}}
## como el suavizado sólo se hace en la dimension de i, no en la de j, esto es como suavizar los promedios de y.ij.1
## en la dimensión i, es decir que calculamos:
y.ij.1<-matrix(y.ij.1,ngrupos,3) ## matriz 89 x 3. Cada fila es un sujeto.
y.i.1<-rowMeans(y.ij.1)
## suavizamos:
## SUM_i{y.i.1-beta0-beta1(x-xi))^2 K_h(x-xi)}

## Paso 2. modelo mixto (otra vez) ahora con los datos centrados segun la estimacion anterior:
y.ij.2<-y.ij-m.1(xi)



