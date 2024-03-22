# Archivo para simulación 1
# Función coseno

# Parámetros iniciales
# B : número de réplicas del experimento
# n : tamaño de muestra
library(nlme)
library(MASS)
library(locpol)
source("toolssim.R")
#source("toolsreal.R", encoding="utf-8")
source("toolsT.R",   encoding="utf-8")



B = 500
n = 100

datos = 0
m = matrix(0,ncol=B,nrow=n)
b = matrix(0,ncol=3,nrow=B)
sigp = 0
steps = 0
hpl = 0
nh = 5
hmin = 0.01
hmax = 0.05

p.e=matrix(0,ncol=B,nrow=n)
# paquetes necesarios
for (r in 1:10) {
  datos = do.sample2(n)
  #hp = pluginBw(datos[,2],datos[,3],1,gaussK)
  zdotdot = sum(sum(datos[,3]+datos[,4]+datos[,5]))/(n*3)
  datos = cbind(datos, datos[,3]-zdotdot)
  datos = cbind(datos, datos[,4]-zdotdot)
  datos = cbind(datos, datos[,5]-zdotdot)
  datos_r = cbind(datos[,1],datos[,2],datos[,6],datos[,7],datos[,8])
  mxi = matrix(0,ncol=1,nrow=n)
  mxi = dbeta(datos[,2],2,2)
  # Doing the Backfit algorithm
  #
  # datos is the table with the data
  #
  # h 
  #
  # stopThreshold is the stop condition for when the m's are too similar to the
  #               previous step. It find the average difference, and if it is
  #               smaller, it stops.
  #
  # minSteps=2
  redh = seq(hmin,hmax,length=nh)
  
  mise = matrix(0,ncol=nh,nrow=n)
  bise = matrix(0,ncol=3,nrow=nh)
  sigpise = 0
  p.eise =  matrix(0,ncol=nh,nrow=n)
  stepsise = 0
  error = 0
  for (u in 1:nh){
    results = list()
    results = backfit.sim(300,datos_r,h= redh[u],stopThreshold = 0.0001, minSteps = 2)
    mise[,u] = results[[1]] +zdotdot# <- ( AquÃƒÂ­ estÃƒÂ¡n los resultados )
    bise[u,] = results[[2]]
    p.eise[,u] = results[[5]]
    sigpise[u] = results[[3]]
    stepsise[u] = results[[4]]
    
    error[u] =  (sum(mise[,u]-mxi)^2)/n
    
  }
  k = which.min(error)
  
  
  
  m[,r] = mise[,k] # <- ( AquÃƒÂ­ estÃƒÂ¡n los resultados )
  b[r,] = bise[k,]
  p.e[,r] = p.eise[,k]
  sigp[r] = sigpise[k]
  steps[r] = stepsise[k]
  hpl[r] = redh[k]
}


