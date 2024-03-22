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

p.e=matrix(0,ncol=B,nrow=n)
# paquetes necesarios
for (r in 447:B) {
  datos = do.sample2(n)
  hp = pluginBw(datos[,2],datos[,3],1,gaussK)
  hp = hp/6
  zdotdot = sum(sum(datos[,3]+datos[,4]+datos[,5]))/(n*3)
  datos = cbind(datos, datos[,3]-zdotdot)
  datos = cbind(datos, datos[,4]-zdotdot)
  datos = cbind(datos, datos[,5]-zdotdot)
  datos_r = cbind(datos[,1],datos[,2],datos[,6],datos[,7],datos[,8])
  
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
  results = list()
  results = backfit.sim(300,datos_r,h= hp,stopThreshold = 0.0001, minSteps = 2)
  
  m[,r] = results[[1]] +zdotdot# <- ( AquÃÂ­ estÃÂ¡n los resultados )
  b[r,] = results[[2]]
  p.e[,r] = results[[5]]
  sigp[r] = results[[3]]
  steps[r] = results[[4]]
  hpl[r] = hp
}

