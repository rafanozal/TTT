library("nlme")
library("locpol")
#source("toolsBasic.R", encoding="utf-8")
source("toolsM.R", encoding="utf-8")

{
  
  BACKFIT_FILENAME   = "lognoise_d05g05.txt"
  DATA_FOLDER        = file.path(paste(getwd(),"/../data/", sep = ""))
  BACKFIT_FILEPATH   = file.path(paste(DATA_FOLDER, "/",BACKFIT_FILENAME,   sep = ""))
  RESULT_FOLDER      = file.path(paste(getwd(),"/../out/", sep = ""))
  
  backfitDF  = read.table(BACKFIT_FILEPATH,   fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  
}
nh=11

range = max(backfitDF[,2]) - min(backfitDF[,2])
n = length(backfitDF[,2])
hmin = range /(2*(n-1))
hmax = range/2
hValues   = 10^seq(log10(hmin),log10(hmax),length = nh) 
n= 89
# Run the model
# -----------------------------------
{
  
  # Doing some tests here to see that everything works as intended
  # Check that the y.. is zero
  
  hp= pluginBw(backfitDF[,2],backfitDF[,3],1,gaussK)
  hp=hp*2
  zdotdot = sum(sum(backfitDF[,3]+backfitDF[,4]+backfitDF[,5]))/(n*3)
  backfitDF = cbind(backfitDF, backfitDF[,3]-zdotdot)
  backfitDF = cbind(backfitDF, backfitDF[,4]-zdotdot)
  backfitDF = cbind(backfitDF, backfitDF[,5]-zdotdot)
  loc = 1:n
  datos_r = cbind(loc,backfitDF[,2],backfitDF[,6],backfitDF[,7],backfitDF[,8])
  # Doing the Backfit algorithm
  #
  # 40: Number of steps, in this case 40 steps.
  #
  # backfitDF is the table with the data
  #
  # 0.2 is the h, For h = 1 or h = 0.5 it seems to diverge
  #
  # stopThreshold is the stop condition for when the m's are too similar to the
  #               previous step. It find the average difference, and if it is
  #               smaller, it stops.
  #
  # minSteps if you want to force at least some steps.
  
  results = backfitting(200, datos_r, h=hp, stopThreshold = 0.0001, minSteps = 2)
  
  mDF = results[[1]] +zdotdot# <- ( AquÃ­ estÃ¡n los resultados )
  bDF = results[[2]]
  pDF = results[[3]]
  
}

win.graph() # Abrimos el primer dispositivo
layout(matrix(1:3,nrow=1,byrow=TRUE))
# Primera fila
plot(backfitDF[,2],backfitDF[,5],ylim=c(-26,-18.5),xlab="",ylab="")
resultado1 = cbind(backfitDF[,2],mDF[,95])
resultado1_o = resultado1[order(resultado1[,1]),]
lines(resultado1_o[,1],resultado1_o[,2]+bDF[95,3],col="black")

plot(backfitDF[,2],backfitDF[,4],ylim= c(-26,-18.5),xlab="",ylab="")
resultado2 = cbind(backfitDF[,2],mDF[,95])
resultado2_o = resultado2[order(resultado2[,1]),]
lines(resultado2_o[,1],resultado2_o[,2]+bDF[95,2],col="black")

plot(backfitDF[,2],backfitDF[,3],ylim=c(-26,-18.5),xlab="",ylab="")
resultado3 = cbind(backfitDF[,2],mDF[,95])
resultado3_o = resultado3[order(resultado3[,1]),]
lines(resultado3_o[,1],resultado3_o[,2]+bDF[95,1],col="black")
bDF[95,]
sd(pDF[,95])
#Tercera fila 
{
  
  BACKFIT_FILENAME   = "lognoise_d1g05.txt"
  DATA_FOLDER        = file.path(paste(getwd(),"/../data/", sep = ""))
  BACKFIT_FILEPATH   = file.path(paste(DATA_FOLDER, "/",BACKFIT_FILENAME,   sep = ""))
  RESULT_FOLDER      = file.path(paste(getwd(),"/../out/", sep = ""))
  
  backfitDF  = read.table(BACKFIT_FILEPATH,   fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  
}
n= 89
# Run the model
# -----------------------------------
{
  
  # Doing some tests here to see that everything works as intended
  # Check that the y.. is zero
  
  hp= pluginBw(backfitDF[,2],backfitDF[,3],1,gaussK)
  hp=hp*5
  zdotdot = sum(sum(backfitDF[,3]+backfitDF[,4]+backfitDF[,5]))/(n*3)
  backfitDF = cbind(backfitDF, backfitDF[,3]-zdotdot)
  backfitDF = cbind(backfitDF, backfitDF[,4]-zdotdot)
  backfitDF = cbind(backfitDF, backfitDF[,5]-zdotdot)
  loc = 1:n
  datos_r = cbind(loc,backfitDF[,2],backfitDF[,6],backfitDF[,7],backfitDF[,8])
  # Doing the Backfit algorithm
  #
  # 40: Number of steps, in this case 40 steps.
  #
  # backfitDF is the table with the data
  #
  # 0.2 is the h, For h = 1 or h = 0.5 it seems to diverge
  #
  # stopThreshold is the stop condition for when the m's are too similar to the
  #               previous step. It find the average difference, and if it is
  #               smaller, it stops.
  #
  # minSteps if you want to force at least some steps.
  
  results = backfitting(200, datos_r, h=hp, stopThreshold = 0.0001, minSteps = 2)
  
  mDF = results[[1]] +zdotdot# <- ( AquÃ­ estÃ¡n los resultados )
  bDF = results[[2]]
  pDF = results[[3]]
  
}

win.graph() # Abrimos el primer dispositivo
layout(matrix(1:3,nrow=1,byrow=TRUE))
# Primera fila
plot(backfitDF[,2],backfitDF[,5],ylim=c(-29,-21),xlab="",ylab="")
resultado1 = cbind(backfitDF[,2],mDF[,103])
resultado1_o = resultado1[order(resultado1[,1]),]
lines(resultado1_o[,1],resultado1_o[,2]+bDF[103,3],col="black")

plot(backfitDF[,2],backfitDF[,4],ylim=c(-29,-21),xlab="",ylab="")
resultado2 = cbind(backfitDF[,2],mDF[,103])
resultado2_o = resultado2[order(resultado2[,1]),]
lines(resultado2_o[,1],resultado2_o[,2]+bDF[103,2],col="black")

plot(backfitDF[,2],backfitDF[,3],ylim=c(-29,-21),xlab="",ylab="")
resultado3 = cbind(backfitDF[,2],mDF[,103])
resultado3_o = resultado3[order(resultado3[,1]),]
lines(resultado3_o[,1],resultado3_o[,2]+bDF[103,1],col="black")
bDF[98,]
sd(pDF[,98])
{
  
  BACKFIT_FILENAME   = "lognoise_d05g1.txt"
  DATA_FOLDER        = file.path(paste(getwd(),"/../data/", sep = ""))
  BACKFIT_FILEPATH   = file.path(paste(DATA_FOLDER, "/",BACKFIT_FILENAME,   sep = ""))
  RESULT_FOLDER      = file.path(paste(getwd(),"/../out/", sep = ""))
  
  backfitDF  = read.table(BACKFIT_FILEPATH,   fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  
}
n= 89
# Run the model
# -----------------------------------
{
  
  # Doing some tests here to see that everything works as intended
  # Check that the y.. is zero
  
  hp= pluginBw(backfitDF[,2],backfitDF[,3],1,gaussK)
  hp=hp
  zdotdot = sum(sum(backfitDF[,3]+backfitDF[,4]+backfitDF[,5]))/(n*3)
  backfitDF = cbind(backfitDF, backfitDF[,3]-zdotdot)
  backfitDF = cbind(backfitDF, backfitDF[,4]-zdotdot)
  backfitDF = cbind(backfitDF, backfitDF[,5]-zdotdot)
  loc = 1:n
  datos_r = cbind(loc,backfitDF[,2],backfitDF[,6],backfitDF[,7],backfitDF[,8])
  # Doing the Backfit algorithm
  #
  # 40: Number of steps, in this case 40 steps.
  #
  # backfitDF is the table with the data
  #
  # 0.2 is the h, For h = 1 or h = 0.5 it seems to diverge
  #
  # stopThreshold is the stop condition for when the m's are too similar to the
  #               previous step. It find the average difference, and if it is
  #               smaller, it stops.
  #
  # minSteps if you want to force at least some steps.
  
  results = backfitting(200, datos_r, h=hp, stopThreshold = 0.0001, minSteps = 2)
  
  mDF = results[[1]] +zdotdot# <- ( AquÃ­ estÃ¡n los resultados )
  bDF = results[[2]]
  pDF = results[[3]]
  
}


#Segunda fila 
win.graph() # Abrimos el primer dispositivo
layout(matrix(1:3,nrow=1,byrow=TRUE))
plot(backfitDF[,2],backfitDF[,5],ylim=c(-29,-21),xlab="",ylab="")
resultado1 = cbind(backfitDF[,2],mDF[,104])
resultado1_o = resultado1[order(resultado1[,1]),]
lines(resultado1_o[,1],resultado1_o[,2]+bDF[104,3],col="black")

plot(backfitDF[,2],backfitDF[,4],ylim=c(-29,-21),xlab="",ylab="")
resultado2 = cbind(backfitDF[,2],mDF[,104])
resultado2_o = resultado2[order(resultado2[,1]),]
lines(resultado2_o[,1],resultado2_o[,2]+bDF[104,2],col="black")

plot(backfitDF[,2],backfitDF[,3],ylim=c(-29,-21),xlab="",ylab="")
resultado3 = cbind(backfitDF[,2],mDF[,104])
resultado3_o = resultado3[order(resultado3[,1]),]
lines(resultado3_o[,1],resultado3_o[,2]+bDF[104,1],col="black")

bDF[106,]
sd(pDF[,106])

# Cuarta fila 
{
  
  BACKFIT_FILENAME   = "lognoise_d1g1.txt"
  DATA_FOLDER        = file.path(paste(getwd(),"/../data/", sep = ""))
  BACKFIT_FILEPATH   = file.path(paste(DATA_FOLDER, "/",BACKFIT_FILENAME,   sep = ""))
  RESULT_FOLDER      = file.path(paste(getwd(),"/../out/", sep = ""))
  
  backfitDF  = read.table(BACKFIT_FILEPATH,   fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  
}
n= 89
# Run the model
# -----------------------------------
{
  
  # Doing some tests here to see that everything works as intended
  # Check that the y.. is zero
  
  hp= pluginBw(backfitDF[,2],backfitDF[,3],1,gaussK)
  hp=hp*3
  zdotdot = sum(sum(backfitDF[,3]+backfitDF[,4]+backfitDF[,5]))/(n*3)
  backfitDF = cbind(backfitDF, backfitDF[,3]-zdotdot)
  backfitDF = cbind(backfitDF, backfitDF[,4]-zdotdot)
  backfitDF = cbind(backfitDF, backfitDF[,5]-zdotdot)
  loc = 1:n
  datos_r = cbind(loc,backfitDF[,2],backfitDF[,6],backfitDF[,7],backfitDF[,8])
  # Doing the Backfit algorithm
  #
  # 40: Number of steps, in this case 40 steps.
  #
  # backfitDF is the table with the data
  #
  # 0.2 is the h, For h = 1 or h = 0.5 it seems to diverge
  #
  # stopThreshold is the stop condition for when the m's are too similar to the
  #               previous step. It find the average difference, and if it is
  #               smaller, it stops.
  #
  # minSteps if you want to force at least some steps.
  
  results = backfitting(200, datos_r, h=hp, stopThreshold = 0.0001, minSteps = 2)
  
  mDF = results[[1]] +zdotdot# <- ( AquÃ­ estÃ¡n los resultados )
  bDF = results[[2]]
  pDF = results[[3]]
  
}
win.graph() # Abrimos el primer dispositivo
layout(matrix(1:3,nrow=1,byrow=TRUE))
# Primera fila
plot(backfitDF[,2],backfitDF[,5],ylim=c(-30,-21),xlab="",ylab="")
resultado1 = cbind(backfitDF[,2],mDF[,105])
resultado1_o = resultado1[order(resultado1[,1]),]
lines(resultado1_o[,1],resultado1_o[,2]+bDF[105,3],col="black")

plot(backfitDF[,2],backfitDF[,4],ylim=c(-30,-21),xlab="",ylab="")
resultado2 = cbind(backfitDF[,2],mDF[,105])
resultado2_o = resultado2[order(resultado2[,1]),]
lines(resultado2_o[,1],resultado2_o[,2]+bDF[105,2],col="black")

plot(backfitDF[,2],backfitDF[,3],ylim=c(-30,-21),xlab="",ylab="")
resultado3 = cbind(backfitDF[,2],mDF[,105])
resultado3_o = resultado3[order(resultado3[,1]),]
lines(resultado3_o[,1],resultado3_o[,2]+bDF[105,1],col="black")
bDF[106,]
sd(pDF[,106])
