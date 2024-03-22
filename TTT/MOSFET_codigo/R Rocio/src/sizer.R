# "Reset" the R session
# -----------------------------------------------------------------------------
{
  
  # You can restart R session via command line, but the code will be stuck here
  # then. You need to restart R manually if you want to drop all packages 
  # properly.
  
  rm(list = ls(all.names = TRUE)) # Clear all objects includes hidden objects.
  gc()                            # Free up memory and report the memory usage.
  
}
library("nlme")
library("locpol")
#source("toolsBasic.R", encoding="utf-8")
source("toolsM.R", encoding="utf-8")
#leemos los datos
{
  
  BACKFIT_FILENAME   = "lognoise_d05g05.txt"
  DATA_FOLDER        = file.path(paste(getwd(),"/../data/", sep = ""))
  BACKFIT_FILEPATH   = file.path(paste(DATA_FOLDER, "/",BACKFIT_FILENAME,   sep = ""))
  RESULT_FOLDER      = file.path(paste(getwd(),"/../out/", sep = ""))
  
  backfitDF  = read.table(BACKFIT_FILEPATH,   fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  
}

nh=11
R=100
range = max(backfitDF[,2]) - min(backfitDF[,2])
n = length(backfitDF[,2])
hmin = range /(2*(n-1))
hmax = range/2
alpha=0.05
hValues   = 10^seq(log10(hmin),log10(hmax),length = nh) 
zdotdot = sum(sum(backfitDF[,3]+backfitDF[,4]+backfitDF[,5]))/(n*3)
backfitDF = cbind(backfitDF, backfitDF[,3]-zdotdot)
backfitDF = cbind(backfitDF, backfitDF[,4]-zdotdot)
backfitDF = cbind(backfitDF, backfitDF[,5]-zdotdot)
loc = 1:n
datos_r = cbind(loc,backfitDF[,2],backfitDF[,6],backfitDF[,7],backfitDF[,8])
derivada = data.frame(matrix(NA, ncol = nh, nrow = n))
desvtipi = data.frame(matrix(NA, ncol = nh, nrow = n))

for (redh in 1:nh) {
  
    
    # Doing some tests here to see that everything works as intended
    # Check that the y.. is zero
    
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
    
    results = backfit.sizer(200, datos_r, h=hValues[redh], stopThreshold = 0.0001, minSteps = 2)
    
    dmDF = results[[1]]  #<- ( AquÃ­ estÃ¡n los resultados )
    mDF = results[[2]]
    bDF = results[[3]]
    pDF = results[[4]]
    e1 = backfitDF[,3]- mDF - pDF - bDF[,1]
    e2 = backfitDF[,4]- mDF - bDF[,2]- pDF
    e3 = backfitDF[,5]- mDF - bDF[,3]- pDF
    varder = data.frame(matrix(NA, ncol = R, nrow = n))
    for (rr in 1:R){
      ind_eb = sample(1:n,size=n,replace=TRUE)
      y1b = 0
      y2b=0
      y3b=0
      for (i in 1:n){
      y1b[i] = backfitDF[ind_eb[i],3]+e1[ind_eb[i]]
      y2b[i] = backfitDF[ind_eb[i],4]+e1[ind_eb[i]]
      y3b[i] = backfitDF[ind_eb[i],5]+e1[ind_eb[i]]
      
      }
      zdotdotb = sum(sum(y1b+y2b+y3b))/(n*3)
      
      loc = 1:n
      datos_b = cbind(loc,backfitDF[,2],y1b-zdotdotb,y2b-zdotdotb,y3b-zdotdotb)
      result = backfit.sizer(200, datos_b, h=hValues[redh], stopThreshold = 0.0001, minSteps = 2)
      varder[,rr] = result[[1]]
    }
    sddmDF=0
    for (i in 1:n){
    sddmDF[i] = sd(varder[i,])
  }
  derivada[,redh] = dmDF
  desvtipi[,redh] = sddmDF
}
myZScore = qnorm(1 - alpha/2)


CIDataLeftInterval   =derivada - myZScore*desvtipi
CIDataRightInterval  = derivada + myZScore*desvtipi
tam = dim(desvtipi)
CidataL = matrix(0,nrow=tam[1],ncol=tam[2]) 
CidataR = matrix(0,nrow=tam[1],ncol=tam[2]) 
for (redh in 1:nh){
temp1= approx(backfitDF[,2],CIDataLeftInterval[,redh],xout=seq(min(backfitDF[,2]),max(backfitDF[,2]),length=n),ties=mean) 
CidataL[,redh] = temp1$y
temp2 = approx(backfitDF[,2],CIDataRightInterval[,redh],xout=seq(min(backfitDF[,2]),max(backfitDF[,2]),length=n),ties=mean)
CidataR[,redh] = temp2$y
}

ESS = matrix(0,nrow=tam[1],ncol=tam[2]) 
myKernel0 =0
for (redh in 1:nh){
  myKernelESS = generateKernelH0Cube(backfitDF[,2],backfitDF[,2],h=hValues[redh],kernel="gaussian")
  myKernel0 = kernelHFunction(0,h=hValues[redh],kernel="gaussian")
  ESS[,redh]= colSums(myKernelESS)/myKernel0
}
CIL = matrix(0,nrow=tam[1],ncol=tam[2]) 
CIR = matrix(0,nrow=tam[1],ncol=tam[2]) 
ESSf = matrix(0,nrow=tam[1],ncol=tam[2]) 
for (redh in 1:nh ){
  xxL = cbind(backfitDF[,2],CIDataLeftInterval[,redh] )
  xxR = cbind(backfitDF[,2],CIDataRightInterval[,redh] )
  xxE = cbind(backfitDF[,2],ESS[,redh])
  xxL_o = xxL[order(backfitDF[,2]),]
  xxR_o = xxR[order(backfitDF[,2]),]
  ESS_o = xxE[order(backfitDF[,2]),]
  CIL[,redh] = xxL_o[,2]
  CIR[,redh] = xxR_o[,2]
  ESSf[,redh] = ESS_o[,2]
}
mapout = matrix(3,nrow=tam[1],ncol=tam[2]) # purple 
mapout[0 < CidataL & CidataR>0] = 1 #blue
mapout[0 > CidataR & CidataL<0] = 4 #red
mapout[ESSf<=5] = 2 #grey
n.colors=4
mycolors=c("blue","grey","purple","red")
win.graph()
image(seq(min(backfitDF[,2]),max(backfitDF[,2]),length=n),log10(hValues),mapout,xlab="",ylab="",col=mycolors,xlim=c(min(backfitDF[,2]),max(backfitDF[,2])),
      ylim=c(log10(min(hValues)),log10(max(hValues))),breaks=((1:(n.colors+1))-0.5))





mapout = matrix(3,nrow=tam[1],ncol=tam[2]) # purple 
mapout[0 < CIL] = 1 #blue
mapout[0 > CIR] = 4 #red
mapout[ESSf<=5] = 2 #grey
n.colors=4
mycolors=c("blue","grey","purple","red")
win.graph()
image(seq(min(backfitDF[,2]),max(backfitDF[,2]),length=n),log10(hValues),mapout,xlab="",ylab="",col=mycolors,xlim=c(min(backfitDF[,2]),max(backfitDF[,2])),
      ylim=c(log10(min(hValues)),log10(max(hValues))),breaks=((1:(n.colors+1))-0.5))
