
num <- Sys.getenv("SGE_TASK_ID")
num<-as.integer(num)

#rm(list=ls())
source('colors_grid.R')
source("tools.R",    encoding="utf-8")

source('test_ttt.R', encoding="utf-8")

K<-1:5
Ns<-c(50,100,500,1000)
param<-expand.grid(K,Ns)
# 
# # 
tarea<-param[num,1]
N<-param[num,2]
M<-200

fil.res<-paste('Error_TypeI_N_',N,'_tarea_',tarea,'.txt',sep="")

#####################################

ini<-(tarea-1)*M+1;fin<-ini+(M-1)
for(i in ini:fin)
{
  ss<-i;set.seed(ss)
  myData<-rexp(N)
  
  if(N==50){xgrid=21; hMin=0.1; hMax=1}
  if(N==100){xgrid=51; hMin=0; hMax=1}
  if(N==500){xgrid=401; hMin=0; hMax=1}
  if(N==1000){xgrid=401; hMin=0;hMax=0.9}

  res.myData<-test.ttt(myData=myData,xgrid=xgrid,bootstrapSample = 500,hMin=hMin,hMax=hMax)
  
  
  colors<-summaryColors(res.myData)
  
  typeI.0<-sum(colors[1:2,2])  ##proportion of  not brown (dark gray) pixels in Sizer-0 
  typeI.2<-sum(colors[9:10,2]) ##proportion of not green (dark gray) pixels in Sizer-2
  res<-c(i,typeI.0,typeI.2)
  
  if (i==ini){write(c('i','TypeI_0','TypeI_2'),file=fil.res,append=F,ncol=length(res))}
  write(res,file=fil.res,append=T,ncol=length(res))
  
}
  
  

