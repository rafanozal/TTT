
num <- Sys.getenv("SGE_TASK_ID")
num<-as.integer(num)

#rm(list=ls())
source('colors_grid.R')
source("tools.R",    encoding="utf-8")

source('test_ttt.R', encoding="utf-8")

K<-1:5
models<-c(1,2,3)
param<-expand.grid(K,models)
# 
# # 
tarea<-param[num,1]
model<-param[num,2]

N<-50;M<-200


fil.res<-paste('Power_H1_',model,'_tarea',tarea,'_N_',N,'.txt',sep="")

###funcion para simular desde un modelo
if (model==1){simula<-function(N)return(rgamma(N,shape=5,scale=1/5))}
if (model==2){simula<-function(N)
             {data<-double(N)
              for(i in 1:N)
              {data[i]<-min(c(rweibull(1,shape=3,scale=2.5),rweibull(1,shape=2,scale=2.5),rweibull(1,shape=.5,scale=2.5)))}
              return(data)}
             }

if (model==3){simula<-function(N){return(rlnorm(N,meanlog=-0.5,sdlog=1))}}
###1=IFR; 2=BFR;3=UFR
#####################################

ini<-(tarea-1)*M+1;fin<-ini+(M-1)
for(i in ini:fin)
{
  ss<-i;set.seed(ss)
  myData<-simula(N)
  
 # if(N==100){xgrid=51}else{xgrid=401}
xgrid<-21
  res.myData<-test.ttt(myData=myData,xgrid=xgrid,hMin=0.1,bootstrapSample = 250)
  
  
  colors<-summaryColors(res.myData)
  
  # no.typeII.0[i]<-sum(colors[3,2])  ##proportion of  not brown (gray) pixels in Sizer-0 #NO TIENE SENTIDO EN ERROR TYPE II
  no.typeII.2<-sum(colors[9:10,2]) ##proportion of not green (gray) pixels in Sizer-2
  res<-c(i,no.typeII.2)
  
  if (i==ini){write(c('i','test.power'),file=fil.res,append=F,ncol=length(res))}
  write(res,file=fil.res,append=T,ncol=length(res))
  
}
  
  

