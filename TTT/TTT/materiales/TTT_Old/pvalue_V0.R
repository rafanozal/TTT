# For a given raw data from a SiZer maps, tells you the percentage for each color in each SiZer map
rm(list=ls())
# Add the needed libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(lubridate)
library(latex2exp)
library(ggpubr)

source("toolsV9.R",    encoding="utf-8")
source("otroVar.r",    encoding="utf-8")
source('tttV9_otro.r', encoding="utf-8")


summaryColors <- function(sizerData){
  
  
  # Find the total of pixel for that map
  
  # totalPs = unique(sizerData$p0)
  # totalHs = unique(sizerData$h)
  
  totalPixels = nrow(sizerData)
  
  # Create the dataframe for each color and each percentage.
  # There are 4 colors in each of the three maps
  
  colorsZero   = c("yellow", "olivedrab", "camel",  "grey")
  colorsFirst  = c("red",    "blue",      "purple", "grey")
  colorsSecond = c("orange", "cyan",      "green",  "grey")
  
  totalColors  = c(colorsZero, colorsFirst, colorsSecond)
  
  mySummary = expand.grid(color = totalColors, percentage = 0)
  
  # Calculate the percentage of each color
  
  # -- Zero SiZer
 totalPixelsZero<-totalPixels-sum(sizerData$ColorCodeZero == "grey")
  percentageGreen  = sum(sizerData$ColorCodeZero == "yellow")/totalPixelsZero
  percentageLemon  = sum(sizerData$ColorCodeZero == "olivedrab")/totalPixelsZero
  percentageBrown  = sum(sizerData$ColorCodeZero == "camel")/totalPixelsZero
 # percentageGrey0  = sum(sizerData$ColorCodeZero == "grey")/totalPixels
  
  # -- First SiZer
  totalPixelsFirst<-totalPixels-sum(sizerData$ColorCodeFirst == "grey")
  percentageRed    = sum(sizerData$ColorCodeFirst == "red")/totalPixelsFirst
  percentageBlue   = sum(sizerData$ColorCodeFirst == "blue")/totalPixelsFirst
  percentagePurple = sum(sizerData$ColorCodeFirst == "purple")/totalPixelsFirst
#  percentageGrey1  = sum(sizerData$ColorCodeFirst == "grey")/totalPixels
  
  # -- Second SiZer
  totalPixelsSecond<-totalPixels-sum(sizerData$ColorCodeSecond == "grey")
  percentageOrange = sum(sizerData$ColorCodeSecond == "orange")/totalPixelsSecond
  percentageCyan   = sum(sizerData$ColorCodeSecond == "cyan")/totalPixelsSecond
  percentageVerde  = sum(sizerData$ColorCodeSecond == "green")/totalPixelsSecond
 # percentageGrey2  = sum(sizerData$ColorCodeSecond == "grey")/totalPixels
  
  # Write it into the dataframe and return it
  
  mySummary[1,2]  = percentageGreen
  mySummary[2,2]  = percentageLemon
  mySummary[3,2]  = percentageBrown
  mySummary[4,2]  = NA #percentageGrey0
  mySummary[5,2]  = percentageRed
  mySummary[6,2]  = percentageBlue
  mySummary[7,2]  = percentagePurple
  mySummary[8,2]  = NA #percentageGrey1
  mySummary[9,2]  = percentageOrange
  mySummary[10,2] = percentageCyan
  mySummary[11,2] = percentageVerde
  mySummary[12,2] = NA #percentageGrey2
  
  return(mySummary)
  
}

pvalue.H0exp<-function(myData,M)
{
  N<-length(myData)
  if(N>100){xgrid=401}else{xgrid=N}
  
  res.myData<-ttt(myData=myData,xgrid=xgrid,saveCSV=FALSE)
  colors<-summaryColors(res.myData)
  cyan.myData<-(colors[10,2]) ####how many pixels are not green
  orange.myData<-(colors[9,2])
#  M<-1000 #number of bootstrap samples
  statist.boot<-matrix(NA,1,2)
  i<-1
  while(i<M)
  {H0.dat<-rexp(N);xgrid=xgrid;H0Sizer<-ttt(myData=H0.dat,xgrid=xgrid,saveCSV=FALSE)
  H0colors<-summaryColors(H0Sizer)
  statist.boot<-rbind(statist.boot,c(H0colors[9:10,2]))  
  i<-i+1
  }
  statist.boot<-statist.boot[-1,]
  colnames(statist.boot)<-c('orange','cyan')
  #pvalue.twosided<-sum(statist.boot>statist.myData)/M
  pvalue.orange<-sum(statist.boot[,1]>orange.myData)/M
  pvalue.cyan<-sum(statist.boot[,2]>cyan.myData)/M
  
  return(list(orange=pvalue.orange,cyan=pvalue.cyan,
              statist.boot=statist.boot,statist.myData=c(orange.myData,cyan.myData)))
}
# 
# granu<-read.table('granulocytic1.txt',h=T)[,1]
# aarset<-read.table('aarset.txt',h=T)[,1]
# N<-length(aarset);xgrid<-N;p.aarset<-pvalue.H0exp(myData=aarset,M=1000)
# 
# N<-length(granu);xgrid<-N;p.granu<-pvalue.H0exp(myData=granu,M=1000)
# datos<-rexp(1000);N<-1000;xgrid<-401
# pruebaExp<-pvalue.H0exp(datos,M=500)

####################
####################  Type I error
####################

typeI.error<-function(H0='Exp',N,M=1000)
{
  #N= 100, 500, 1000
  typeI.0<-double(M);typeI.2<-double(M)
  
 for(i in 1:M)
 {
  myData<-rexp(N)
  
  if(N==100){xgrid=51}else{xgrid=401}
  
  res.myData<-ttt(myData=myData,xgrid=xgrid,saveCSV=FALSE)
  
  colors<-summaryColors(res.myData)
 typeI.0[i]<-sum(colors[1:2,2])  ##number of pixels not brown (gray) in Sizer-0
 typeI.2[i]<-sum(colors[9:10,2]) ##number of pixels not green (gray) in Sizer-2

 }
 return(list(typeI.0,typeI.2))
  
}

#############################
###### type II error and test power

typeII.error<-function(H1=model,N,M=1000)
{
  #N= 100, 500, 1000
 typeII.2<-double(M)
 if (model=='IFR'){simula<-function(N)return(rgamma(N,shape=5,scale=1/5))}
 if (model=='BFR'){simula<-function(N)
                           {data<-double(N)
                           for(i in 1:N)
                          {data[i]<-min(c(rweibull(1,shape=3,scale=2.5),rweibull(1,shape=2,scale=2.5),rweibull(1,shape=.5,scale=2.5)))}
                           return(data)
                           }}
 
  if (model=='UFR'){simula<-function(N){return(rlnorm(N,meanlog=-0.4,sdlog=1))}}
 
  for(i in 1:M)
  {
     myData<-simula(N)
    
    if(N==100){xgrid=51}else{xgrid=401}
    
    res.myData<-ttt(myData=myData,xgrid=xgrid,saveCSV=FALSE)
    
    colors<-summaryColors(res.myData)
   # typeII.0[i]<-sum(colors[3,2])  ##number of pixels brown (gray) in Sizer-0 #NO TIENE SENTIDO EN ERROR TYPE II
    typeII.2[i]<-sum(colors[11,2]) ##number of pixels green (gray) in Sizer-2
    
  }
  return(list(typeII.2))
  
}

