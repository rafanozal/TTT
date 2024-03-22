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



pvalue.H0exp<-function(myData,M){
  
  N<-length(myData)

  res.myData<-ttt(myData=myData,xgrid=xgrid,saveCSV=FALSE)
  
  colors<-summaryColors(res.myData)
  
  cyan.myData<-(colors[10,2]) ####how many pixels are not green
  
  orange.myData<-(colors[9,2])
  
#  M<-1000 #number of bootstrap samples
  
  statist.boot<-matrix(NA,1,2)
  
  
  for (i in 1:M) {
    
    
    
  }
  
  i<-1
  while(i<M){
    
    H0.dat<-rexp(N);
    xgrid=xgrid;
    H0Sizer<-ttt(myData=H0.dat,xgrid=xgrid,saveCSV=FALSE)
    H0colors<-summaryColors(H0Sizer)
    statist.boot<-rbind(statist.boot,c(H0colors[9:10,2]))  
    i<-i+1
    
  }
  
  statist.boot<-statist.boot[-1,]
  colnames(statist.boot)<-c('orange','cyan')
  
  #pvalue.twosided<-sum(statist.boot>statist.myData)/M
  
  pvalue.orange<-sum(statist.boot[,1]>orange.myData)/M
  pvalue.cyan<-sum(statist.boot[,2]>cyan.myData)/M
  
  toReturn = list(orange=pvalue.orange, cyan=pvalue.cyan,
                  statist.boot=statist.boot, statist.myData=c(orange.myData,cyan.myData))
  
  return()
}

# 
# granu<-read.table('granulocytic1.txt',h=T)[,1]
# aarset<-read.table('aarset.txt',h=T)[,1]
# N<-length(aarset);xgrid<-N;p.aarset<-pvalue.H0exp(myData=aarset,M=1000)
# 
# N<-length(granu);xgrid<-N;p.granu<-pvalue.H0exp(myData=granu,M=1000)
# datos<-rexp(1000);N<-1000;xgrid<-401
# pruebaExp<-pvalue.H0exp(datos,M=500)



