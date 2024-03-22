# For a given raw data from a SiZer maps, tells you the percentage for each color in each SiZer map
rm(list=ls())
#setwd('C:/Users/Usuario/Dropbox/TTT final')
# Add the needed libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(lubridate)
library(latex2exp)
library(ggpubr)
# 



source("tools.R",    encoding="utf-8")
source('ttt_0.r', encoding="utf-8")
source('ttt.r', encoding="utf-8")

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
  # totalPixelsZero<-totalPixels#-sum(sizerData$ColorCodeZero == "grey")
  
  percentageGreen  = sum(sizerData$ColorCodeZero == "olivedrab")/totalPixels    # < 0 white
  percentageLemon  = sum(sizerData$ColorCodeZero == "yellow")/totalPixels       # > 0 black
  percentageBrown  = sum(sizerData$ColorCodeZero == "camel")/totalPixels        # = 0 dark-gray
  percentageGrey0  = sum(sizerData$ColorCodeZero == "grey")/totalPixels
  
  # -- First SiZer
  #totalPixelsFirst<-totalPixels#-sum(sizerData$ColorCodeFirst == "grey")
  percentageRed    = sum(sizerData$ColorCodeFirst == "red")/totalPixels        # <0 white
  percentageBlue   = sum(sizerData$ColorCodeFirst == "blue")/totalPixels       # > 0 black
  percentagePurple = sum(sizerData$ColorCodeFirst == "purple")/totalPixels     # = 0 dark-gray
  percentageGrey1  = sum(sizerData$ColorCodeFirst == "grey")/totalPixels
  
  # -- Second SiZer
  # totalPixelsSecond<-totalPixels#-sum(sizerData$ColorCodeSecond == "grey")
  percentageOrange = sum(sizerData$ColorCodeSecond == "orange")/totalPixels     # <0 white
  percentageCyan   = sum(sizerData$ColorCodeSecond == "cyan")/totalPixels       # > 0 black
  percentageVerde  = sum(sizerData$ColorCodeSecond == "green")/totalPixels      # = 0 dark-gray
  percentageGrey2  = sum(sizerData$ColorCodeSecond == "grey")/totalPixels
  
  # Write it into the dataframe and return it
  
  mySummary[1,2]  = percentageGreen
  mySummary[2,2]  = percentageLemon
  mySummary[3,2]  = percentageBrown
  mySummary[4,2]  = percentageGrey0
  mySummary[5,2]  = percentageRed
  mySummary[6,2]  = percentageBlue
  mySummary[7,2]  = percentagePurple
  mySummary[8,2]  = percentageGrey1
  mySummary[9,2]  = percentageOrange
  mySummary[10,2] = percentageCyan
  mySummary[11,2] = percentageVerde
  mySummary[12,2] = percentageGrey2
  
  return(mySummary)
  
}

pvalue.H0exp<-function(N,M,nboot,hMin,xgrid,myMethod)
{
  
  statist.boot<-matrix(NA,1,5)

  fil.res<-paste('SiZer_pval_N_',N,'_nboot_',nboot,'_hMin_',hMin,'_method_',myMethod,'_xgrid_',xgrid,'.txt',sep="")
  i<-1
  while(i<(M+1))
  {ss<-i+1;set.seed(ss)
    H0.dat<-rexp(N);#xgrid=xgrid;
    H0Sizer<-ttt(myData=H0.dat,bootstrapSample = nboot,hMin=hMin,xgrid=xgrid,
                                           myMethod=myMethod,savePlots = FALSE, saveCSV = FALSE, saveLog = FALSE)
 
  H0colors<-summaryColors(H0Sizer)
 # res2<-c(i,H0colors[9:10,2]) ## SiZer-2
  res<-c(i,H0colors[1:2,2],H0colors[9:10,2]) ##
  statist.boot<-rbind(statist.boot,res)  
  
  if (i==1){write(c('i','olivedrab','yellow','orange','cyan'),file=fil.res,append=F,ncol=length(res))}
  write(res,file=fil.res,append=T,ncol=length(res))
  
  i<-i+1
  }
 
  return(summary(statist.boot))
}
# 
 
 res<-pvalue.H0exp(N=43,M=1000,nboot=1000,hMin=0.1,xgrid=21,myMethod='quadratic')
