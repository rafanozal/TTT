# For a given raw data from a SiZer maps, tells you the percentage for each color in each SiZer map
#rm(list=ls())
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
#source('ttt.r', encoding="utf-8")
source('test_ttt.r', encoding="utf-8")

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

