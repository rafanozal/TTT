# "Reset" the R session
# -----------------------------------------------------------------------------
{
  
  # You can restart R session via command line, but the code will be stuck here
  # then. You need to restart R manually if you want to drop all packages 
  # properly.
  
  rm(list = ls(all.names = TRUE)) # Clear all objects includes hidden objects.
  gc()                            # Free up memory and report the memory usage.
  
}


# Set working directory to file location
# To do this manually, go to:
# Session -> Set Working Directory -> To Source File Location
# -----------------------------------------------------------------------------
{
  this.dir = "C:/R/src"
  setwd(this.dir)
}

# Add the needed libraries
# -----------------------------------------------------------------------------
library(reshape2) # melt
library("nlme")
library("locpol")
#source("toolsBasic.R", encoding="utf-8")
source("toolsM.R", encoding="utf-8")

# Load datafiles
# -----------------------------------------------------------------------------
{
  
  BACKFIT_FILENAME   = "newdata_g4_backfit.txt"
  DATA_FOLDER        = file.path(paste(getwd(),"/../data/", sep = ""))
  BACKFIT_FILEPATH   = file.path(paste(DATA_FOLDER, "/",BACKFIT_FILENAME,   sep = ""))
  RESULT_FOLDER      = file.path(paste(getwd(),"/../out/", sep = ""))
  
  backfitDF  = read.table(BACKFIT_FILEPATH,   fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  
}

# Run the model
# -----------------------------------
{
  
  # Doing some tests here to see that everything works as intended
  # Check that the y.. is zero
  {
    yDF = backfitDF[3:5]
    print( " Result for Y.. (it should be zero) ")
    print( paste0("Y.. :", sum(yDF))  )
  }
  hp= pluginBw(backfitDF[,2],backfitDF[,3],1,gaussK)
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
  
  results = backfitting(200, backfitDF, h=hp, stopThreshold = 0.0001, minSteps = 2)
  
  mDF = results[[1]] # <- ( AquÃ­ estÃ¡n los resultados )
  bDF = results[[2]]
  pDF = results[[3]]
  
}


# Datos no centrados

{
  
  BACKFIT_FILENAME   = "g4_nocentrado.txt"
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

