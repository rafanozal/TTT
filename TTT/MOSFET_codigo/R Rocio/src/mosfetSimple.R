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
  this.dir = dirname(parent.frame(2)$ofile)
  setwd(this.dir)
}

# Add the needed libraries
# -----------------------------------------------------------------------------
library(reshape2) # melt
source("toolsBasic.R",    encoding="utf-8")
source("toolsMosfet.R",   encoding="utf-8")

# Load datafiles
# -----------------------------------------------------------------------------
{
  
  BACKFIT_FILENAME   = "datos_backfit.txt"
  DATA_FOLDER        = file.path(paste(getwd(),"/../data/", sep = ""))
  BACKFIT_FILEPATH   = file.path(paste(DATA_FOLDER, "/",BACKFIT_FILENAME,   sep = ""))
  RESULT_FOLDER      = file.path(paste(getwd(),"/../out/", sep = ""))
  
  backfitDF  = read.table(BACKFIT_FILEPATH,   fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  
}

# Transform data
#
# -- Bias condition:
#        Set proper Vd and Vg values
#        For each row, is set like this in the file "vd05vg05"
#        Transform to two new columns called Vd and Vg 
#
# -- DIE coordinates:
#        Set proper Y X coordinates values
#        For each row, is set like this in the file "0_0"
#        Transform to two new columns called x and y
# -----------------------------------------------------------------------------
{
  
  # Get the total rows
  totalRows = nrow(backfitDF)

  # Init all values to 0
  # -- Vd / Vg
  backfitDF$Vd = 0
  backfitDF$Vg = 0      
  # -- X / Y coordinates
  backfitDF$x = 0
  backfitDF$y = 0
  
  # Get the bias value for each row in that particular DF
  for (j in 1:totalRows) {
    
    # Get the raw string values
    currentLocation = backfitDF[j,1]
    currentBias     = backfitDF[j,6]
    
    # Transform the strings into numerical values
    coorValues = getXYValues(currentLocation)
    biasValues = getBiasValues(currentBias)
    
    # Set the values into the dataframe
    backfitDF$x[j]  = coorValues[[1]][1]
    backfitDF$y[j]  = coorValues[[2]][1]
    
    backfitDF$Vd[j] = biasValues[[1]][1]
    backfitDF$Vg[j] = biasValues[[2]][1]
  }
  
  # Drop the bias and location string variables once you are finish
  backfitDF$DIE_Location    = NULL
  backfitDF$Bias_conditions = NULL
  
}
     
# Run the model
# -----------------------------------
{
  
  # Doing some tests here to see that everything works as intended
  # Check that the y.. is zero
  {
    yDF = getYIJTable(backfitDF)
    print( " Result for Y.. (it should be zero) ")
    print( paste0("Y.. :", sum(yDF))  )
  }
  
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
  results = backfitting(40, backfitDF, 0.2, stopThreshold = 0.00001, minSteps = 5)
  
  mDF = results[[1]] # <- ( Aquí están los resultados )
  bDF = results[[2]]
  
}
