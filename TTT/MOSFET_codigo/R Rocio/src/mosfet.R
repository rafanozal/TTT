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
source("toolsPlotting.R", encoding="utf-8")
source("constants.R",     encoding="utf-8")


# Load datafiles
# -----------------------------------------------------------------------------
{
  logNoiseDF = read.table(LOG_NOISE_FILEPATH, fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  noNoiseDF  = read.table(NO_NOISE_FILEPATH,  fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  noNullDF   = read.table(NO_NULL_FILEPATH,   fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  largoDF    = read.table(LARGO_FILEPATH,     fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  backfitDF  = read.table(BACKFIT_FILEPATH,   fileEncoding = "UTF-8", stringsAsFactors = FALSE)

}

# Transform data
# -----------------------------------------------------------------------------
{

  # Largo DF is on a melted data format, we need to make the DF to coincide
  # with the rest of the DF columns format
  {
    
    # Count how many rows we are suppose to have
    F100RowsDF   = largoDF[largoDF$Freq == 100,]
    totalNewRows = nrow(F100RowsDF)
    
    # Create the dataframe were we are going to write everything
    newLargoDF = data.frame(matrix(NA, nrow = totalNewRows, ncol = 6))
    colnames(newLargoDF) = colnames(logNoiseDF)
    
    # The original file has way too much data redundancy, the Vth, DIE and bias
    # coincide for each frequency, so we just copy those in the same order
    newLargoDF$Vth             = F100RowsDF$Vth
    newLargoDF$DIE_Location    = F100RowsDF$DIE
    newLargoDF$Bias_conditions = F100RowsDF$biasCond
    
    # Now we just need to copy the Frequencies
    # TODO: Check that the log noise is the actual value that is suppose to go here
    newLargoDF$F10e2Hz = F100RowsDF$lognoise
    newLargoDF$F10e3Hz = largoDF[largoDF$Freq == 1000 ,2]
    newLargoDF$F10e4Hz = largoDF[largoDF$Freq == 10000,2]
    
  }

  # Add all the datafarames into a list
  {
    allDFList      = vector("list", length = 5)
    allDFList[[1]] = logNoiseDF
    allDFList[[2]] = noNoiseDF
    allDFList[[3]] = noNullDF
    allDFList[[4]] = newLargoDF
    allDFList[[5]] = backfitDF
  }

  # Save the number of rows for each dataframe
  {
    totalDFs        = length(allDFList)  
    totalRowsVector = rep(0,totalDFs)
    
    for(i in 1:totalDFs){
      
      totalRowsVector[i] = nrow(allDFList[[i]])
      
    }
    
  }

  # Set proper Vd and Vg values
  # Set proper Y X coordinates values
  # Set the noise value in Formula 2.2
  {
    
    # For each DF
    for(i in 1:totalDFs){
      
      currentDF        = allDFList[[i]]  
      currentTotalRows = totalRowsVector[i]

      # Init all values to 0
      # -- Vd / Vg
      currentDF$Vd = 0
      currentDF$Vg = 0      
      # -- X / Y coordinates
      currentDF$x = 0
      currentDF$y = 0
      
      # Get the bias value for each row in that particular DF
      for (j in 1:currentTotalRows) {
        
        # Get the raw string values
        currentLocation = currentDF[j,1]
        currentBias     = currentDF[j,6]
        
        # Transform the strings into numerical values
        coorValues = getXYValues(currentLocation)
        biasValues = getBiasValues(currentBias)
        
        # Set the values into the dataframe
        currentDF$x[j]  = coorValues[[1]][1]
        currentDF$y[j]  = coorValues[[2]][1]
        
        currentDF$Vd[j] = biasValues[[1]][1]
        currentDF$Vg[j] = biasValues[[2]][1]
      }
      
      # Drop the bias and location string variables once you are finish
      currentDF$DIE_Location    = NULL
      currentDF$Bias_conditions = NULL
      
      # Find the noise value, Formula in 2.2
      # TODO: Wrong stimate, correct
      myEstimate      = 0
      currentDF$Noise = currentDF$F10e2Hz + currentDF$F10e3Hz + currentDF$F10e4Hz + currentDF$Vd + currentDF$Vg + myEstimate
      
      # Change the columns order, the final order is this:
      # Vth, F10e2, F10e3, F10e4, Noise, Vd, Vg, x, y
      currentDF = currentDF[,c(1,4,3,2,9,5,6,7,8)]
      
      # Write the DF back to the list
      # R is a horrible language that doesn't allow for pointers
      allDFList[[i]] = currentDF
      
    }
  }

}

# R is a horrible language that doesn't allow for pointers. Set the DF back again
# This is not needed, but makes debugging easy
{
  logNoiseDF = allDFList[[1]]
  noNoiseDF  = allDFList[[2]]
  noNullDF   = allDFList[[3]]
  newLargoDF = allDFList[[4]]
  backfitDF  = allDFList[[5]]
}

# Create a simple DF to test the results for debugging
{
  simpleDF = newLargoDF
  simpleDF = simpleDF[c(1,2,3,4),]
  simpleDF[,1] = c( 1, 2, 3, 4)
  simpleDF[,2] = c(-2,-1, 0, 1)
  simpleDF[,3] = c(-1, 0, 1, 2)
  simpleDF[,4] = c( 0, 1, 2, 3)
  simpleDF[,5] = c( 0.1, 0.2, 0.3, 0.4)
}

# Select a DF to run the backfit algorithm
{
  selectedDF        = backfitDF
  #selectedDF        = newLargoDF
  #selectedDF        = simpleDF
  selectedTotalRows = ncol(selectedDF)
  variablesNames    = colnames(selectedDF)
}

# Getting Indexes
# -----------------------------------------------------------------------------
{
  # Get the indexes of the variables
  # You don't need to do this, but makes the code easy to read later
  xCoordinateIndex = grep("^x$",       colnames(selectedDF))
  yCoordinateIndex = grep("^y$",       colnames(selectedDF))
  noiseIndex       = grep("^Noise$",   colnames(selectedDF))
  VthIndex         = grep("^Vth$",     colnames(selectedDF))
  f102Index        = grep("^F10e2Hz$", colnames(selectedDF))
  f103Index        = grep("^F10e3Hz$", colnames(selectedDF))
  f104Index        = grep("^F10e4Hz$", colnames(selectedDF))
  
}

  
# Doing plots
#    -- Check for normality of all variables (nothing is normal)
#    -- Show the waffle plot
# -----------------------------------------------------------------------------
{

  # QQPlots for testing normality
  # -- There is nothing normal here
  print("QQ")
  for(i in 1:selectedTotalRows){
  
    doQQPlot(selectedDF, i , RESULT_FOLDER, plotTitle = variablesNames[i],
             plotSubtitle="(Rounded Shapiro-Wilk’s test, p-value > 0.05 => NORMALITY)")
    
  }
  
  # Scatterplot for noise
  
  # Waffle plot
  print("Waffle")
  doWafflePlot(selectedDF, VthIndex, xCoordinateIndex, yCoordinateIndex, RESULT_FOLDER,
              plotTitle = "Waffle with the noise values")
  
}


# Run the model
# -----------------------------------------
{
  
  # Doing some tests here to see that everything works as intended
  # Check that the y.. is zero
  {
    yDF = getYIJTable(selectedDF)
    print( paste0("Y.. :", sum(yDF))  )
  }
  
  # Doing the final step
  # 40: Number of steps, in this case 40 steps.
  # selectDF which is the table with the data that you want
  # 0.2 is the h, For h = 1 or h = 0.5 it seems to diverge
  results = backfitting(40, selectedDF, 0.2)
  
  mDF = results[[1]] # <- ( Aquí están los resultados )
  bDF = results[[2]]
  
}


