

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

