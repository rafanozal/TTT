# For a given raw data from a SiZer maps, tells you the percentage for each color in each SiZer map


hypothesisColors <- function(sizerData){
  
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
  
  percentageGreen  = sum(sizerData$ColorCodeZero == "yellow")/totalPixels
  percentageLemon  = sum(sizerData$ColorCodeZero == "olivedrab")/totalPixels
  percentageBrown  = sum(sizerData$ColorCodeZero == "camel")/totalPixels
  percentageGrey0  = sum(sizerData$ColorCodeZero == "grey")/totalPixels
  
  # -- First SiZer
  
  percentageRed    = sum(sizerData$ColorCodeFirst == "red")/totalPixels
  percentageBlue   = sum(sizerData$ColorCodeFirst == "blue")/totalPixels
  percentagePurple = sum(sizerData$ColorCodeFirst == "purple")/totalPixels
  percentageGrey1  = sum(sizerData$ColorCodeFirst == "grey")/totalPixels
  
  # -- Second SiZer
  
  percentageOrange = sum(sizerData$ColorCodeSecond == "orange")/totalPixels
  percentageCyan   = sum(sizerData$ColorCodeSecond == "cyan")/totalPixels
  percentageVerde  = sum(sizerData$ColorCodeSecond == "green")/totalPixels
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