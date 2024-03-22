# Add the needed libraries
# ---- Basics
library(ggplot2)       # Basic ggplot2 library
library(ggnewscale)    # Allows for multiple color scales in one plot
library(RColorBrewer)  # Color settins and palettes
library(shadowtext)    # Drop shadows in text for better viewing
library(latex2exp)     # Allows latex expressions in ggplot
library(qqplotr)       # QQ plots with bands
library(forcats)       # Reverse order of factors


library(grid)
library(gridExtra)     # grid.arrange()


# Get the extension of a filePath
#
# (string) filePath is any form of path as in "../asfaf/asdfas/asf/myFile.txt"
#
# return extension including the dot, if exist
#
# ie:
# 
# "../asfaf/asdfas/asf/myFile.txt"   -> ".txt"
# "../asfaf/asdfas/asf/myFile"       -> ""
# "../asfaf/asdfas/asf/myFile."      -> "."
# "../asfaf/asdfas/asf/my.Fi.le.txt" -> ".txt"
#
getFileExtension <- function(filePath){
  
  myExtension = ""
  
  texFilePathSplitted       = strsplit(filePath, "/")[[1]] # Get rid of the filePath and get only the fileName
  texRelativeLocationPath   = texFilePathSplitted[length(texFilePathSplitted)]
  
  texFilePathExtension      = strsplit(texRelativeLocationPath, "\\.")[[1]]
  
  if(length(texFilePathExtension) > 1){
    
    myExtension = paste(".", texFilePathExtension[length(texFilePathExtension)], sep='')
    
  }
  
  return(myExtension)
  
}

# Check if a string is a filename or a folder name
#
# Simply check that it ends with a .xxxxx extension or not
# Right now it only checks for PNGs since is the only image format that matters
# (aside from vector images, which doesn't work in ggplot2)
#
# Return
# 
#     TRUE  = This is a filepath of a folder (/home/user/doc)
#     FALSE = This is a filepath of a file   (/home/user/doc/myImage.png)
checkIfFolder <- function (myFileString){
  
  result = TRUE
  
  # If the file string is NULL, we don't have a file
  if(is.null(myFileString)){
    
    result = FALSE
    
  }
  else{
    
    # If the file string doesn't ends in PNG, we don't have a file
    if( grepl("\\.png$", myFileString) == FALSE ){
      
      result = FALSE
    }
    
  }
  
  return (result)
  
}


# Get parameters for the theme depending on your theme
# - white
# - simple  / same as white without legend
# - blank   / without anything
# - default / bad naming, this means default ggplot2 theme
getThemeParameters <- function(themeName){
  
  # Init the default theme
  myTheme = list(rep(NA,10))
  
  myTheme[[1]] = element_blank()                  # 1: Background color
  myTheme[[2]] = element_line(colour = "black")   # 2: Axys X and Y Line color
  myTheme[[3]] = element_line(colour = "grey70")  # 3: Axys Y , major breaks
  myTheme[[4]] = element_blank()                  # 4: Axys X , major breaks
  myTheme[[5]] = "right"                          # 5: Legend
  
  
  # If you didn't gave a NULL theme
  if(!is.null(themeName)){
    
    if(themeName == "white"){ # White theme is the default too
      myTheme[[1]] = element_blank()
      myTheme[[2]] = element_line(colour = "black")
      myTheme[[3]] = element_line(colour = "grey70")
      myTheme[[4]] = element_blank()
      myTheme[[5]] = "right"
    }
    
    if(themeName == "blank"){
      myTheme[[1]] = element_blank()
      myTheme[[2]] = element_blank()
      myTheme[[3]] = element_blank()
      myTheme[[4]] = element_blank()
      myTheme[[5]] = "none"
    }
    
    if(themeName == "default"){
      myTheme[[1]] = element_rect(colour = "grey")   # 1: Background color
      myTheme[[2]] = element_blank()                 # 2: Axys X and Y Line color
      myTheme[[3]] = element_line(colour = "white")  # 3: Axys Y , major breaks
      myTheme[[4]] = element_line(colour = "white")  # 4: Axys X , major breaks
      myTheme[[5]] = "right"                         # 5: Legend
    }
    
    if(themeName == "simple"){
      myTheme[[1]] = element_blank()
      myTheme[[2]] = element_line(colour = "black")
      myTheme[[3]] = element_line(colour = "grey70")
      myTheme[[4]] = element_blank()
      myTheme[[5]] = "none"
    }
    
  }
  
  return(myTheme)
  
}



# Gives an automatic filepath to save your file.
#
# This is very useful when you want to generate filepath for doing a lot of
# plots, analysis or whatever, and give consistent naming rules.
#
# -- If your filePath is a file, do nothing and return the same filepath.
#
# -- If your filePath is a folder, it make up a filepath for you inside that
#    folder with these rules:
#
#     <File_Type> + <Table Name> + <List of names from the relevant variables> +
#     <extension type>
#
#     Example:
#
#     "/../tables/BMIPlot_myPatientsDataFrame_Sex_BMI.png"
# 
#    "/../tables/LatexTable_myPatientsDataFrame.tex"
#
# String filePath       ; Either a filePath with a file or a folder (read above)
#
# DataFrame myDataFrame ; (NULL) The data frame that you want to use for
#                                automatic name giving
#
# String tableName      ; (NULL) The name of myDataFrame. You can't get this
#                                automatically with deparse inside the function.
#                                so you need to call deparse before the function
#                                if you want the same name as the variable, or,
#                                you can override the name of the dataframe by
#                                giving whatever name you want here.
#
# String fileType       ; (NULL) What kind of file you are trying to generate.
#                                It recommended that you don't leave this to
#                                NULL as default. The options are:
#
#                                NULL = You have no idea of what you are
#                                       generating, but the function will try
#                                       to generate a meaningful name anyway
#                                       with prefix "Unknown_File_Type"
#
#                                images = You are making an image, probably a
#                                          plot from the plotting library.
#
#                                             "AbsBarplot"
#                                             "RelBarplot"
#                                             "LongAbsBarplot"
#                                             "BMIPlot"
#                                             "Boxplot"
#                                             "RechabilityBoxplot"
#                                             "Histogram"
#                                             "Density"
#                                             "CategoricalHistogram"
#                                             "CategoricalDensity"
#                                             "pValuesHeatmap"
#                                             "QQ"
#                                             "Scatterplot"
#                                             "Tableplot"
#                                             "SimulaitonLineplot"
#
#                                latex = You are making a latex file. It
#                                          could be a file that contain an
#                                          image or a table, but a .tex file
#                                          nevertheless.
#
#                                             "LatexTable"
#                                             "LatexImage"
#
#                                tables = A .txt table. Probably you run
#                                         already all the relevant functions
#                                         to make this txt human friendly.
#
#                                             "txtTable"
#
#
#                                Finally, if you didn't gave NULL, but you gave
#                                an option that is not in the list of option,
#                                you will get a "Invalid_File_Type" prefix.

# Int variableIndexX    ; (NULL) The column index of myDataFrame that you want
#                                to use for the automatic name giving. You can
#                                have up to 3 variables depending of the type
#                                of plot or file that you are making

# bool rootPath ; (FALSE) ignore the extension and return the filePath with no
#                         extension. Usefull to generate automatic names for
#                         the same analysis where you are going to have a .png,
#                         .tex, and so on.

automaticFilePath <- function(filePath, rootPath = FALSE,
                              myDataFrame = NULL, tableName = NULL,
                              fileType = NULL, variableIndex1 = NULL,
                              variableIndex2 = NULL, variableIndex3 = NULL){
  
  # This is the final string variable to return
  finalPath = filePath
  
  # If we have a complete filepath, do nothing.
  # Otherwise, gives an automatic name accordingly
  haveFile = checkIfFolder(filePath)
  if(haveFile == FALSE){
    
    # Defaults
    fileExtension = ""
    
    # If you didn't go for NULL filetype, check if the filetype you have is a
    # valid one and give the proper file extension.
    if(!is.null(fileType)){
      
      # (these list are sorted alphabetically)
      
      validImages = c("AbsBarplot",
                      "BMIPlot",
                      "Boxplot",
                      "CategoricalBoxplot",
                      "CategoricalHistogram",
                      "CategoricalDensity",
                      "CombinedRelBarplot",
                      "Density",
                      "Histogram",
                      "LongAbsBarplot",
                      "pValuesHeatmap",
                      "QQ",
                      "RechabilityBoxplot",
                      "RelBarplot",
                      "SimulationLinePlot",
                      "Scatterplot",
                      "Tableplot",
                      "Waffle")
      
      validAnalysis = c("XiTable")
      
      validLatex  = c("LatexImage",
                      "LatexTable")
      
      validTables = c("txtTable")
      
      if(fileType %in% validImages)   fileExtension = ".png"
      if(fileType %in% validLatex)    fileExtension = ".tex"
      if(fileType %in% validTables)   fileExtension = ".txt"
      if(fileType %in% validAnalysis) fileExtension = ".tex"
      
      # If you fail all the check, you made a human error
      if(fileExtension == "") fileType = "Invalid_File_Type"
      
    }
    else{
      
      fileType      = "Unknown_File_Type" # Gives you a very generic plot name
      # if you don't provide one
    }
    
    # The name of your dataframe variable is given as a parameter
    #
    # Note that if you do this:
    #
    #    tableName = deparse(substitute(myDataFrame))
    #
    # It will not copy the original dataframe name, but the name of which
    # whatever function you are calling this function. So you need to pass the
    # tableName manually. Also note that R doesn't need to pass variables as
    # reference for performance. (or so they say so at CRAN)
    #
    # In any case, give an automatic name if you don't have one
    if(is.null(tableName)){
      
      tableName = "NoTableNameGiven"
      
    }
    
    
    # Name of each of the variables
    name1 = colnames(myDataFrame)[variableIndex1] # It doesn't matter if it is NULL, you get an empty string
    name2 = colnames(myDataFrame)[variableIndex2]
    name3 = colnames(myDataFrame)[variableIndex3]
    
    # Depending of your file type, you return one type of string or another
    # In any case, put everything together in this variable
    fileName = ""
    
    # Unknown and Invalid
    if(fileType == "Unknown_File_Type")    fileName = paste(fileType, "_", tableName, "_", name1, "_", name2, "_", name3, sep="")
    if(fileType == "Invalid_File_Type")    fileName = paste(fileType, "_", tableName, "_", name1, "_", name2, "_", name3, sep="")
    # TXT
    # Latex
    if(fileType == "LatexTable")           fileName = paste(fileType, "_", tableName,                                     sep="")
    # Images
    #     0 Variables (Specials)
    if(fileType == "RechabilityBoxplot")   fileName = paste(fileType, "_", tableName,                                     sep="")
    if(fileType == "TablePlot")            fileName = paste(fileType, "_", tableName,                                     sep="")
    if(fileType == "SimulationLinePlot")   fileName = paste(fileType, "_", tableName,                                     sep="")
    #     1 Variable
    #         Categoricals
    if(fileType == "AbsBarplot")           fileName = paste(fileType, "_", tableName, "_", name1,                         sep="")
    if(fileType == "LongAbsBarplot")       fileName = paste(fileType, "_", tableName, "_", name1,                         sep="")
    if(fileType == "RelBarplot")           fileName = paste(fileType, "_", tableName, "_", name1,                         sep="")
    #         Numericals
    if(fileType == "Boxplot")              fileName = paste(fileType, "_", tableName, "_", name1,                         sep="")
    if(fileType == "Density")              fileName = paste(fileType, "_", tableName, "_", name1,                         sep="")
    if(fileType == "Histogram")            fileName = paste(fileType, "_", tableName, "_", name1,                         sep="")
    if(fileType == "QQ")                   fileName = paste(fileType, "_", tableName, "_", name1,                         sep="")
    #     2 Variables
    #         Categorical + Categorical
    if(fileType == "CombinedRelBarplot")   fileName = paste(fileType, "_", tableName, "_", name1, "_", name2,             sep="")
    if(fileType == "XiTable")              fileName = paste(fileType, "_", tableName, "_", name1, "_", name2,             sep="")
    #         Categorical + Numerical
    if(fileType == "BMIPlot")              fileName = paste(fileType, "_", tableName, "_", name1, "_", name2,             sep="")
    if(fileType == "CategoricalBoxplot")   fileName = paste(fileType, "_", tableName, "_", name1, "_", name2,             sep="")
    if(fileType == "CategoricalHistogram") fileName = paste(fileType, "_", tableName, "_", name1, "_", name2,             sep="")
    if(fileType == "CategoricalDensity")   fileName = paste(fileType, "_", tableName, "_", name1, "_", name2,             sep="")
    if(fileType == "pValuesHeatmap")       fileName = paste(fileType, "_", tableName, "_", name1, "_", name2,             sep="")
    #         Numerical   + Numerical
    if(fileType == "Scatterplot")          fileName = paste(fileType, "_", tableName, "_", name1, "_", name2,             sep="")
    #     3 Variables
    #         Numerical + Numerical + Categorical
    if(fileType == "CombineRegression")    fileName = paste(fileType, "_", tableName, "_", name1, "_", name2, "_", name3, sep="")
    if(fileType == "Waffle")               fileName = paste(fileType, "_", tableName, "_", name1, "_", name2, "_", name3, sep="")
    
    # Merge all the information into a variable and return it
    # -- If you have an actual filepath, use that (TODO: FIX THIS!! This has nothing to do with roots)
    if(rootPath == TRUE){
      finalPath = paste(filePath, fileName, sep="")
    }
    # -- If you have a folder
    else{
      
      # If you have one of the special analysis where you generate several
      # tables and images, each with it own name
      if(fileType == "XiTable"){
        
        finalPath      = vector("list", length = 5)
        finalPath[[1]] = paste0(filePath, fileName, "_summary",    fileExtension)    
        finalPath[[2]] = paste0(filePath, fileName, "_absolute",   fileExtension)    
        finalPath[[3]] = paste0(filePath, fileName, "_relative",   fileExtension)    
        finalPath[[4]] = paste0(filePath, fileName, "_difference", fileExtension)    
        finalPath[[5]] = paste0(filePath, fileName, "_pvalues",    fileExtension)    
        
      }
      # In any other case, return what you need
      else{
        finalPath = paste(filePath, fileName, fileExtension, sep="")    
      }
    }
  }
  
  return(finalPath)
  
}



