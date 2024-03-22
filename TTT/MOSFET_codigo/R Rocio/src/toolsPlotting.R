
# -----------------------------------------------
#     QQ
# -----------------------------------------------
{
  
  doQQPlot <- function(tableBase, variableIndex, plotFilePath,
                       plotTitle = NULL, plotSubtitle = NULL, plotCaption = NULL, plotXLabel = NULL, plotYLabel = NULL,
                       plotTheme = NULL,
                       overrideTableName = NULL
                       ){
    
    # Init variables
    myPlotType = "QQ"
    myTableName = deparse(substitute(tableBase))
    
    # If you need to override your table name, then do it
    if(!is.null(overrideTableName)){
      myTableName = overrideTableName
      
    }
    
    
    # Get an automatic name if you don't have a proper one
    imageFilePath = ""
    
    myFileExtension = getFileExtension(plotFilePath)
    
    if(myFileExtension == ".png") imageFilePath = plotFilePath
    
    else imageFilePath = automaticFilePath(plotFilePath, myDataFrame = tableBase,
                                           tableName = myTableName, fileType = myPlotType,
                                           variableIndex1 = variableIndex)
    
    # Get the theme information
    themeData = getThemeParameters(plotTheme)
    
    # Get info about different categories and variables
    variableName   = colnames(tableBase)[variableIndex]
    
    # Run the normality test
    # (Shapiro-Wilk’s test., p-value > 0.05 => NORMALITY)
    # First we need to check if all the values are the same
    # In that case is not normally distributed
    normalValue   = 0
    areAllEqual   = sd(tableBase[,variableIndex], na.rm = TRUE)
    if(areAllEqual != 0)
      normalValue   = shapiro.test(tableBase[,variableIndex])$p.value
    
    plotLabelNormalValue = round(normalValue,2)
    
    # Do the plot
    ggplot(tableBase, aes(sample = tableBase[,variableIndex] )) +
      
      # The 95% confident bands
      geom_qq_band(bandType = "pointwise", alpha = 0.3) +
      
      # The normality line
      stat_qq_line(color="blue") +
      
      # Each point of data
      stat_qq_point() +
      
      # Add the normality text result
      geom_text(aes(x=-Inf,y=Inf,hjust=0,vjust=1, label=plotLabelNormalValue), color="red", parse = TRUE) +
      
      # Create titles and subtitles
      labs(title    = plotTitle,
           subtitle = plotSubtitle,
           caption  = plotCaption,
           color    = variableName,
           x = plotXLabel, y = plotYLabel) +
      
      # Apply the theme
      theme(panel.background   = themeData[[1]],
            axis.line          = themeData[[2]],
            panel.grid.major.y = themeData[[3]],
            panel.grid.major.x = themeData[[4]])
    
    ggsave(imageFilePath, width = 8)
    
    return(normalValue)
  }
  
}


# -----------------------------------------------
#     Waffle
# -----------------------------------------------
{
  
  doWafflePlot <- function(tableBase, variableIndex, coordinateXIndex, coordinateYIndex, plotFilePath,
                           plotTitle = NULL, plotSubtitle = NULL, plotCaption = NULL, plotXLabel = NULL, plotYLabel = NULL,
                           plotTheme = NULL,
                           overrideTableName = NULL
  ){
    
    # Init variables
    myPlotType = "Waffle"
    myTableName = deparse(substitute(tableBase))
    
    # If you need to override your table name, then do it
    if(!is.null(overrideTableName)){
      myTableName = overrideTableName
      
    }
    
    
    # Get an automatic name if you don't have a proper one
    imageFilePath = ""
    
    myFileExtension = getFileExtension(plotFilePath)
    
    if(myFileExtension == ".png") imageFilePath = plotFilePath
    
    else imageFilePath = automaticFilePath(plotFilePath, myDataFrame = tableBase,
                                           tableName = myTableName, fileType = myPlotType,
                                           variableIndex1 = variableIndex,
                                           variableIndex2 = coordinateXIndex,
                                           variableIndex3 = coordinateYIndex)
    
    # Get the theme information
    themeData = getThemeParameters(plotTheme)
    
    # Get info about different categories and variables
    variableName   = colnames(tableBase)[variableIndex]
    
    # Do the plot
    ggplot(tableBase, aes(x = tableBase[,coordinateXIndex] , y = tableBase[,coordinateYIndex])) +

      # Each dot in the waffle
      geom_point(aes(colour = tableBase[,variableIndex]), size = 3) +
      
      # Create titles and subtitles
      labs(title    = plotTitle,
           subtitle = plotSubtitle,
           caption  = plotCaption,
           color    = variableName,
           x = plotXLabel, y = plotYLabel) +
      
      # Apply the theme
      theme(panel.background   = themeData[[1]],
            axis.line          = themeData[[2]],
            panel.grid.major.y = themeData[[3]],
            panel.grid.major.x = themeData[[4]])
    
    ggsave(imageFilePath, width = 8)

  }
  
}
