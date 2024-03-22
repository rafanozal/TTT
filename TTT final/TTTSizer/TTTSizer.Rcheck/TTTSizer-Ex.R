pkgname <- "TTTSizer"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('TTTSizer')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("getDataFromFile")
### * getDataFromFile

flush(stderr()); flush(stdout())

### Name: getDataFromFile
### Title: TTTSizer
### Aliases: getDataFromFile

### ** Examples

getDataFromFile("/home/me/myNumbers.txt")




cleanEx()
nameEx("hello")
### * hello

flush(stderr()); flush(stdout())

### Name: hello
### Title: Hello, World!
### Aliases: hello

### ** Examples

hello()



cleanEx()
nameEx("kernelFunction")
### * kernelFunction

flush(stderr()); flush(stdout())

### Name: kernelFunction
### Title: TTTSizer
### Aliases: kernelFunction

### ** Examples

kernelFunction(1.5)
          kernelFunction(0, kernel="biweight")




cleanEx()
nameEx("kernelHFunction")
### * kernelHFunction

flush(stderr()); flush(stdout())

### Name: kernelHFunction
### Title: TTTSizer
### Aliases: kernelHFunction

### ** Examples

kernelHFunction(1,2,"gaussian")




cleanEx()
nameEx("ttt")
### * ttt

flush(stderr()); flush(stdout())

### Name: ttt
### Title: For a given array of data, creates several TTT-SiZer plots -
###   Plot with the raw data using density a density plot.  - Plot with the
###   phi vector. - The Theta / Phi_h vector / Family plots for Theta0,
###   Theta1, Theta2 - All theta plot together in the same image - The
###   SiZer 0, SiZer 1, and SiZer 2 plots - All SiZer plot together in the
###   same image with the family plot Theta0 - The SiZer 0, SiZer 1, and
###   SiZer 2 plots with the z quantiles instead of strict categorical
###   pixels - The ESS (Effective Sample Space)
### Aliases: ttt

### ** Examples


This should take about 40 seconds if you run it with a CPU made from a toaster:

xgrid = 50
ygrid = 11
myData  = getRandomData(50,"gamma", 1/5, 5)

ttt(myData, xgrid, ygrid, kernel = "gaussian", myMethod = "quadratic", variance = "bootstrap",  ESSLimit = 5, saveCSV = TRUE, saveLog = FALSE)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
