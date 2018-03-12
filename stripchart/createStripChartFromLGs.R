

options(stringsAsFactors = FALSE, warn = 1)

write("Running createStripChartFromLGs.R.", stdout())

if(!require(beanplot))
{
  install.packages('beanplot', repos='http://cran.us.r-project.org')
}

library(beanplot)

#  mybeanplot.R replaces the beanplot function with a customized version, modified by Ben.
source("mybeanplot.R")


# ########################################################
# Input Parameters:

location <- "/home/benrancourt/Downloads/stripcharts"

#outputFileNameSansPostfix <- "plot3"

#outputFileNameSansPostfix <- "Plot-R38-m5,006-density"
#outputFileNameSansPostfix <- "Plot-R38-m5,186-density"
#outputFileNameSansPostfix <- "Plot-R45-m7,727-density"
#outputFileNameSansPostfix <- "Plot-LP-Merge-m9,511-density"

#outputFileNameSansPostfix <- "Plot-LP-Merge-m9,511-density-withRGA639-overlay"
#outputFileNameSansPostfix <- "Plot-LP-Merge-m9,511-density-withRGA639-underneath"
#outputFileNameSansPostfix <- "Plot-LP-Merge-m9,511-density-withRLK290-overlay"
#outputFileNameSansPostfix <- "Plot-LP-Merge-m9,511-density-withRLK290-underneath"

#outputFileNameSansPostfix <- "Plot-R38-m5,186-density-withRGA380-overlay"
outputFileNameSansPostfix <- "Plot-R38-m5,186-density-withRGA380-underneath"

# We assume the input csv file colums are organized: 
# Column 1: LG, Column 2: position
#inputFileName <- "R38-m5,006.csv"
inputFileName <- "R38-m5,186.csv"
#inputFileName <- "R45-m7,727.csv"
#inputFileName <- "LP-Merge-m9,511.csv"

hasOverlayData <- TRUE
#overlayInputFileName <- "RGA639.csv" # to overlay the 9,511 data
#overlayInputFileName <- "RLK290.csv" # to overlay the 9,511 data
overlayInputFileName <- "RGA380-from-R38-g5,186.csv" # to overlay the R38 5,186 data

hasSecondOverlayData <- FALSE
#secondOverlayInputFileName <- "RGA639.csv"
#secondOverlayInputFileName <- "RLK290.csv"

plotOverlayUnderneath <- TRUE # True = underneath, FALSE = overlay

lgLabelPrefix <- "LG"

if (!hasOverlayData)
{
  firstColor <- "black"
} else 
{
  firstColor <- "lightblue"
}
secondColor <- "red"
thirdColor <- "green"

includeLegend <- FALSE

# legend names are just taken from the input file names, unless you specify custom names here
if (hasOverlayData)
{
  legendNameMainData <- strsplit(inputFileName,".csv",fixed=TRUE)[[1]]
  legendNameOverlay <- strsplit(overlayInputFileName,".csv",fixed=TRUE)[[1]]
  if (hasSecondOverlayData)
  {
    legendNameSecondOverlay <- strsplit(secondOverlayInputFileName,".csv",fixed=TRUE)[[1]]
  }
}

# End of Input Parameters.
# ########################################################


#
# Function mapData
#
mapData <- function(path,dataFileName,isOverlay=FALSE,overlayUnderneath=FALSE,tickColor="black")
{

  input <- read.csv(paste(path.expand(path),dataFileName, sep="/"),header=TRUE)

  # We assume the input csv file colums are organized: Column 1: LG, Column 2: position

  LGsInData <- sort(unique(input[ , 1]), decreasing=TRUE)

  write("LG's in data: ", stdout())
  write(LGsInData, stdout())
  
  positions <- list()
  lgNames <- c()

  minPositions <- c()
  maxPositions <- c()
  bandwidths <- c()

  for (i in 1:length(LGsInData))
  {
    pos <- input[ input[,1]==LGsInData[i], 2] # grab all the input data for the current LG
    positions <- c(positions, list(pos)) # add pos as a new vector element to the positions list by first creating a list of one with pos in it, and then concatenating it to positions
    lgNames[i] <- paste0(lgLabelPrefix,LGsInData[i])
    minPositions[i] <- min(pos)
    maxPositions[i] <- max(pos)
  }

  names(positions) <- lgNames

  # Select Usage info for beanplot:
  #     beanplot(..., bw = "SJ-dpi", kernel = "gaussian", cut = 3, cutmin = -Inf, 
  #              cutmax = Inf, grownage = 10, what = c(TRUE, TRUE, TRUE, TRUE), 
  #              add = FALSE, col, axes = TRUE, log = "auto", handlelog = NA, 
  #              ll = 0.16, wd = NA, maxwidth = 0.8, maxstripline = 0.96, 
  #              method = "stack", names, overallline = "mean", beanlines = overallline, 
  #              horizontal = FALSE, side = "no", jitter = NULL, beanlinewd = 2, 
  #              frame.plot = axes, border = NULL, innerborder = NA, at = NULL, 
  #              boxwex = 1, ylim = NULL, xlim = NULL, show.names = NA)
  #
  #  cutmin: the low-ends of the beans are cut below ‘mincut*bw’. Defaults
  #          to ‘cut’.
  #
  #    what: a vector of four booleans describing what to plot. In the
  #          following order, these booleans stand for the total average
  #          line, the beans, the bean average, and the beanlines. For
  #          example, ‘what=c(0,0,0,1)’ produces a ‘stripchart’
  #
  #     col: the colors to be used. A vector of up to four colors can be
  #          used. In the following order, these colors stand for the area
  #          of the beans (without the border, use ‘border’ for that
  #          color), the lines inside the bean, the lines outside the
  #          bean, and the average line per bean. Transparent colors are
  #          supported. ‘col’ can be a list of color vectors, the vectors
  #          are then used per bean, and reused if necessary.
  #
  #     ll: the length of the beanline per point found.
  #
  #     maxwidth: the maximum width of a bean.
  #
  #     maxstripline: the maximum length of a beanline.
  #
  #    side: the side on which the beans are plot. Default is ‘"no"’, for
  #          symmetric beans. ‘"first"’, ‘"second"’ and ‘"both"’ are also
  #          supported.
  #
  # frame.plot: if true, plots a frame.
  #
  #  border: the color for the border around a bean. ‘NULL’ for
  #          ‘par("fg")’, ‘NA’ for no border. If border is a vector, it
  #          specifies the color per bean, and colors are reused if
  #          necessary.
  #
  # innerborder: a color (vector) for the border inside the bean(s).
  #          Especially useful if ‘side="both"’. Use ‘NA’ for no inner
  #          border. Colors are reused if necessary. The inner border is
  #          drawn as the last step in the drawing process.
  #
  #      at: the positions at which a bean should be drawn.

  # The las, bty and all cex parameters refer to font orientation and sizes, and are documented in help(par)

  #myFramePlot <- FALSE
  #myBty <- "l"

  myWhat <- c(F,T,F,T)
  myInnerborder <- NA #"black"
  myAdd <- FALSE
  mySide <- "second"
  myMaxwidth <- 1.8
  myLl <- 1.8
  if (overlayUnderneath)
  {
    myLl <- 0.8 # this sets the line height, so needs to be adjusted even if this layer isn't the overlay
    myMaxwidth <- 1 # this sets the height of the density plot, so needs to be adjusted even if this layer isn't the overlay
  }

  # adjust some of the input parameters for beanplot if this is an overlay layer to be plotted 
  # on top of an already existing plot
  if (isOverlay)
  {
    myWhat <- c(F,F,F,T)
    myInnerborder <- NA
    #myFramePlot <- FALSE
    myAdd <- TRUE
    mySide <- "second"

    if (overlayUnderneath)
    {
      mySide <- "first"
    }
  }
  
  beanplot(positions, kernel="gaussian", horizontal=TRUE, add=myAdd, col=c("lemonchiffon",tickColor,tickColor), what=myWhat, method="overplot", grownage=40,ll=myLl, side=mySide, innerborder=myInnerborder, maxwidth=myMaxwidth, cutmin=minPositions, cutmax=maxPositions, bw="SJ-dpi", las=1, cex.axis=1.1)


  if (!isOverlay)
  {
    # draw the long base lines explicitly
    for (i in 1:length(LGsInData))
    {
      #lines( c(i,i), c(minPositions[i],maxPositions[i]), lwd=1) # vertical
      lines( c(minPositions[i],maxPositions[i]), c(i,i), lwd=1) # horizontal
    }

    ## draw custom lines/ticks for each data point
    #for (i in 1:length(positions))
    #{
    #  pos <- positions[[i]]
    #  for (k in 1:length(pos)) 
    #  {
    #    lines( c(pos[k],pos[k]), c(i-0.3,i+0.3), lwd=1)
    #  }
    #}
  }

} # end function mapData


#
# Function generatePlot
#
generatePlot <- function()
{
  # set graph margins: bottom,left,top,right
  par(mar=c(2,3,0,0)+0.4)

  # draw the main layer
  mapData(location, inputFileName, isOverlay=FALSE, overlayUnderneath=plotOverlayUnderneath, tickColor=firstColor)

  # draw the second layer, if it exists
  if (hasOverlayData)
  {
    mapData(location, overlayInputFileName, isOverlay=TRUE, overlayUnderneath=plotOverlayUnderneath, tickColor=secondColor)
  }

  # draw the third layer, if it exists
  if (hasSecondOverlayData)
  {
    mapData(location, secondOverlayInputFileName, isOverlay=TRUE, overlayUnderneath=plotOverlayUnderneath, tickColor=thirdColor)
  }

  if (includeLegend)
  {
    if (hasOverlayData && !hasSecondOverlayData)
    {
      legend("topright", inset=0.02, legend=c(legendNameMainData, legendNameOverlay), col=c(firstColor, secondColor), lty=c(1,1), box.lty=0, cex=1.1)
    } else if (hasOverlayData && hasSecondOverlayData) 
    {
      legend("topright", inset=0.02, legend=c(legendNameMainData, legendNameOverlay, legendNameSecondOverlay), col=c(firstColor, secondColor, thirdColor), lty=c(1,1,1), box.lty=0, cex=1.1)
    } else
    {
      legend("topright", inset=0.02, legend=c(legendNameMainData), col=c(firstColor), lty=c(1), box.lty=0, cex=1.1)
    }
  }

} # end function generatePlot



myHeightPng <- 800
myHeightPdf <- 8
myHeightTiff <- 8
if (hasOverlayData && plotOverlayUnderneath)
{
  myHeightPng <- 900
  myHeightPdf <- 9
  myHeightTiff <- 9
}

# create png image
write("-----------------", stdout())
write("Creating png.", stdout())
png(filename=paste(path.expand(location),paste0(outputFileNameSansPostfix,"-web.png"),sep="/"), width=1200, height=myHeightPng)
generatePlot()
dev.off()

# create pdf image (useful if you need to scale it extra large)
write("-----------------", stdout())
write("Creating pdf (scalable).", stdout())
pdf(file=paste(path.expand(location),paste0(outputFileNameSansPostfix,"-scalable.pdf"), sep="/"), width=12, height=myHeightPdf)
generatePlot()
dev.off()

tiffFileName <- paste(path.expand(location),paste0(outputFileNameSansPostfix,"-print-600dpi-12inch.tiff"), sep="/")

# create tiff image
write("-----------------", stdout())
write("Creating large tiff (12 inches wide at 600 dpi).", stdout())
tiff(file=tiffFileName, width=12, height=myHeightTiff, units = 'in', res=600, compression = "lzw", pointsize=10)
generatePlot()
dev.off()

smallerTiffFileName <- paste(path.expand(location),paste0(outputFileNameSansPostfix,"-print-600dpi-6inch.tiff"), sep="/")
# Now create a scaled down tiff at 6x6 inches. I had to export it first as a 12x12 inch tiff (above) 
#(7200 pixels at 600 dpi) because if I try to export it as a 6x6 inch tiff the lines are too thick 
# and it doesn't look good. To created a version scaled to 6x6 inches at 600 dpi, use the following 
# command-line command:
#    convert inputFile.tiff -resize 3600x3600 outputFile.tiff
# (ImageMagick has to be installed for this to work, which you can install simply by installing
# Inkscape via (on CentOS):
#    sudo yum install inkscape
write("-----------------", stdout())
write("Creating scaled down tiff (6 inches wide at 600 dpi).", stdout())
system(paste0("convert ",tiffFileName," -resize 3600x2700 ",smallerTiffFileName))

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())


