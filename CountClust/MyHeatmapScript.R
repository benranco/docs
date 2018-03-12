#########################################################
options(stringsAsFactors = FALSE, warn = 1)
#options(scipen = 999) # disable scientific notation for large numbers

#write("Running filterForHolly.R.", stdout())

#########################################################
# Initial Input Parameters:

#args <- commandArgs(trailingOnly = TRUE)

path <- "/home/benrancourt/Desktop/holly/RPKM-omega-theta-withRowNames"
heatmapSubfolder <- "heatmaps"
pdfSubfolder <- "pdf"
pngSubfolder <- "png"

inputFileNameBeforeNumber <- "RPKM-omega"
inputFileNameAfterNumber <- "-withRowNames.csv"

for(i in 15:35) {
  inputFileName <- paste0(inputFileNameBeforeNumber,i,inputFileNameAfterNumber)

  input <- read.csv(paste(path.expand(path),inputFileName,sep="/"),header=TRUE) ##, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.)

  rNames <- input[,1]
  input <- input[,-1]

  row.names(input) <- rNames
  inputMatrix <- as.matrix(input)

  heatmapFileName <- paste0(inputFileNameBeforeNumber,i,"-heatmap.png")
  png(filename=paste(path,heatmapSubfolder,pngSubfolder,heatmapFileName,sep="/"), width=1200, height=1000)
  heatmap(inputMatrix)
  dev.off()

  heatmapFileName <- paste0(inputFileNameBeforeNumber,i,"-heatmap.pdf")
  pdf(file=paste(path,heatmapSubfolder,pdfSubfolder,heatmapFileName,sep="/"), width=12, height=10)
  heatmap(inputMatrix)
  dev.off()

}


