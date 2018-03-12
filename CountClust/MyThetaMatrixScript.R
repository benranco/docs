#########################################################
options(stringsAsFactors = FALSE, warn = 1)
#options(scipen = 999) # disable scientific notation for large numbers

#write("Running filterForHolly.R.", stdout())

#########################################################
# Initial Input Parameters:

#args <- commandArgs(trailingOnly = TRUE)

path <- "/home/benrancourt/Desktop/holly"
thetaSubFolder <- "RPKM-omega-theta-withRowNames/theta"
topicsSubFolder <- "RPKM-omega-theta-withRowNames/theta/topics"

inputFileName <- "RPKM-theta6-withRowNames.csv"
interestingColNames <- c("X1","X3","X4")

outputFileNameBase <- "RPKM-theta6-topic"
outputFileNameEnd <- ".csv"

#originalDataFileName <- "RPKM_matrix_8101_8112_8206_8214.csv"

#origData <- read.csv(paste(path.expand(path),originalDataFileName,sep="/"),header=TRUE) 
#origDataRowNames <- origData[,1]

input <- read.csv(paste(path.expand(path),thetaSubFolder,inputFileName,sep="/"),header=TRUE) 
# this assumes the original data is ordered the same as the theta matrix:
#row.names(input) <- origDataRowNames 


for(i in 1:length(interestingColNames)) {
  
  output <- input[ , c(names(input)[[1]], names(input)[[2]], interestingColNames[i]) ]
  output <- output[ order(output[,interestingColNames[i]]), ]

  outputFileName <- paste0(outputFileNameBase,interestingColNames[i],outputFileNameEnd)

  write.csv(output, paste(path.expand(path), topicsSubFolder, outputFileName,sep="/"), row.names=FALSE)

}









