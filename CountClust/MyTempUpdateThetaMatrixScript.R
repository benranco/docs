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

inputFileNameBase <- "RPKM-theta"
inputFileNameEnd <- "-withRowNames.csv"

outputFileNameBase <- "RPKM-theta"
outputFileNameEnd <- "-withRowNames.csv"

originalDataFileName <- "RPKM_matrix_8101_8112_8206_8214.csv"

origData <- read.csv(paste(path.expand(path),originalDataFileName,sep="/"),header=TRUE) 
origDataRowNames <- origData[,1]


for(i in 4:14) {
  inputFileName <- paste0(inputFileNameBase,i,inputFileNameEnd)

  input <- read.csv(paste(path.expand(path),thetaSubFolder,inputFileName,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.) 

  rNames <- input[,1]
  input <- input[,-1]

  #row.names(input) <- rNames

  input <- cbind(rNames, origDataRowNames, input)

  outputFileName <- paste0(outputFileNameBase,i,outputFileNameEnd)

  write.csv(input, paste(path.expand(path), topicsSubFolder, outputFileName,sep="/"), row.names=FALSE)
}












