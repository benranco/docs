options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)

write("Running generateAllPairsInEachSet.R.", stdout())


# ########################################################
# Input Parameters:

path <- "/home/ben/pfc/junjun/ass4/step2"

# The input file should be a csv file containing two columns (header line is expected but column names are unimportant): 
#   col1==identifies which set/cluster an item is to be included in when pairing all items in each set/cluster.
#   col2==the individual items
inputFile <- "Envelop158-clusters-seqIds.csv"

# The output is a simple text file listing one pair per row, with the pairs delimited by the specified delimiter
outputFile <- "Envelop158-clusters-allPairs.tab"
outputDelimiter <- "\t"


# ########################################################
# Execution Code:

data <- read.csv(paste(path.expand(path),inputFile, sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

setNames <- unique(data[,1])

totalPairs <- 0
firstWrite <- TRUE

for (name in setNames) {
  
  setItems <- data[(data[1]==name),2]
  write(paste0("Creating ", choose(length(setItems), 2), " pairs in set: ", name), stdout())
  totalPairs <- totalPairs + choose(length(setItems), 2)
  
  pairsMatrix <- combn(setItems,2)
  
  pairsVector <- character(ncol(pairsMatrix))
  for (i in 1:length(pairsVector)) {
    pairsVector[i] <- paste(pairsMatrix[,i], collapse=outputDelimiter)    
  }
  
  if (firstWrite) {
    firstWrite <- FALSE
    write(pairsVector, file=paste(path.expand(path),outputFile,sep="/"), append=FALSE )
  } else {
    write(pairsVector, file=paste(path.expand(path),outputFile,sep="/"), append=TRUE )
  }
  
} # end outer for-loop

write(paste0("The total number of pairs written to the ouput file should be: ",totalPairs), stdout())

write(paste0("FINISHED."), stdout())


