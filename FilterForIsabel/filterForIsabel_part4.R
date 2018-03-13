options(stringsAsFactors = FALSE, warn = 1)
options(scipen = 999) # disable scientific notation for large numbers

args <- commandArgs(trailingOnly = TRUE)


write("Running filterForIsabel_part4.R.", stdout())

# #######################################################################

#basePath <- "/run/media/benrancourt/27a1f1ff-6590-4f38-aa2c-13f6f1e0e771/ben/SNPpipeline-isabel"
basePath <- args[1]

#pathToReports <- paste(path.expand(basePath),"reports/reports-backup/done-exceptMissingData/filtered",sep="/")
pathToReports <- paste(path.expand(basePath),args[2],sep="/")

#inputReportName <- "filtered_by_depth_recoded.csv"
#inputReportName <- "test.csv"
inputReportName <- args[3]

#outputFastaName <- "filtered_by_depth_recoded.fasta"
#outputFastaName <- "test.fasta"
outputFastaName <- args[4]


# #######################################################################

inputReport <- read.csv(paste(path.expand(pathToReports),inputReportName,sep="/"),header=TRUE)
# check if first col is not CHROM, and remove it if so (it's just a left-over old index column if it exists)
if (colnames(inputReport)[1] != "CHROM")
{
  inputReport <- inputReport[,-1]
}

# we don't need the first three columns to make the fasta file:
inputReport <- inputReport[,c(-1,-2,-3)]

# transpose:
inputReport <- t(inputReport)

sampleNames <- row.names(inputReport)
# we just want to use a portion of the sample name
sampleNames <- strsplit(sampleNames,"\\.")
func <- function(vec) {
  vec[length(vec)-1]
}
sampleNames <- sapply(sampleNames, func) # take just the second last element of each vector returned by strsplit
sampleNames <- paste0(">",sampleNames)

sequences <- character(nrow(inputReport))

for (i in 1:nrow(inputReport)) {
  seq <- paste0(inputReport[i,],collapse="")
  # Add a \n every 60 characters. Explanation: 
  # (.{60}) means match any character (.) up to 60 times, \\1 means back-reference the 
  # last match and use it in the substitution string, \n adds a \n to the substitution string.
  seq <- gsub("(.{60})","\\1\n",seq)
  sequences[i] <- seq
}

#interleave the sampleNames and sequences
outputData <- as.vector(rbind(sampleNames,sequences))

write(outputData, file=paste(path.expand(pathToReports), outputFastaName,sep="/"), sep="", append=FALSE)

