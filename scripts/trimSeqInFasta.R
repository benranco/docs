options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)


write("Running trimSeqInFasta.R.", stdout())

##########################################################
# What this script does:
#
# - accept a fasta file as input
# - accept a .csv file as input, with four columns holding the folloing data: 
#     col1=seqName, col2=Hit length, col3=Query start, col4=Query end
# - for each row in the .csv file, find the corresponding sequenceId
#   in the fasta file, and save a subsequence of it beginning at Query start 
#   and ending at Query end.
# - write an updated fasta file containing only the trimmed subsequences.
# NOTE: duplicate sequence ids will be given unique names in the output fasta.
# NOTE: sequence ids that are not found in the input fasta will not be included.
# NOTE: other error checking and warnings may occur. See the console output.
#
##########################################################
# Input Parameters. Edit these as required:

# The path to the folder containing the file
path <- "/home/ben/pfc/junjun/ass1"

# Input fasta file
#fastaFile <- "LG12-LPss-527cds.fasta"
# Output fasta file name
#outputFastaFile <- "LG12-LPss-527cds-trimmed.fasta"

# Input file name (it must be a .csv file)
#csvFile <- "LG12-seq480-for trim.csv"


# Input fasta file
fastaFile <- "Cr4cM5-Seg491nt-for extr-exonic496.fasta"
# Output fasta file name
outputFastaFile <- "Cr4cM5-Seg491nt-for extr-exonic496-trimmed.fasta"

# Input file name (it must be a .csv file)
csvFile <- "Cr4cM5-Totalseq491-region496.csv"

##############################################################
# Excecution Code. There should be no need to edit anything 
# below this line.

if(!require(seqinr)) {
  install.packages('seqinr', repos='http://cran.us.r-project.org')
}
library(seqinr)


write(paste0("Input .csv file: ",csvFile), stdout())
write(paste0("Input fasta file: ",fastaFile), stdout())

fastaData <- read.fasta(file = paste(path.expand(path),fastaFile,sep="/"), as.string = FALSE, forceDNAtolower=FALSE)


naStrings <- c("NA","na",""," ","-")
csvData <- read.csv(paste(path.expand(path),csvFile,sep="/"), na.strings=naStrings,header=TRUE, check.names=FALSE)  # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

# If there are any sequence names in the .csv file that aren't found in the fasta 
# file, remove these rows before proceeding.
rowsToRemove <- logical(nrow(csvData)) # initializes all to false
for(csvRow in 1:nrow(csvData)) {
  name <- csvData[csvRow,1]
  rowsToRemove[csvRow] <- is.null(fastaData[name][[1]]) # sets to true if data for name is null
}
if (length(which(rowsToRemove)) >= 1) {
  write("WARNING: Removing the following sequence ids from the input id list because they were not found in the fasta file: ", stdout())
  write(paste0(csvData[which(rowsToRemove),1],collapse="\n"), stdout())
  csvData <- csvData[!rowsToRemove,] # remove all rowsToRemove
}

# If there are multiple occurences of some sequence names, make sure each occurence is
# given a new unique name to use for the output fasta data.
newSequenceNames <- csvData[,1]
#for (i in 1:length(newSequenceNames)) {
for (i in 1:nrow(csvData)) {
  name <- newSequenceNames[i]
  occurences <- grep(name, newSequenceNames) # returns a list of indices of occurences of that name
  if (length(occurences) > 1) {
    for (k in 1:length(occurences)) {
      newSequenceNames[occurences[k]] <- paste0(name,"-f",k)
    }
  }
}

# character vector to store the trimmed fasta data in:
trimmedFastaData <- vector(mode="list", length=nrow(csvData))
names(trimmedFastaData) <- newSequenceNames



for(csvRow in 1:nrow(csvData)) {
  seq <- csvData[csvRow, 1]
  hitLength <- as.numeric(csvData[csvRow, 2])
  qStart <- as.numeric(csvData[csvRow, 3])
  qEnd <- as.numeric(csvData[csvRow, 4])

  seqLength <- length(fastaData[seq][[1]])
  sequence <- fastaData[seq][[1]]
  
  if (is.na(qStart) ) { 
    qStart <- 1 
  }
  if (is.na(qEnd) ) { 
    qEnd <- seqLength 
  }  

  if (qStart < 1 || qStart > seqLength) {
    qStart <- 1
    write(paste0("WARNING: For seq=",seq," the Query start is beyond the bounds of the sequence length. Setting Query start to 1."), stdout())
  }
  
  if (qEnd < 1 || qEnd > seqLength) {
    qEnd <- seqLength
    write(paste0("WARNING: For seq=",seq," the Query end is beyond the bounds of the sequence length. Setting Query end to the sequence length: ",seqLength,"."), stdout())
  }
  
  if (qEnd < qStart) {
    qEnd <- seqLength
    write(paste0("WARNING: For seq=",seq," the Query end is less than the Query start. Setting Query end to the sequence length: ",seqLength,"."), stdout())
  }
  
#  write(paste0("seq: ",seq), stdout())
#  write(paste0("hitLength: ",hitLength), stdout())
#  write(paste0("qStart: ",qStart), stdout())
#  write(paste0("qEnd: ",qEnd), stdout())
#  write(paste0("seqLength: ",seqLength), stdout())
#  write(paste0("value at qStart: ",sequence[qStart]), stdout())
#  write(paste0("value at qEnd: ",sequence[qEnd]), stdout())
#  write(paste0("sequence: \n\n",paste0(sequence,collapse="")), stdout())
  
  sequence <- sequence[qStart:qEnd]
  trimmedFastaData[csvRow][[1]] <- sequence
  
  if (!is.na(hitLength) && hitLength != length(sequence)) {
    write(paste0("WARNING: For seq=",seq," the trimmed sequence length is ",length(sequence)," but the recorded Hit length is ",hitLength,"."), stdout())
  }
  
#  write(paste0("trimmed sequence: \n\n",paste0(sequence,collapse="")), stdout()) 
}

write.fasta(trimmedFastaData, names(trimmedFastaData), file.out=paste(path.expand(path),outputFastaFile,sep="/"), as.string = FALSE)

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())


