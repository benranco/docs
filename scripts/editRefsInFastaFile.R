options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)


write("Running editRefsInFastaFile.R.", stdout())

##########################################################
# What this script does:
#
# - accept a fasta file as input
# - accept a .csv file as input, with four columns holding the folloing data: 
#     col1=sequenceId, col2=position, col3+reference, col3=replaceWith
# - for each row in the .csv file, search for the corresponding sequenceId+position+ref
#   in the fasta file, and replace the ref value with the replaceWith value. The replaceWith
#   value will also be surrounded by square brackets, eg. if A/G, it will be [A/G].
# - write the updated fasta file to the output file name.
#
##########################################################
# Input Parameters. Edit these as required:

# The path to the folder containing the file
path <- "/home/benrancourt/Desktop/junjun/33sampleRNA/ModifyFastaFileRefs"

# Input fasta file
#fastaFile <- "Total-seq61-from-TRI-Pipelie10.fasta"
# Output fasta file name
#outputFastaFile <- "Total-seq61-from-TRI-Pipelie10-edited.fasta"

# Input file name (it must be a .csv file)
#csvFile <- "Total-SNP94-for-put-SNP-in-seq.csv"


# Input fasta file
#fastaFile <- "ADDITIONAL-seq20nt from-TRINITY-Pipelie10.fasta"
# Output fasta file name
#outputFastaFile <- "ADDITIONAL-seq20nt from-TRINITY-Pipelie10-edited.fasta"

# Input file name (it must be a .csv file)
#csvFile <- "ADDITIONAL-seq20nt from-TRINITY-Pipelie10-variants.csv"


# Input fasta file
fastaFile <- "Total-seq81-from-TRI-Pipelie10.fasta"
# Output fasta file name
outputFastaFile <- "Total-seq81-from-TRI-Pipelie10-edited.fasta"

# Input file name (it must be a .csv file)
csvFile <- "Total-SNP114-for-put-SNP-in-seq.csv"


##############################################################
# Excecution Code. There should be no need to edit anything 
# below this line.

if(!require(seqinr)) {
  install.packages('seqinr', repos='http://cran.us.r-project.org')
}
library(seqinr)


csvData <- read.csv(paste(path.expand(path),csvFile,sep="/"),header=TRUE, check.names=FALSE)  # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

fastaData <- read.fasta(file = paste(path.expand(path),fastaFile,sep="/"), as.string = FALSE, forceDNAtolower=FALSE)
sequenceNames <- names(fastaData)

#fastaData["TRINITY_DN427_c0_g2_i11"][[1]][1:10]
#fastaData["TRINITY_DN427_c0_g2_i11"][[1]][3] <- "[A/G]"

for(csvRow in 1:nrow(csvData)) {
  seq <- csvData[csvRow, 1]
  pos <- csvData[csvRow, 2]
  ref <- csvData[csvRow, 3]
  replace <- paste0("[",csvData[csvRow, 4],"]")

  seqLength <- length(fastaData[seq][[1]])

  if ( seqLength >= pos) {
    if (fastaData[seq][[1]][pos] == ref) {
      fastaData[seq][[1]][pos] <- replace
    }
    else {
      write(paste0("WARNING: For seq=",seq,", pos=",pos,", ref=",ref,", the fasta file actually has ", fastaData[seq][[1]][pos]," instead of ",ref,". Skipping replacing it with ",replace,"."), stdout())
    }
  }
  else {
    if ( seqLength == 0) {
      write(paste0("WARNING: For seq=",seq," the fasta file doesn't have an entry."), stdout())
    } 
    else {
      write(paste0("WARNING: For seq=",seq,", pos=",pos,", ref=",ref,", the seqquence length in the fasta file is only ", seqLength,"."), stdout())
    }
  }
}

write.fasta(fastaData, sequenceNames, file.out=paste(path.expand(path),outputFastaFile,sep="/"), as.string = FALSE)

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())


