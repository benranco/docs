options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)


write("Running putFastaFilesIntoCsv.R.", stdout())

##########################################################
# Input Parameters. Edit these as required:

# The path to the folder containing the file
path <- "/home/ben/pfc/junjun/ass5"

#file1 <- "test-nt.fasta"
#file2 <- "test-aa.fasta"
file1 <- "Totalmapped9,645-nt.fasta"
file2 <- "Totalmapped9,645-aa.fasta"

file1Type <- "DNA" # either "DNA" or "AA"
file2Type <- "AA" # either "DNA" or "AA"

idColHeader <- "Seq_ID"
file1ColHeader <- "DNA_Sequence"
file2ColHeader <- "Protein_Sequence"

outputFileName <- "Totalmapped9,645-dna+protein.csv"


##############################################################
# Excecution Code. There should be no need to edit anything 
# below this line.

if(!require(seqinr)) {
  install.packages('seqinr', repos='http://cran.us.r-project.org')
}
library(seqinr)

write(paste0("Processing..."), stdout())

input1 <- read.fasta(file = paste(path.expand(path),file1,sep="/"), as.string = TRUE, seqtype=file1Type, forceDNAtolower=FALSE, set.attributes=FALSE, seqonly=FALSE, strip.desc=TRUE, whole.header=FALSE)

input1Df <- data.frame(names(input1), unlist(unname(input1)))
names(input1Df) <- c(idColHeader, file1ColHeader)


input2 <- read.fasta(file = paste(path.expand(path),file2,sep="/"), as.string = TRUE, seqtype=file2Type, forceDNAtolower=FALSE, set.attributes=FALSE, seqonly=FALSE, strip.desc=TRUE, whole.header=FALSE)

input2Df <- data.frame(names(input2), unlist(unname(input2)))
names(input2Df) <- c(idColHeader, file2ColHeader)


output <- merge(input1Df, input2Df, by=idColHeader, all=TRUE)

write.csv(output, file=paste(path.expand(path),outputFileName,sep="/"), row.names=FALSE )


write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())


