options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]

write("Running selectSeqsFromFasta.R.", stdout())

# The path to the folder containing the files. All file paths should be relative to this path, 
# even if they are not in the folder itself.
#path <- "./test"
path <- args[1]

# List of sequence ids, one per line:
#seqIdFile <- "WWP-longestSeqs-1000.txt"
seqIdFile <- args[2]

# Input fasta file
#fastaFile <- "../../trin-assemblies-renamedSeqIds/WWP18-nt675,367-Trinity.fasta"
fastaFile <- args[3]

# Output fasta file name
#outputFastaFile <- "WWP-fromR-1000.fasta"
outputFastaFile <- args[4]



##############################################################
# Excecution Code. There should be no need to edit anything 
# below this line.

if(!require(seqinr)) {
  install.packages('seqinr', repos='http://cran.us.r-project.org')
}
library(seqinr)


write(paste0("Input seqIdFile: ",seqIdFile), stdout())
write(paste0("Input fasta file: ",fastaFile), stdout())

write(paste0("Reading the seqIdFile."), stdout())
seqIds <- read.table(file=paste(path.expand(path),seqIdFile,sep="/"), header=FALSE)

write(paste0("Reading the fasta file..."), stdout())
fastaData <- read.fasta(file = paste(path.expand(path),fastaFile,sep="/"), as.string = TRUE, strip.desc=TRUE,  forceDNAtolower=FALSE)

write(paste0("Selecting only the sequences from the input sequence ids."), stdout())
selectedFastaData <- fastaData[seqIds[,1]]

write(paste0("Writing output file: ",outputFastaFile), stdout())
write.fasta(selectedFastaData, getAnnot(selectedFastaData), file.out=paste(path.expand(path),outputFastaFile,sep="/"), as.string = FALSE)
# I noticed that when write.fasta's as.string = FALSE while read.fasta's as.string = TRUE, it writes each sequence out 
# on just one line.

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())


