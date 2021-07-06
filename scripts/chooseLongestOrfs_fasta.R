options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)

write("Running chooseLongestOrfs_fasta.R.", stdout())

# Call it like this: 
# Rscript chooseLongestOrfs_fasta.R $outDir $allGeneIds $fastaData $pepData $trinityFile $finalOutputFasta

# The path to the folder containing the files. All file paths should be relative to this path, 
# even if they are not in the folder itself.
#path <- "./test3"
path <- args[1]

# List of gene ids, one per line:
#geneIdFile <- "test3-wwp-allGeneIds.txt"
geneIdFile <- args[2]

# Selected data from the fasta file (from Trinity) sequence headers: geneId,seqId,seqLength
#fastaDataCsvFile <- "test3-wwp-fastaData.csv"
fastaDataCsvFile <- args[3]

# Selected data from the pep file (from TransDecoder) sequence headers: geneId,seqId,seq_p,orfLength
#pepDataCsvFile <- "test3-wwp-pepData.csv"
pepDataCsvFile <- args[4]

# Input fasta file (from Trinity)
#fastaFile <- "wwp-n20.fasta"
fastaFile <- args[5]

# Output fasta file name
#outputFastaFile <- "selected.fasta"
outputFastaFile <- args[6]

#outFilePrefix <- "test3-wwp-"
outFilePrefix <- args[7]

##############################################################
# Excecution Code. There should be no need to edit anything 
# below this line.

if(!require(seqinr)) {
  install.packages('seqinr', repos='http://cran.us.r-project.org')
}
library(seqinr)

out_longestOrfSeqsAll <- paste0(outFilePrefix,"longestOrfSeqs-all.txt")
out_longestOrfSeqsFromPep <- paste0(outFilePrefix,"longestOrfSeqs-fromPep.txt")
out_longestDNASeqsFromFasta <- paste0(outFilePrefix,"longestDNASeqs-fromFasta.txt")

write(paste0("Input geneIdFile: ",geneIdFile), stdout())
write(paste0("Input fastaDataCsvFile: ",fastaDataCsvFile), stdout())
write(paste0("Input pepDataCsvFile: ",pepDataCsvFile), stdout())
write(paste0("Input fasta file: ",fastaFile), stdout())

write(paste0("Reading the geneIdFile."), stdout())
geneIds <- read.table(file=paste(path.expand(path),geneIdFile,sep="/"), header=FALSE)

write(paste0("Reading the fasta info csv file."), stdout())
fastaInfo <- read.csv(paste(path.expand(path),fastaDataCsvFile,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

write(paste0("Reading the pep info csv file."), stdout())
pepInfo <- read.csv(paste(path.expand(path),pepDataCsvFile,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.


###################################################

# initialize chosenSeqs as dataframe containing a character vector (seqId) and logical vector (hasOrf) of same length as geneIds
# for each gene
#   initialize variables longestLen and longestOrfSeq to 0/NA
#   get all rowNums in pepInfo that have the geneId
#   for each of these rowNums
#     if its orfLength is the longest orf so far, save the length and the seqId
#     if its orfLegnth ties the longest orf so far:
#       get the length of the DNA sequence of the existing longestOrfSeq from fastaInfo
#       get the length of the DNA sequence of the contending longestOrfSeq from fastaInfo
#       if the contender seqLength is the longest of the two, save its seqId (orf length stays the same)
#
#   if the variable longestOrfSeq is still NA:
#     set chosenSeqs[geneNum,"hasOrf"] <- FALSE
#     get all rowNums in fastaInfo that have the geneId
#     for each of these rowNums
#       if its seqLength is the longest so far, save as longestOrfSeq
#   else 
#     set chosenSeqs[geneNum,"hasOrf"] <- TRUE
#
#   save longestOrfSeq in the matching spot in chosenSeqs
# end for-loop
#
# write out three text files: out_longestOrfSeqsAll, out_longestOrfSeqsFromPep, out_longestDNASeqsFromFasta

write(paste0(Sys.time(), " - Beginning comparisons..."), stdout())

chosenSeqs <- data.frame(seqId=character(nrow(geneIds)), hasOrf=logical(nrow(geneIds)))

for (geneNum in 1:nrow(geneIds)) {

  longestLen <- 0
  longestOrfSeq <- NA
  
  gene <- geneIds[geneNum,]
  
  pep <- pepInfo[pepInfo$geneId==gene,]
  
  if (nrow(pep) > 0) {
    for (rowNum in 1:nrow(pep)) {
      if (pep[rowNum,"orfLength"] > longestLen) {
        longestLen <- pep[rowNum,"orfLength"]
        longestOrfSeq <- pep[rowNum,"seqId"]
      }
      else if (pep[rowNum,"orfLength"] == longestLen) {
        oldDnaLength <- fastaInfo[fastaInfo$seqId == longestOrfSeq, "seqLength"]
        newDnaLength <- fastaInfo[fastaInfo$seqId == pep[rowNum,"seqId"], "seqLength"]
        if (newDnaLength > oldDnaLength) {
          longestOrfSeq <- pep[rowNum,"seqId"]
        }
      }    
    } # end inner for-loop
  }
  
  if (is.na(longestOrfSeq)) {
    # chosenSeqs[geneNum,"hasOrf"] <- FALSE # this line isn't needed because it was all initialzed to FALSE
    fa <- fastaInfo[fastaInfo$geneId==gene,]    
    longestOrfSeq <- fa[which.max(fa[,"seqLength"]), "seqId"]
  } 
  else { 
    chosenSeqs[geneNum,"hasOrf"] <- TRUE
  }
  
  chosenSeqs[geneNum,"seqId"] <- longestOrfSeq
  
} # end outer for-loop

write(paste0(Sys.time(), " - done."), stdout())

write(paste0("Writing selected sequence ids."), stdout())

write.table(chosenSeqs[,"seqId"], file=paste(path.expand(path),out_longestOrfSeqsAll,sep="/"), row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(chosenSeqs[chosenSeqs$hasOrf,"seqId"], file=paste(path.expand(path),out_longestOrfSeqsFromPep,sep="/"), row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(chosenSeqs[!chosenSeqs$hasOrf,"seqId"], file=paste(path.expand(path),out_longestDNASeqsFromFasta,sep="/"), row.names=FALSE, col.names=FALSE, quote=FALSE)


###################################################

write(paste0(Sys.time(), " - Reading the fasta file..."), stdout())
fastaData <- read.fasta(file = paste(path.expand(path),fastaFile,sep="/"), as.string = TRUE, strip.desc=TRUE,  forceDNAtolower=FALSE)

write(paste0(Sys.time(), " - Grabbing only the selected sequences from the fasta file."), stdout())
selectedFastaData <- fastaData[chosenSeqs[,"seqId"]]

write(paste0(Sys.time(), " - done."), stdout())

write(paste0("Writing output fasta file file: ",outputFastaFile), stdout())
write.fasta(selectedFastaData, getAnnot(selectedFastaData), file.out=paste(path.expand(path),outputFastaFile,sep="/"), as.string = FALSE)
# I noticed that when write.fasta's as.string = FALSE while read.fasta's as.string = TRUE, it writes each sequence out 
# on just one line.

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())


