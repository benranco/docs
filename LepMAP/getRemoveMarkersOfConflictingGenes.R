


options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)


write("Running getRemoveMarkersOfConflictingGenes.R.", stdout())




# ########################################################
# Input Parameters:

path <- "/home/benrancourt/Downloads/JunJun-ResponseToPaperFeedback/JulyProcessing/removeMarkersOfConflictingGenes"

inputFile <- "r38-less10%-SNP-IDs.csv"
inputFileConflictingGenes <- "r38-less10%-SNP-IDs-ConflictingGenes.csv"

outputFile <- "r38-less10%-SNP-IDs-conflictingGenesRemoved.csv"


# separatorForGeneNameVsMarkerId: This parameter relates to the marker ids in the original
# .csv data file that you supplied. The marker id column in your input .csv file contains
# marker ids, which are composed of 
# a gene id joined to a string that uniquely identifies each marker of a gene. There should
# be a separator character between the gene portion of the id and the marker portion of the
# id. For example, the marker id "*M100356-1324K" is composed of two parts: the gene part
# (*M100356), and the marker part (1324K), separated by a "-". Common separators are "-" and
# "_". All id's in your input .csv file must use the same separator character to
# differentiate the gene portion from the marker portion. This script uses the separator 
# character to extract the gene id from the whole marker id string. It then uses the gene
# id when tabulating it's final output data. Please indicate what your input .csv file uses
# as the separator character.
separatorForGeneNameVsMarkerId <- "_"

# missingDataIndicator: This parameter relates to the data in the original
# .csv data file that you supplied. Some of the data in your original input .csv will
# likely be missing, in which case it might be represented as NA or "-" or maybe even be
# empty, or something like that to indicate that it is missing.
# Indicate what your .csv file uses to represent missing data.
missingDataIndicator <- "-"


# End of input parameters.
##############################################################

if(!require(stringr))
{
  install.packages('stringr', repos='http://cran.us.r-project.org')
}

library(stringr)


# data cells in the input .csv files (both fathers and seeds) which contain these 
# strings will be interpreted as NA/missing data:
myNaString <- "-"

input <- read.csv(paste(path.expand(path),inputFile,sep="/"),header=TRUE, na.strings=myNaString, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

inputConflictingGenes <- read.csv(paste(path.expand(path),inputFileConflictingGenes,sep="/"),header=TRUE, na.strings=myNaString, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

conflictGenes <- inputConflictingGenes[,1]
conflictGenesConcatenated <- paste(conflictGenes, collapse="|")

# create a new column in inputWithAllMarkers that just has the gene names, derived from the marker_id column.
# library(stringr) required for str_split_fixed.
# This assumes that column 1 of input is the marker id
gene <- str_split_fixed(input[,1],separatorForGeneNameVsMarkerId,2)[ ,1]
input <- cbind(gene, input)

rowsToRemove <- logical(nrow(input)) # initializes all to FALSE

for (i in 1:nrow(input)) {
  curGene <- input[i,"gene"]

  for (k in 1:length(conflictGenes)) {
    if (curGene == conflictGenes[k]) {
      rowsToRemove[i] <- TRUE
      break
    }
  }
}

rowsToKeep <- !rowsToRemove
output <- input[rowsToKeep, -1 ] # also remove the gene column that we created


write.csv(output, paste(path.expand(path), outputFile,sep="/"), row.names=FALSE, na=myNaString)

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())



