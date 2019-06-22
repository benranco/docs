


options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)


write("Running getAllMarkersThatMatchFinalUniqueGenes.R.", stdout())




# ########################################################
# Input Parameters:

path <- "/home/benrancourt/Downloads/JunJun-ResponseToPaperFeedback/JulyProcessing/R38-cutoff10percentMissingData/firstIteration"

inputFileUniqueGenes <- "R-38-M32178-LessThan10percentMissing-postLepMAP2-allChr-uniqueGenes-withPercentMissingData.csv"

inputFileAllMarkers <- "R-38-M32178-LessThan10percentMissing.csv"

outputFile <- "R-38-M32178-LessThan10percentMissing-postLepMAP2-allChr-uniqueGenesWithAllMarkers.csv"


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

inputWithUniqueGenes <- read.csv(paste(path.expand(path),inputFileUniqueGenes,sep="/"),header=TRUE, na.strings=myNaString, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

inputWithAllMarkers <- read.csv(paste(path.expand(path),inputFileAllMarkers,sep="/"),header=TRUE, na.strings=myNaString, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

geneNames <- inputWithUniqueGenes[,"gene"]
geneNamesConcatenated <- paste(geneNames, collapse="|")

colnames(inputWithAllMarkers)[1] <- "marker_id"

# create a new column in inputWithAllMarkers that just has the gene names, derived from the marker_id column.
# library(stringr) required for str_split_fixed:
gene <- str_split_fixed(inputWithAllMarkers$marker_id,separatorForGeneNameVsMarkerId,2)[ ,1]
inputWithAllMarkers <- cbind(gene, inputWithAllMarkers)

rowsToKeep <- logical(nrow(inputWithAllMarkers)) # initializes all to FALSE

for (i in 1:nrow(inputWithAllMarkers)) {
  curGene <- inputWithAllMarkers[i,"gene"]

  for (k in 1:length(geneNames)) {
    if (curGene == geneNames[k]) {
      rowsToKeep[i] <- TRUE
      break
    }
  }
}

output <- inputWithAllMarkers[rowsToKeep, ]


write.csv(output, paste(path.expand(path), outputFile,sep="/"), row.names=FALSE, na=myNaString)

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())



