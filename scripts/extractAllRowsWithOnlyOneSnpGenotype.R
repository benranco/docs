


options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)


write("Running extractAllRowsWithOnlyOneSnpGenotype.R.", stdout())


# ########################################################
# Input Parameters:

path <- "/home/benrancourt/Desktop/junjun/region45-SNPpipelineOct2018/reports/LepMAP2/extractSingleSnpRows"

#inputFile <- "R-38-M32178-38329_p0.01-postLepMAP2-lod10and6-allChr.csv"
#numMetaColumns <- 5
#outputFile <- "R-38-M32178-38329_p0.01-postLepMAP2-lod10and6-allChr-withPercentMissingData.csv"

inputFile <- "R45-for-Cr4-co-seg-extraction.csv"
numMetaColumns <- 6
outputFileMostlyUniformSnpRows <- "R45-for-Cr4-co-seg-extraction-mostlyUniformSnpRows.csv"
outputFileAllOtherRows <- "R45-for-Cr4-co-seg-extraction-allOtherRows.csv"

# data cells in the input .csv files (both fathers and seeds) which contain these 
# strings will be interpreted as NA/missing data:
myNaString <- c("-")

# a number indicating how many of the second SNP (if there's more than one) are allowed
numAltSnpAllowed <- 2

##############################################################


input <- read.csv(paste(path.expand(path),inputFile,sep="/"),header=TRUE, na.strings=myNaString, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

firstDataCol <- numMetaColumns + 1

rowsWithOneSnpType <- logical(length=nrow(input)) # initializes all elements to FALSE


for(rowNum in 1:nrow(input)) {
  rowData <- as.matrix(input[rowNum,firstDataCol:ncol(input)]) 
  # count and tabulate the SNP types in the row, in descending order of frequency (should be only one or two types)
  tabulatedRowData <- sort(table(rowData, useNA = "no"),decreasing=TRUE)
  tabulatedRowDataNames <- names(tabulatedRowData)

  # if there's only one type of non-NA data
  if( length(tabulatedRowDataNames) == 1) {
    rowsWithOneSnpType[rowNum] <- TRUE
  }
  else if( length(tabulatedRowDataNames) > 1 && tabulatedRowData[2] <= numAltSnpAllowed) {
    rowsWithOneSnpType[rowNum] <- TRUE
  }
}

onlyRowsWithOneSnp <- input[rowsWithOneSnpType, ]
allOtherRows <- input[!rowsWithOneSnpType, ]

write.csv(onlyRowsWithOneSnp, paste(path.expand(path), outputFileMostlyUniformSnpRows,sep="/"), row.names=FALSE, na=myNaString)
write.csv(allOtherRows, paste(path.expand(path), outputFileAllOtherRows,sep="/"), row.names=FALSE, na=myNaString)

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())



