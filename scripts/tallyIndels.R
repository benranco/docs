


options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)


write("Running tallyIndels.R.", stdout())

##########################################################
# This script basically does this:
#
#- read in MAF table
#- keep all rows whose REF value's length > 1
#- print out how many that is
#- reformat the alleles in the data columns to remove the "/" at the end 
#- for each row:
#  - save the alleles of H and A to new columns
#  - save the totals of H and A to new columns  
#  - change the most common allele to H
#  - change second most common allele to A 
#  - (if the most common and second most common allele are tied, then the one the matches the REF is set to H)
#  - (change any additional remaining alleles to NA)
#- write the final outout to .csv file (outputting all NA's as -)
#
##########################################################

# ########################################################
# Input Parameters:

#path <- "/home/benrancourt/Desktop/junjun/region45-SNPpipelineNov2018-withIndels-withFupansRef/reports"
path <- args[1]
reportsSubDir <- "reports"

inputFile <- "MAF_cutoff_report.csv"

outputFile <- "MAF_cutoff_report_indels-keepingAdditionalAlleles.csv"
alternateOutputFile <- "MAF_cutoff_report_indels.csv"

# data cells in the input .csv files (both fathers and seeds) which contain these 
# strings will be interpreted as NA/missing data. 
inputNaStrings <- c("NA","-")

# All NA values will be represented as this string in the output file:
outputNaString <- c("-")

# a number indicating how many of the second SNP (if there's more than one) are allowed
numAltSnpAllowed <- 2

##############################################################

# function to get rid of the last character of the given string
dropLastCharacterOfString <- function(x) {
  substring(x,1,nchar(x)-1)
}

# function to replace a given string x if it matches string toMatch
replaceMatchingString <- function(x, toMatch, replacement) {
  if (!is.na(x) && x == toMatch) {
    x <- replacement
  }
  x
}

# function to replace a given string x with replacement if it doesn't match toMatch1 or toMatch2
replaceNonMatchingString <- function(x, toMatch1, toMatch2, replacement) {
  if (!is.na(x) && x != toMatch1 && x != toMatch2) {
    x <- replacement
  }
  x
}


mafReport <- read.csv(paste(path.expand(path),reportsSubDir,inputFile,sep="/"),header=TRUE, na.strings=inputNaStrings, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

# check if first col is not CHROM, and remove it if so (it's just a left-over old index column if it exists)
if (colnames(mafReport)[1] != "CHROM")
{
  mafReport <- mafReport[,-1]
}

firstDataCol <- 4
# check if COMBINED col exists, and adjust starting col num 
if (colnames(mafReport)[4] == "COMBINED")
{
  firstDataCol <- 5
}


# keep only those rows whose REF value has more than one character in it (i.e. it's an indel):
mafReport <- mafReport[ (nchar(mafReport$REF) > 1), ]

write(paste0("Number of indels: ",nrow(mafReport)), stdout())

# initialize new data vectors:
H_type <- character(nrow(mafReport))
A_type <- character(nrow(mafReport))
H_total <- numeric(nrow(mafReport))
A_total <- numeric(nrow(mafReport))

dataAsMatrix <- as.matrix(mafReport[ ,firstDataCol:ncol(mafReport)]) 


for(rowNum in 1:nrow(dataAsMatrix)) {

  # reformat the data values: get rid of the "/" at the end of each allele in the data:
  dataAsMatrix[rowNum, ] <- sapply(dataAsMatrix[rowNum, ], dropLastCharacterOfString, USE.NAMES=FALSE)
  
  # count and tabulate the SNP types in the row, in descending order of frequency (should be only one or two types)
  tabulatedRowData <- sort(table(dataAsMatrix[rowNum, ], useNA = "no"),decreasing=TRUE)
  tabulatedRowDataNames <- names(tabulatedRowData)
 
  if( length(tabulatedRowDataNames) >= 1) {
    H_type[rowNum] <- tabulatedRowDataNames[1]
    H_total[rowNum] <- tabulatedRowData[1]
  }
  if( length(tabulatedRowDataNames) > 1) {
    # if the second most common allele is less frequent than the first
    if( tabulatedRowData[2] < tabulatedRowData[1]) {
      A_type[rowNum] <- tabulatedRowDataNames[2]
      A_total[rowNum] <- tabulatedRowData[2]
    }
    # if the second most common allele is tied with the first, we set whichever allele matches the
    # Ref to H and the other to A:
    else if( tabulatedRowData[2] == tabulatedRowData[1]) {
      if( mafReport[rowNum,'REF'] == tabulatedRowDataNames[1]) {
        A_type[rowNum] <- tabulatedRowDataNames[2]
        A_total[rowNum] <- tabulatedRowData[2]
      }
      else {
        A_type[rowNum] <- tabulatedRowDataNames[1]
        A_total[rowNum] <- tabulatedRowData[1]
        H_type[rowNum] <- tabulatedRowDataNames[2]
        H_total[rowNum] <- tabulatedRowData[2]
      }
    }
  }

  # replace all occurences of the H allele in the row with "H":
  dataAsMatrix[rowNum, ] <- sapply(dataAsMatrix[rowNum, ], replaceMatchingString, toMatch=H_type[rowNum], replacement="H", USE.NAMES=FALSE)
  # replace all occurences of the A allele in the row with "A":
  dataAsMatrix[rowNum, ] <- sapply(dataAsMatrix[rowNum, ], replaceMatchingString, toMatch=A_type[rowNum], replacement="A", USE.NAMES=FALSE)

} # end for-loop

# make a second copy of the data for the next loop:
alternateDataAsMatrix <- dataAsMatrix

for(rowNum in 1:nrow(alternateDataAsMatrix)) {
  # replace all other data that doesn't match either of the H or the A alleles in the row with NA:
  alternateDataAsMatrix[rowNum, ] <- sapply(alternateDataAsMatrix[rowNum, ], replaceNonMatchingString, toMatch1="H", toMatch2="A", replacement=NA, USE.NAMES=FALSE)
}

# put together our new output data.frame:
output <- cbind(mafReport[,1:3], H_type, A_type, H_total, A_total, dataAsMatrix)
alternateOutput <- cbind(mafReport[,1:3], H_type, A_type, H_total, A_total, alternateDataAsMatrix)

write.csv(output, paste(path.expand(path), reportsSubDir, outputFile,sep="/"), row.names=FALSE, na=outputNaString)
write.csv(alternateOutput, paste(path.expand(path), reportsSubDir, alternateOutputFile,sep="/"), row.names=FALSE, na=outputNaString)

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())



