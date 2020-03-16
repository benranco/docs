############################################################################################
# This script converts a .csv report from our SNPpipeline (such as filled_report.csv or 
# MAF_cutoff_report.csv or any other report that follows this format) file from our SNPpipeline 
# output into PLINK format, according to the information at: 
# http://zzz.bwh.harvard.edu/plink/data.shtml
#
# Before converting the csv report to PLINK format, it first removes all rows that have more
# than two alleles, because PLINK can only handle biallelic data.
#
# The specific type of PLINK format files it creates is "transposed fileset" (.tped and .tfam), 
# which presents the same data as the default PLINK .ped and .map files but organised in a 
# way which is closer to the structure of our MAF_cutoff_report.csv.
# 
# To convert our MAF_cutoff_report to (transposed PLINK format) .tped and .tfam format: 
# 
# For the .tped file, for each row in our MAF_cutoff_report:
# col1 - add an initial chromosome column, set to 1 (it needs to be 1-22 for it to not be 
#        excluded by KING as non-autosomal)
# col2 - the input Maf report's CHROM column concatenated with the POS and REF columns is used 
#        for the rs# or snp identifier
# col3 - add a new Genetic distance column (set to 0 for all rows for our purposes)
# col4 - the POS column is used as-is for the Base-pair position (bp units)
# col5->end - the remaining genetic data makes up the rest of the row, but converted according 
# to these examples: 
#     "A/" becomes "A A", "A/A" becomes "A A", "A/C" becomes "A C", NA becomes "0 0".
# 
# For the .tfam file, for each sample column, create a row with six columns:
# col1 - family name = whatever we decide, the same for all rows
# col2 - individual id = the name of the column (minus the .tab)
# col3 - paternal id = 0
# col4 - maternal id = 0 or "mom"
# col5 - sex = 0 (for unknown)
# col6 - phenotype = -9 (for missing)
# 
# Once these .tped and .tfam PLINK files have been created, this script then uses PLINK to convert 
# them to binary PLINK format. The PINK command it uses is of this form:
# plink --tped inputfile.tped --tfam inputfile.tfam --make-bed --out outputfile
# 
# Info on the KING software can be found here:
# http://people.virginia.edu/~wc9c/KING/
# 
############################################################################################
# Input Parameters:
args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE, warn = 1)


# path: The directory path of the folder that the input .csv file is in. This is also the 
# folder that the output .linkage file will be saved to. 
#path <- "/home/benrancourt/Desktop/junjun/reports-Conifer_192gDNA/KINGdata_MAFreport"
path <- "/home/benrancourt/Desktop/junjun/reports-Conifer_192gDNA/KINGdata_filledReport"

#inputFilename <- "MAF_cutoff_report-first500.csv"
inputFilename <- "filled_report--addedMother-removedNonBiallelicRows.csv"

#outputFilenameRoot <- "MAF_cutoff_report-first500"
outputFilenameRoot <- "filled_report"

# can be anything
familyName <- "192gDNA"

# choose (TRUE or FALSE) whether to search for non-biallelic rows in the input data and remove them.
removeNonBiallelicRows <- FALSE

# choose (TRUE or FALSE) whether to create a new sample called "mom" based on the REF column,
# which will be considered the mother of all the other samples:
createMomSample <- FALSE

# missingDataIndicator: This parameter relates to the data in the original
# .csv data file that you supplied. Some of the data in your original input .csv will
# likely be missing, in which case it might be represented as NA or "-" or maybe even be
# empty, or something like that to indicate that it is missing.
# Indicate what your .csv file uses to represent missing data.
#missingDataIndicator <- "-"
missingDataIndicator <- NA

# full path to the plink executable, which will be used to convert the plink data to binary plink format:
plinkExecutableLocation <- "/home/benrancourt/Desktop/KING/plink-1.07-x86_64/plink"


# a list of samples (column names) to remove from the data:
samplesToRemove <- c("NS.1125.002.IDT_i7_19...IDT_i5_19.MR504.tab", "NS.1125.002.IDT_i7_20...IDT_i5_20.MR507.tab", "NS.1125.002.IDT_i7_43...IDT_i5_43.MR843.tab", "NS.1125.002.IDT_i7_67...IDT_i5_67.MR1034.tab", "NS.1125.003.IDT_i7_100...IDT_i5_100.MS130.tab", "NS.1125.003.IDT_i7_109...IDT_i5_109.MS301.tab", "NS.1125.003.IDT_i7_113...IDT_i5_113.MS331.tab", "NS.1125.003.IDT_i7_125...IDT_i5_125.MS450.tab", "NS.1125.003.IDT_i7_140...IDT_i5_140.MS540.tab", "NS.1125.003.IDT_i7_148...IDT_i5_148.MS591.tab", "NS.1125.003.IDT_i7_157...IDT_i5_157.MS621.tab", "NS.1125.003.IDT_i7_158...IDT_i5_158.MS622.tab", "NS.1125.003.IDT_i7_161...IDT_i5_161.MS633.tab", "NS.1125.003.IDT_i7_170...IDT_i5_170.MS665.tab")
# OLD samples to remove:
##samplesToRemove <- c("NS.1125.002.IDT_i7_19...IDT_i5_19.MR504.tab", "NS.1125.002.IDT_i7_43...IDT_i5_43.MR843.tab", "NS.1125.002.IDT_i7_67...IDT_i5_67.MR1034.tab", "NS.1125.003.IDT_i7_100...IDT_i5_100.MS130.tab", "NS.1125.003.IDT_i7_113...IDT_i5_113.MS331.tab", "NS.1125.003.IDT_i7_157...IDT_i5_157.MS621.tab", "NS.1125.003.IDT_i7_158...IDT_i5_158.MS622.tab", "NS.1125.003.IDT_i7_161...IDT_i5_161.MS633.tab")

#########################################################################
# Function borrowed from report_gen_p2.R of the SNPpipeline:
# Function: isRowNonBiallelic - determines if the given row data has more than two alleles.
# input:  a list of all the unique data items in the report. Each data item is formatted:
#         "allele/allele" or possibly (for legacy purposes) "allele/". An allele is usually
#         only one character long, but in the case where indels were not removed in the first
#         part of the pipeline, an allele can have two or more characters.
# returns: TRUE if the row has more than two alleles, FALSE if it only has two (or one) allele
isRowNonBiallelic <- function(namesOfGenotypesInRow)
{
  alleleIndex <- 1
  individualAlleles <- character(length(namesOfGenotypesInRow) * 2) # because, for example "A/A" has two alleles, "A/G" has two. If a type is of format "A/", this function will pretend it is in format "A/A" and count the single A twice. This is okay to do, because the count is only used internally to this function.
  for ( type in namesOfGenotypesInRow )
  {
    allelesInType <- strsplit(type, "/")[[1]]
    individualAlleles[alleleIndex] <- allelesInType[1]
    alleleIndex <- alleleIndex + 1
    individualAlleles[alleleIndex] <- ifelse (length(allelesInType) > 1, allelesInType[2], allelesInType[1]) # Because for haploid data, it is in format "A/" instead of "A/A".
    alleleIndex <- alleleIndex + 1
  }
  # if there are 3 or more unique alleles in the row, return TRUE, else return FALSE
  return ( ifelse( length(names(table(individualAlleles, useNA = "no"))) >= 3, TRUE, FALSE ) )
}



#########################################################################
# Execution Code:

if (is.na(missingDataIndicator)) {
  missingDataIndicator = "NA"
}

inputReport <- read.csv(paste(path.expand(path),inputFilename,sep="/"), na.strings=missingDataIndicator)

if (exists("samplesToRemove") && !is.null(samplesToRemove) && !is.na(samplesToRemove) && length(samplesToRemove > 0)) {
  message("removing these samples (columns) from the report: ")
  message(paste0(samplesToRemove, collapse=",\n  "))
  # remove the columns identified in samplesToRemove
  inputReport <- inputReport[ , -which(colnames(inputReport) %in% samplesToRemove)]
}

if(!("COMBINED" %in% colnames(inputReport))) {
  startCol <- 4
} else {
  startCol <- 5
}


####
if (removeNonBiallelicRows == TRUE || createMomSample == TRUE) {

  if (removeNonBiallelicRows == TRUE) {
    message("first removing all non-biallelic rows from the report.")
  }
  if (createMomSample == TRUE) {
    message("also creating a mother sample by inference.")  
  }
  mom <- character(nrow(inputReport))

  originalNumRows <- nrow(inputReport)
  rowsToRemove <- logical(originalNumRows) # initializes all to FALSE
  numColsInReport <- ncol(inputReport)

  for (curRow in 1:nrow(inputReport)) {
    rowData <- as.matrix(inputReport[curRow,startCol:numColsInReport]) 
    tabulatedRowData <- table(rowData, useNA = "no")
    tabulatedRowDataNames <- names(tabulatedRowData)

    if (removeNonBiallelicRows == TRUE) {
      # remove row if it has 3 or more unique alleles (eg. a row with A/A,A/C,C/C is kept, but 
      # a row with A/A,A/C,C/T is removed). We only keep biallelic data.
      if( isRowNonBiallelic(tabulatedRowDataNames) ) {
        rowsToRemove[curRow] <- TRUE
      }
    }

    # TODO: There is a weakness in this code section, because it assumes the data is always 
    # homozygous (e.g always in format "A/" or "A/A", never "A/G". But since the data this 
    # script was written for is homozygous, this is safe. If the script is repurposed for 
    # heterozygous input data then we might want to revamp this part to count individual alleles
    # in the row and construct the momValue from that. But this would be far more complex and would
    # involve a major overhaul of the function isRowNonBiallelic to serve this purpose in addition to
    # its original purpse.
    if (createMomSample == TRUE) {
      momValue <- NA
      if (length(tabulatedRowData) == 1) {
        val <- strsplit(tabulatedRowDataNames[1], "/")[[1]][1]
        momValue <- paste(val, val, sep="/")
      } else if (length(tabulatedRowData) > 1) {
        tabulatedRowData <- sort(tabulatedRowData, decreasing=TRUE)
        tabulatedRowDataNames <- names(tabulatedRowData)

        if (tabulatedRowData[1] / (tabulatedRowData[1] + tabulatedRowData[2]) >= 0.7 ) {
          # if the most frequent of the top two values occurs >= 70% of the time, make momValue homozygous
          val <- strsplit(tabulatedRowDataNames[1], "/")[[1]][1]
          momValue <- paste(val, val, sep="/")
        } else {
          # make momValue heterozygous
          val1 <- strsplit(tabulatedRowDataNames[1], "/")[[1]][1]
          val2 <- strsplit(tabulatedRowDataNames[2], "/")[[1]][1]
          momValue <- paste(val1, val2, sep="/")
        }
      }
      
      mom[curRow] <- momValue
    }


  } # end for-loop

  if (removeNonBiallelicRows == TRUE) {
    inputReport <- inputReport[!rowsToRemove, ]
    newNumRows <- nrow(inputReport)
    message(paste0("removed ", (originalNumRows - newNumRows), " non-biallelic rows, ", newNumRows, " rows remaining."))
  }
  if (createMomSample == TRUE) {
    if (removeNonBiallelicRows == TRUE) {      
      mom <- mom[!rowsToRemove]
    }
    inputReport$mom <- mom
  }
}

# if the input report was modified in a way that might take a long time to redo, save a copy of it to file
if (removeNonBiallelicRows == TRUE || createMomSample == TRUE ) {

  editedCsvFileName <- strsplit(inputFilename, ".csv")[[1]][1]
  suffix <- "-"
  if (createMomSample == TRUE) {
    suffix <- paste0(suffix, "-addedMother")
  }
  if (removeNonBiallelicRows == TRUE) {
    suffix <- paste0(suffix, "-removedNonBiallelicRows")
  }
  suffix <- paste0(suffix, ".csv")
  editedCsvFileName <- paste0(editedCsvFileName, suffix)
  message(paste0("saving the edited csv report as ",editedCsvFileName))
  write.csv(inputReport, paste(path.expand(path), editedCsvFileName, sep = "/"), row.names=FALSE)

}

####
message("converting the report to PLINK format.")


# For the .tfam file, for each sample column, create a row with six columns:
# col1 - family name = whatever we decide, the same for all rows
# col2 - individual id = the name of the column (minus the .tab)
# col3 - paternal id = 0
# col4 - maternal id = 0 or "mom"
# col5 - sex = 0 (for unknown)
# col6 - phenotype = -9 (for missing)

# assume the first three columns of the MAF_cutoff_report are meta columns:
sampleNames <- colnames(inputReport)[4:length(colnames(inputReport))]
#We actually want to keep the .tab#sampleNames <- sub(".tab", "", sampleNames, fixed=TRUE)

family <- character(length(sampleNames))
family[1:length(family)] <- familyName

father <- character(length(sampleNames))
father[1:length(father)] <- "0"

mother <- character(length(sampleNames))
if (createMomSample == TRUE) {
  mother[1:length(mother)] <- "mom"
  mother[length(mother)] <- "0" # the last sample column in inputReprot is sample id "mom" created by this script earlier, so it has no mother
} else {
  mother[1:length(mother)] <- "0"
}

sex <- character(length(sampleNames))
sex[1:length(sex)] <- "0"

phenotype <- character(length(sampleNames))
phenotype[1:length(phenotype)] <- "-9"

tfam <- data.frame(family, sampleNames, father, mother, sex, phenotype)

# For the .tped file, for each row in our MAF_cutoff_report:
# col1 - add an initial chromosome column, set to 1 (it needs to be 1-22 for it to not be 
#        excluded by KING as non-autosomal)
# col2 - the input Maf report's CHROM column concatenated with the POS and REF columns is used 
#        for the rs# or snp identifier
# col3 - add a new Genetic distance column (set to 0 for all rows for our purposes)
# col4 - the POS column is used as-is for the Base-pair position (bp units)
# col5->end - the remaining genetic data makes up the rest of the row, but converted according 
# to these examples: 
#     "A/" becomes "A A", "A/A" becomes "A A", "A/C" becomes "A C", NA becomes "0 0".

col1 <- character(nrow(inputReport))
col1[1:length(col1)] <- "1"

snpId <- paste(inputReport[,1], inputReport[,2], inputReport[,3], sep="--")

distance <- character(nrow(inputReport))
distance[1:length(distance)] <- "0"

for (colNum in startCol:ncol(inputReport)) {

  # when doing lots of replacements in a dataframes, it's much faster to operate on the column 
  # separately and then replace the whole column in the dataframe with the updated one.
  curCol <- inputReport[,colNum]

  for (rowNum in 1:length(curCol)) {
    if (!is.na(curCol[rowNum])) {
      alleles <- strsplit(curCol[rowNum], "/")[[1]]
      if (length(alleles)==1) {
        curCol[rowNum] <- paste(alleles[1],alleles[1], sep=" ")
      }
      else if (length(alleles)==2) {
        curCol[rowNum] <- paste(alleles[1],alleles[2], sep=" ")
      }
    }
  } 
  
  inputReport[,colNum] <- curCol
}

inputReport[is.na(inputReport)] <- "0 0"

inputReport <- cbind(col1, snpId, distance, inputReport[,2], inputReport[,4:ncol(inputReport)])


tfamName <- paste0(outputFilenameRoot,".tfam")
tpedName <- paste0(outputFilenameRoot,".tped")

write.table(tfam, file= paste(path.expand(path),tfamName,sep="/"), append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

write.table(inputReport, file= paste(path.expand(path),tpedName,sep="/"), append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# ------

message("done.")

message("Converting to binary plink format. ")

tpedFile <- paste(path, paste0(outputFilenameRoot,".tped"), sep="/")
tfamFile <- paste(path, paste0(outputFilenameRoot,".tfam"), sep="/")
plinkOutFilenameBase <- paste(path, paste0("binary_",outputFilenameRoot), sep="/")

plinklog <- system2(plinkExecutableLocation, c("--tped", tpedFile, "--tfam", tfamFile, "--make-bed", "--out", plinkOutFilenameBase), stdout=TRUE)

message("Done.")
#message(paste0("Done. \n\n",plinklog))

