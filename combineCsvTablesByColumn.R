


options(stringsAsFactors = FALSE, warn = 1)

write("Running convertCsvToCircosInput.R.", stdout())


# ########################################################
# Input parameters:

path <- "/home/benrancourt/Downloads/combine"

inputMAFFile <- "MAF_cutoff_report.csv"
inputPercentSNPFile <- "percentage_snps.csv"
inputHollyFile <- "5x11_Df_families_aa_changes_vs_3345_ref-35-65%_freq-modifiedMappingName.csv"

fullOutputFileName <- "full_output.csv"

reducedOutputFileName <- "reduced_output.csv"

# ########################################################

# using check.names=FALSE for these tables in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.
inputMAF <- read.csv(paste(path.expand(path),inputMAFFile,sep="/"),header=TRUE, check.names=FALSE)
inputPercentSNP <- read.csv(paste(path.expand(path),inputPercentSNPFile,sep="/"),header=TRUE, check.names=FALSE)
inputHolly <- read.csv(paste(path.expand(path),inputHollyFile,sep="/"),header=TRUE, check.names=FALSE)

# - create a column for the allele in inputPercentSNP (and figure out what the allele is)
# - rename columns in inputHolly that are used in the key, so they are the same as the column names in inputPercentSNP
# - remove the columns from inputHolly and inputPercentSNP that won't be kept
# - merge inputHolly and inputPercentSNP by their shared key
# - merge the mergedData with inputMAF, and write result to .csv
# - create a copy that removes all rows with NA in the inputHolly columns, and write result to .csv


# - create a column for the allele in inputPercentSNP (and figure out what the allele is)

Allele <- numeric(nrow(inputPercentSNP))

for (i in 1:nrow(inputPercentSNP))
{
  # find max and second max; if max is not REF, then it is allele, else second max is allele

  indexOfMax <- which(inputPercentSNP[i,c("A","C","T","G")]==inputPercentSNP[i,"max"])
  nameOfMax <- colnames(inputPercentSNP[i,c("A","C","T","G")])[indexOfMax[1]]
  if (length(indexOfMax) > 1)
  { 
    # if there's a tie for max and second max
    nameOfSecondMax <- colnames(inputPercentSNP[i,c("A","C","T","G")])[indexOfMax[2]]  
  }
  else
  {
    indexOfSecondMax <- which(inputPercentSNP[i,c("A","C","T","G")]==inputPercentSNP[i,"max"])
    nameOfSecondMax <- colnames(inputPercentSNP[i,c("A","C","T","G")])[indexOfSecondMax[1]]
  }

  if (inputPercentSNP[i,"REF"] != nameOfMax)
  {
    Allele[i] <- nameOfMax
  }
  else
  {
    Allele[i] <- nameOfSecondMax
  }
}

inputPercentSNP <- cbind(inputPercentSNP, Allele)

# - rename columns in inputHolly that are used in the key, so they are the same as the column names in inputPercentSNP
colnames(inputHolly)[grep("Mapping", colnames(inputHolly))] <- "CHROM"
colnames(inputHolly)[grep("Position", colnames(inputHolly))] <- "POS" # "Position" is found only in "Reference.Position"
colnames(inputHolly)[grep("Reference", colnames(inputHolly))] <- "REF"

# - remove the columns from inputHolly and inputPercentSNP that won't be kept
inputPercentSNP <- inputPercentSNP[,-which(names(inputPercentSNP) %in% c("A","C","T","G","empty"))]
inputHolly <- inputHolly[,-which(names(inputHolly) %in% c("Type","Length","T","G","empty"))]

# - merge inputHolly and inputPercentSNP by their shared key
mergedData <- merge(x=inputPercentSNP, y=inputHolly, all.x=TRUE, by=c("CHROM","POS","REF","Allele"))

# - merge the mergedData with inputMAF, and write result to .csv
mergedData <- merge(x=inputMAF, y=mergedData, all.x=TRUE, by=c("CHROM","POS","REF"))
write.table(mergedData, file= paste(path.expand(path),fullOutputFileName,sep="/"), append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

# - create a copy that removes all rows with NA in the inputHolly columns, and write result to .csv
reducedData <- mergedData[ !(is.na(mergedData$Frequency)), ]
write.table(reducedData, file= paste(path.expand(path),reducedOutputFileName,sep="/"), append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())



