


options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)


write("Running filter42sampleDataset.R.", stdout())


# ########################################################
# Input Parameters:

path <- "/home/benrancourt/Desktop/junjun/reports-42-samples-corrected-without1000rowsfromP35"

inputMafFile <- "MAF_cutoff_report.csv"
inputSnpsFile <- "percentage_snps.csv"

outputMafFile <- "MAF_cutoff_report_filtered.csv"
outputSnpsFile <- "percentage_snps_filtered.csv"

numMetaColumnsMaf <- 3
numSamples <- 42

# a number indicating how many of the second SNP (if there's more than one) are allowed
#numAltSnpAllowed <- 2

##############################################################
# Jun-Jun's instructions:
# Filter data from file ‘MAF_cutoff_report.csv’ (~250K rows) using setting: 
#   MAF >= 5%, 
#   missing data <10%, 
#   no indel 
# (should be ~30,292 rows). 
# Add a additional column as SNP-ID by ‘Contig235_752_T_C/T’.

# Steps:
#   - get true/false column for inputSnps indicating which rows to keep, based on above filtering requirements
#   - filter inputSnps by that true/false column
#   - create new column in inputSnps containing new SNP_ID per above instructions
#   - extract the CHROM and SNP_ID columns from the inputSnps data
#   - merge that with the inputMaf data to filter the inputMaf data.
#   - write both inputSnps and inputMaf to file

inputMaf <- read.csv(paste(path.expand(path),inputMafFile,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

inputSnps <- read.csv(paste(path.expand(path),inputSnpsFile,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.


# ----------
# filter inputSnps

snpsRowsToKeep <- logical(length=nrow(inputSnps)) # initializes all elements to FALSE

for(rowNum in 1:nrow(inputSnps)) {
  if ( inputSnps[rowNum, "MAF"] >= 0.05 
      && inputSnps[rowNum, "empty"] < (0.1*numSamples)
      && nchar(inputSnps[rowNum, "REF"])  == 1 ) {
    # the row passed the filtering requirements
    snpsRowsToKeep[rowNum] <- TRUE
  }
}

inputSnps <- inputSnps[snpsRowsToKeep, ]

# ----------
# create new snp_id column in inputSnps

bases <- inputSnps[ , c(4:7)] # just the columns named: A, C, T, G
basesColNames <- colnames(bases)

snp_id <- character(length=nrow(inputSnps))

for(rowNum in 1:nrow(inputSnps)) {
  snpIdFirstPart <- paste(inputSnps[rowNum, "CHROM"], inputSnps[rowNum, "POS"], inputSnps[rowNum, "REF"], sep="_" )
  
  firstMaxIndex <- which.max(bases[rowNum,])
  firstMax <- basesColNames[firstMaxIndex]
  remainingBaseColNames <- basesColNames[-firstMaxIndex]
  secondMax <- remainingBaseColNames[ which.max(bases[rowNum,-firstMaxIndex]) ]
  snpIdSecondPart <- paste0(firstMax,secondMax)

  snp_id[rowNum] <- paste(snpIdFirstPart, snpIdSecondPart, sep="_")
}

inputSnps <- cbind(snp_id, inputSnps)

# ----------
# filter inputMaf by inputSnps and add snp_id column to it
  
toMerge <- inputSnps[ ,c(1:4)]
inputMaf <- merge(toMerge, inputMaf, by=c("CHROM","POS","REF"), all=FALSE)


write.csv(inputMaf, paste(path.expand(path), outputMafFile,sep="/"), row.names=FALSE)
write.csv(inputSnps, paste(path.expand(path), outputSnpsFile,sep="/"), row.names=FALSE)

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())



