options(stringsAsFactors = FALSE, warn = 1)
options(scipen = 999) # disable scientific notation for large numbers

args <- commandArgs(trailingOnly = TRUE)

write("Running findBadSamples.R.", stdout())

####################################################
# Input Parameters

path <- "/home/benrancourt/Desktop/junjun/reports-Conifer_192gDNA/reports-fromSplittingTheRefIntoGroups/reports-Conifer_192gDNA-fromSplittingRefIntoGroups"

junjunsFileName <- "SNP64553-Similarity.csv"

chiReportName <- "MAF_cutoff_report_chi.csv"

outputSampleQualityFileName <- "sample_quality.csv"


####################################################
# Execution

jjReport <- read.csv(paste(path.expand(path),junjunsFileName,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

chiReport <- read.csv(paste(path.expand(path),chiReportName,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.



# ----------
# filter jjReport by deleting all rows whose values (gene names) in Locus1 and Locus2 columns don't match

jjRowsToKeep <- logical(length=nrow(jjReport)) # initializes all elements to FALSE

for(rowNum in 1:nrow(jjReport)) {
  geneId1 <- jjReport[rowNum, "Locus1"]
  geneId1 <- strsplit(geneId1,"-",fixed=TRUE)[[1]]
  geneId1 <- paste0(geneId1[1:(length(geneId1)-1)],collapse="-") # exclude the POS from the comparison
  geneId2 <- jjReport[rowNum, "Locus2"]
  geneId2 <- strsplit(geneId2,"-",fixed=TRUE)[[1]]
  geneId2 <- paste0(geneId2[1:(length(geneId2)-1)],collapse="-") # exclude the POS from the comparison
  if ( geneId1 == geneId2 ) {
    # the row passed the filtering requirements
    jjRowsToKeep[rowNum] <- TRUE
  }
}

jjReport <- jjReport[jjRowsToKeep, ]

# ----------
# create a list of sample ids (sample column names in the chi report), with numeric value
# associated with each to use as a tally, initialized to 0.

sampleNames <- colnames(chiReport)[-c(1:7)] # the first 7 columns are meta columns, not samples
tally <- numeric(length=length(sampleNames))
names(tally) <- sampleNames

percentMatches <- numeric(length=nrow(jjReport))
numUnMatched <- numeric(length=nrow(jjReport))
locus1_H <- character(length=nrow(jjReport))
locus1_A <- character(length=nrow(jjReport))
locus2_H <- character(length=nrow(jjReport))
locus2_A <- character(length=nrow(jjReport))
locus_mismatches <- character(length=nrow(jjReport))

# ----------
# for each row in jjReport:
#   - find the Locus1 and Locus2 gene ids (one row each) in the chi report
#   - first, make sure that the H type and A type in both rows match each other, and aren't opposite
#   - if they are swapped (H in one == A in two), then take that into account
#   - compare the two rows cell by cell, and record which sample ids don't match in the tally list

numSamples <- ncol(chiReport) - 7

for (rowNum in 1:nrow(jjReport)) {
  id1 <- strsplit(jjReport[rowNum, "Locus1"],"-",fixed=TRUE)[[1]]
  geneId1 <- paste0(id1[1:(length(id1)-1)],collapse="-")
  pos1 <- id1[length(id1)]

  id2 <- strsplit(jjReport[rowNum, "Locus2"],"-",fixed=TRUE)[[1]]
  geneId2 <- paste0(id2[1:(length(id2)-1)],collapse="-")
  pos2 <- id2[length(id2)]
  
  locus1 <- chiReport[chiReport$CHROM == geneId1 & chiReport$POS == pos1, ]
  locus2 <- chiReport[chiReport$CHROM == geneId2 & chiReport$POS == pos2, ]
  
  if (nrow(locus1) == 1 && nrow(locus2) == 1) {
    locus1_H[rowNum] <- locus1[1, "H_Type"]
    locus1_A[rowNum] <- locus1[1, "A_Type"]
    locus2_H[rowNum] <- locus2[1, "H_Type"]
    locus2_A[rowNum] <- locus2[1, "A_Type"]
    
    row1H <- "H"
    row1A <- "A"
    row2H <- "H"
    row2A <- "A"
#    if (locus1_H[rowNum] == locus2_A[rowNum] && locus1_A[rowNum] == locus2_H[rowNum]) {
#      row2H <- "A"
#      row2A <- "H"
#    }

    if ( (locus1_H[rowNum] == locus2_H[rowNum] || locus1_H[rowNum] == locus2_A[rowNum]) &&
         (locus1_A[rowNum] == locus2_A[rowNum] || locus1_A[rowNum] == locus2_H[rowNum]) ) {
      locus_mismatches[rowNum] <- "-"
    }
    else {
      # there is at least one irreconcileable locus mismatch
      locus_mismatches[rowNum] <- "!!!"
    }
    
    numUnmatchedSamples <- 0

    for (colNum in 8:ncol(locus1)) {
      match <- TRUE
      if ( ((locus1[1,colNum] == row1H || locus1[1,colNum] == "-") && 
           !(locus2[1,colNum] == row2H || locus2[1,colNum] == "-")) ) {
        match <- FALSE
      }
      else if ( ((locus1[1,colNum] == row1A || locus1[1,colNum] == "-") && 
                !(locus2[1,colNum] == row2A || locus2[1,colNum] == "-")) ) {
        match <- FALSE
      }

      if (match == FALSE) {
        numUnmatchedSamples <- numUnmatchedSamples + 1
        sample <- colnames(locus1)[colNum]
        tally[sample] <- tally[sample] + 1
      }
    } # end inner for-loop

    #percentMatched <- (numSamples - numUnmatchedSamples) / numSamples
    numUnMatched[rowNum] <- numUnmatchedSamples
    percentMatches[rowNum] <- (numSamples - numUnmatchedSamples) / numSamples

#    if (percentMatched != jjReport[rowNum,"Similarity"]) {
#      write(paste0("!!!!! Percent matched does not equal Similarity for ",jjReport[rowNum, "Locus1"]," and ",jjReport[rowNum, "Locus2"],"."), stdout())
#    }

  }
  else {
    write(paste0("!!!!! Did not find exactly one row matching ",geneId1," ",pos1," or ",geneId2," ",pos2,"."), stdout())
  }
} # end outer for-loop


tallyDataFrame <- data.frame(tally)
colnames(tallyDataFrame) <- "num_mismatches"
write.csv(tallyDataFrame, paste(path.expand(path), outputSampleQualityFileName,sep="/"), row.names=TRUE)

jjReport <- cbind(jjReport,numUnMatched,percentMatches,locus1_H,locus1_A,locus2_H,locus2_A,locus_mismatches)
filteredJJReportName <- strsplit(junjunsFileName,".csv",fixed=TRUE)[[1]][1]
filteredJJReportName <- paste0(filteredJJReportName,"-filtered.csv")
write.csv(jjReport, paste(path.expand(path), filteredJJReportName,sep="/"), row.names=FALSE)


write("========================FINISHED.", stdout())
