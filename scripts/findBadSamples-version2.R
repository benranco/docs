options(stringsAsFactors = FALSE, warn = 1)
options(scipen = 999) # disable scientific notation for large numbers

args <- commandArgs(trailingOnly = TRUE)

write("Running findBadSamples-version2.R.", stdout())

####################################################
# Input Parameters

path <- "/home/benrancourt/Desktop/junjun/reports-Conifer_192gDNA/reports-fromSplittingTheRefIntoGroups/reports-Conifer_192gDNA-fromSplittingRefIntoGroups"

junjunsFileName <- "G1-7_ML.csv"

chiReportName <- "MAF_cutoff_report_chi.csv"

outputFilteredReportFileNamePostFix <- "-filtered-run4.csv"
outputSampleQualityFileNamePostFix <- "-sample_quality-run4.csv"


####################################################
# Execution

jjReportMaster <- read.csv(paste(path.expand(path),junjunsFileName,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

chiReportMaster <- read.csv(paste(path.expand(path),chiReportName,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.


################################
partOne <- function(jjReport, chiReport) {

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
  write(paste0("Number of rows to keep: ",nrow(jjReport)), stdout())
  # ----------

  distance <- numeric(length=nrow(jjReport))
  coupling_or_repulsion <- character(length=nrow(jjReport))
  num_missing <- numeric(length=nrow(jjReport))
  num_matches <- numeric(length=nrow(jjReport))
  num_mismatches <- numeric(length=nrow(jjReport))
  mismatch_percent <- numeric(length=nrow(jjReport))
  mismatched_sampleIdList <- character(length=nrow(jjReport))

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
      
      distance[rowNum] <- abs(as.numeric(pos1) - as.numeric(pos2))

      nMissing <- 0
      nMatched <- 0
      nMismatched <- 0
      mismatchSampleIds <- ""

      for (colNum in 8:ncol(locus1)) {
        if (locus1[1,colNum] == "-" || locus2[1,colNum] == "-") {
          nMissing <- nMissing + 1
        }
        else {
          # only do these tests if neither has missing data

          match <- TRUE
        
          if ( locus1[1,colNum] == row1H && locus2[1,colNum] != row2H ) {
            match <- FALSE
          }
          else if ( locus1[1,colNum] == row1A && locus2[1,colNum] != row2A ) {
            match <- FALSE
          }

          if (match == TRUE) {
            nMatched <- nMatched + 1
          } 
          else {
            # match is false
            nMismatched <- nMismatched + 1          
            sampleName <- colnames(locus1)[colNum]
            sampleName <- strsplit(sampleName,".",fixed=TRUE)[[1]]
            sampleName <- sampleName[length(sampleName)-1]
            mismatchSampleIds <- paste(mismatchSampleIds,sampleName,sep=";")
          }
        }
      } # end inner for-loop

      if (nMatched >= nMismatched) {
        coupling_or_repulsion[rowNum] <- "C"
      }
      else {
        coupling_or_repulsion[rowNum] <- "R"
        # treat the mismatched as matched and vice-versa:
        tempMismatched <- nMismatched
        nMismatched <- nMatched
        nMatched <- tempMismatched
      }

      num_missing[rowNum] <- nMissing
      num_matches[rowNum] <- nMatched
      num_mismatches[rowNum] <- nMismatched  
      mismatch_percent[rowNum] <- nMismatched / (nMismatched + nMatched)
      mismatched_sampleIdList[rowNum] <- mismatchSampleIds

    }
    else {
      write(paste0("!!!!! Did not find exactly one row matching ",geneId1," ",pos1," or ",geneId2," ",pos2,"."), stdout())
    }
  } # end outer for-loop


  jjReport <- cbind(jjReport, distance, num_missing, num_matches, num_mismatches, mismatch_percent, coupling_or_repulsion, mismatched_sampleIdList, locus1_H,locus1_A,locus2_H,locus2_A,locus_mismatches)

  filteredJJReportName <- strsplit(junjunsFileName,".csv",fixed=TRUE)[[1]][1]
  filteredJJReportName <- paste0(filteredJJReportName,outputFilteredReportFileNamePostFix)
  write.csv(jjReport, paste(path.expand(path), filteredJJReportName,sep="/"), row.names=FALSE)
  
  jjReport
} # end function partOne

################################


################################
partTwo <- function(jjReport, chiReport) {

  # ----------
  # create a list of sample ids (sample column names in the chi report), with numeric value
  # associated with each to use as a tally, initialized to 0.

  sampleNames <- colnames(chiReport)[-c(1:7)] # the first 7 columns are meta columns, not samples
  for (i in 1:length(sampleNames)) {
    name <- strsplit(sampleNames[i],".",fixed=TRUE)[[1]]
    name <- name[length(name)-1]
    sampleNames[i] <- name
  }

  tally_no_filter <- numeric(length=length(sampleNames))
  names(tally_no_filter) <- sampleNames

  tally_mismatch_max10percent <- numeric(length=length(sampleNames))
  names(tally_mismatch_max10percent) <- sampleNames

  tally_lod_limit5 <- numeric(length=length(sampleNames))
  names(tally_lod_limit5) <- sampleNames

  tally_distance_max100 <- numeric(length=length(sampleNames))
  names(tally_distance_max100) <- sampleNames

  tally_all_filter <- numeric(length=length(sampleNames))
  names(tally_all_filter) <- sampleNames

  # ----------
  for (rowNum in 1:nrow(jjReport)) {
    mismatchedSamples <- jjReport[rowNum,"mismatched_sampleIdList"]
    mismatchedSamples <- strsplit(mismatchedSamples,";",fixed=TRUE)[[1]]
    
    if (length(mismatchedSamples) > 0) {
      for (i in 1:length(mismatchedSamples)) {
        sample <- mismatchedSamples[i]
        if (nchar(sample) > 0) {
          tally_no_filter[sample] <- tally_no_filter[sample] + 1
          
          meetsMismatch <- FALSE
          meetsLod <- FALSE
          meetsDist <- FALSE

          if (jjReport[rowNum,"mismatch_percent"] <= 0.1) {
            tally_mismatch_max10percent[sample] <- tally_mismatch_max10percent[sample] + 1
            meetsMismatch <- TRUE
          }
          if (jjReport[rowNum,"LOD"] >= 5) {
            tally_lod_limit5[sample] <- tally_lod_limit5[sample] + 1
            meetsLod <- TRUE
          }
          if (jjReport[rowNum,"distance"] <= 100) {
            tally_distance_max100[sample] <- tally_distance_max100[sample] + 1
            meetsDist <- TRUE
          }
          if (meetsMismatch && meetsLod && meetsDist) {
            tally_all_filter[sample] <- tally_all_filter[sample] + 1
          }
        }
      } # end inner for-loop
    }

  } # end outer for-loop

  tallyDataFrame <- data.frame(tally_no_filter)
  tallyDataFrame <- cbind(tallyDataFrame, tally_mismatch_max10percent, tally_lod_limit5, tally_distance_max100, tally_all_filter)

  filename <- paste0( strsplit(junjunsFileName,".csv")[[1]][1], outputSampleQualityFileNamePostFix)
  write.csv(tallyDataFrame, paste(path.expand(path), filename,sep="/"), row.names=TRUE)

} # end function partTwo

################################

jjReportMaster <- partOne(jjReportMaster, chiReportMaster)
write("========================FINISHED PART 1.", stdout())

# for testing purposes:
#jjReportMaster <- read.csv(paste(path.expand(path),"G1-7_ML-filtered-run3.csv",sep="/"),header=TRUE, check.names=FALSE)

partTwo(jjReportMaster, chiReportMaster)
write("========================FINISHED PART 2.", stdout())


