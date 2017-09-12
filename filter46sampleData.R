options(stringsAsFactors = FALSE, warn = 1)
options(scipen = 999) # disable scientific notation for large numbers

write("Running getDepthStats.R.", stdout())


# ########################################################
# Input Parameters:

args <- commandArgs(trailingOnly = TRUE)

path <- "/home/benrancourt/Desktop/junjun/newSNPpipeline/reports-snppipeline-exom14706ref"

#inputMAFcutoffReportName <- "MAF_cutoff_report_1000rows_tweaked.csv"
inputMAFcutoffReportName <- "MAF_cutoff_report.csv"

#outputMAFcutoffReportName <- "MAF_cutoff_report_1000rows_tweaked_filtered.csv"
outputMAFcutoffReportName <- "MAF_cutoff_report_filtered_no_overlap.csv"

accept10percentOverlap <- FALSE

# ########################################################
# Execution code:


mafReport1 <- read.csv(paste(path.expand(path),inputMAFcutoffReportName,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.)

# check if first col is not CHROM, and remove it if so (it's just a left-over old index column if it exists)
if (colnames(mafReport1)[1] != "CHROM")
{
  mafReport1 <- mafReport1[,-1]
}

startingCol <- 4
# check if COMBINED col exists, and adjust starting col num 
if (colnames(mafReport1)[4] == "COMBINED")
{
  startingCol <- 5
}

rsColNames <- c("RS1.tab", "RS2_rep.tab", "RS3.tab", "RS4.tab", "RS5.tab", "RS6.tab", "RS7.tab", "RS8.tab", "RS9.tab", "RS10.tab")
ssColNames <- c("SS1.tab", "SS2.tab", "SS3.tab", "SS4.tab", "SS5.tab", "SS6.tab", "SS7.tab", "SS8.tab", "SS9.tab", "SS10.tab")
lprsColNames <- c("LP-RS1.tab", "LP-RS8.tab", "LP-RS9.tab")
lpssColNames <- c("LP-SS4.tab", "LP-SS5.tab", "LP-SS6.tab")

rsColumns <- mafReport1[ ,rsColNames]
ssColumns <- mafReport1[ ,ssColNames]
lprsColumns <- mafReport1[ ,lprsColNames]
lpssColumns <- mafReport1[ ,lpssColNames]

rsRowsNumNAs <- rowSums(is.na(rsColumns))
ssRowsNumNAs <- rowSums(is.na(ssColumns))


rowsToKeep <- logical(length=nrow(mafReport1)) # all initialized to FALSE

for (rowNum in 1:nrow(mafReport1))
{
  # The following conditions in this for-loop fulfill Jun-Jun's requirement that for each row we:
  # - make sure the RS (1-10) and SS (1-10) groups each have at least 50% not NA
  #     - if either has more than 50% NA, substitute the LP-RS or LP-SS samples for their equivalent 
  #       RS or SS samples
  # - check all 20 samples from the RS and SS sample sets to see if there are two or three genotypes
  #     - if 2 genotypes, the RS group should have one, and the SS group the other, with no overlap,
  #       or an overlap with just one sample (10%)
  #     - if 3 genotypes, the RS group should have two, and the SS group the other, with at most one
  #       sample overlapping (10%)
  # - keep only those rows that fulfill the above requirements

  numRsNAs <- rsRowsNumNAs[rowNum]
  numSsNAs <- ssRowsNumNAs[rowNum]

  if (numRsNAs > 5)
  {
    if ( is.na(rsColumns[rowNum, "RS1.tab"]) && !is.na(lprsColumns[rowNum, "LP-RS1.tab"]) )
    {
      rsColumns[rowNum, "RS1.tab"] <- lprsColumns[rowNum, "LP-RS1.tab"]
      numRsNAs <- numRsNAs - 1
    }
    if ( is.na(rsColumns[rowNum, "RS8.tab"]) && !is.na(lprsColumns[rowNum, "LP-RS8.tab"]) )
    {
      rsColumns[rowNum, "RS8.tab"] <- lprsColumns[rowNum, "LP-RS8.tab"]
      numRsNAs <- numRsNAs - 1
    }
    if ( is.na(rsColumns[rowNum, "RS9.tab"]) && !is.na(lprsColumns[rowNum, "LP-RS9.tab"]) )
    {
      rsColumns[rowNum, "RS9.tab"] <- lprsColumns[rowNum, "LP-RS9.tab"]
      numRsNAs <- numRsNAs - 1
    }
  }

  if (numSsNAs > 5)
  {
    if ( is.na(ssColumns[rowNum, "SS4.tab"]) && !is.na(lpssColumns[rowNum, "LP-SS4.tab"]) )
    {
      ssColumns[rowNum, "SS4.tab"] <- lpssColumns[rowNum, "LP-SS4.tab"]
      numSsNAs <- numSsNAs - 1
    }
    if ( is.na(ssColumns[rowNum, "SS5.tab"]) && !is.na(lpssColumns[rowNum, "LP-SS5.tab"]) )
    {
      ssColumns[rowNum, "SS5.tab"] <- lpssColumns[rowNum, "LP-SS5.tab"]
      numSsNAs <- numSsNAs - 1
    }
    if ( is.na(ssColumns[rowNum, "SS6.tab"]) && !is.na(lpssColumns[rowNum, "LP-SS6.tab"]) )
    {
      ssColumns[rowNum, "SS6.tab"] <- lpssColumns[rowNum, "LP-SS6.tab"]
      numSsNAs <- numSsNAs - 1
    }
  }

  if (numRsNAs <= 5 && numSsNAs <= 5) 
  {
    rsRowData <- as.matrix(rsColumns[rowNum, ]) 
    tabulatedRsRowData <- table(rsRowData, useNA = "no")
    tabulatedRsRowDataNames <- names(tabulatedRsRowData)

    ssRowData <- as.matrix(ssColumns[rowNum, ]) 
    tabulatedSsRowData <- table(ssRowData, useNA = "no")
    tabulatedSsRowDataNames <- names(tabulatedSsRowData)
    
    numRsAlleles <- length(tabulatedRsRowDataNames)
    numSsAlleles <- length(tabulatedSsRowDataNames)
    
    if (!accept10percentOverlap)
    {
      # keep row if RS and SS samples each have only one allele, and not the same allele    
      if (numRsAlleles == 1 && numSsAlleles == 1 && tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[1]])
      {
        rowsToKeep[rowNum] <- TRUE
      }
      # keep row if RS has 2 alleles and SS samples have one allele, and not the same as RS alleles
      else if ( numRsAlleles == 2 && numSsAlleles == 1
              && tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[1]] 
              && tabulatedRsRowDataNames[[2]] != tabulatedSsRowDataNames[[1]] )
      {
        rowsToKeep[rowNum] <- TRUE
      }
    }
    else
    {
      # keep row if RS and SS samples each have only one allele, and not the same allele    
      if (numRsAlleles == 1 && numSsAlleles == 1 && tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[1]])
      {
        rowsToKeep[rowNum] <- TRUE
      }
      else if (numRsAlleles == 1 && numSsAlleles == 2)
      {
        # keep row if RS has 1 allele and SS samples have two alleles, and the first SS allele is not the 
        # same as the RS allele, and the second SS allele matches the RS allele but only occurs once in 
        # the SS data set
        if ( tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[1]] 
                && tabulatedRsRowDataNames[[1]] == tabulatedSsRowDataNames[[2]] 
                && tabulatedSsRowData[[ tabulatedSsRowDataNames[[2]] ]] == 1 )
        {
          rowsToKeep[rowNum] <- TRUE
        }
        # keep row if RS has 1 allele and SS samples have two alleles, and the second SS allele is not the 
        # same as the RS allele, and the first SS allele matches the RS allele but only occurs once in 
        # the SS data set
        else if ( tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[2]] 
                && tabulatedRsRowDataNames[[1]] == tabulatedSsRowDataNames[[1]] 
                && tabulatedSsRowData[[ tabulatedSsRowDataNames[[1]] ]] == 1 )
        {
          rowsToKeep[rowNum] <- TRUE
        }
      }
      else if (numRsAlleles == 2 && numSsAlleles == 1)
      {
        # keep row if RS has 2 alleles and SS samples have one allele, and not the same as RS alleles
        if ( tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[1]] 
                && tabulatedRsRowDataNames[[2]] != tabulatedSsRowDataNames[[1]] )
        {
          rowsToKeep[rowNum] <- TRUE
        }
        # keep row if RS has 2 alleles and SS samples have one allele that matches the second RS allele, 
        # but the second RS allele only occurs once in the RS data set
        else if ( tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[1]] 
                && tabulatedRsRowDataNames[[2]] == tabulatedSsRowDataNames[[1]] 
                && tabulatedRsRowData[[ tabulatedRsRowDataNames[[2]] ]] == 1 )
        {
          rowsToKeep[rowNum] <- TRUE
        }
        # keep row if RS has 2 alleles and SS samples have one allele that matches the first RS allele, 
        # but the first RS allele only occurs once in the RS data set
        else if ( tabulatedRsRowDataNames[[1]] == tabulatedSsRowDataNames[[1]] 
                && tabulatedRsRowDataNames[[2]] != tabulatedSsRowDataNames[[1]] 
                && tabulatedRsRowData[[ tabulatedRsRowDataNames[[1]] ]] == 1 )
        {
          rowsToKeep[rowNum] <- TRUE
        }
      }
      else if (numRsAlleles == 2 && numSsAlleles == 2)
      {
        # keep row if RS has 2 alleles and SS samples have two alleles, and the first SS allele does not match
        # the RS alleles, and the second SS allele only occurs once and matches either of the RS alleles
        if ( tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[1]] 
          && tabulatedRsRowDataNames[[2]] != tabulatedSsRowDataNames[[1]] 
          && (tabulatedRsRowDataNames[[1]] == tabulatedSsRowDataNames[[2]] || tabulatedRsRowDataNames[[2]] == tabulatedSsRowDataNames[[2]])
          && tabulatedSsRowData[[ tabulatedSsRowDataNames[[2]] ]] == 1 )
        {
          rowsToKeep[rowNum] <- TRUE
        }
        # keep row if RS has 2 alleles and SS samples have two alleles, and the second SS allele does not match
        # the RS alleles, and the first SS allele only occurs once and matches either of the RS alleles
        else if ( tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[2]] 
          && tabulatedRsRowDataNames[[2]] != tabulatedSsRowDataNames[[2]] 
          && (tabulatedRsRowDataNames[[1]] == tabulatedSsRowDataNames[[1]] || tabulatedRsRowDataNames[[2]] == tabulatedSsRowDataNames[[1]])
          && tabulatedSsRowData[[ tabulatedSsRowDataNames[[1]] ]] == 1 )
        {
          rowsToKeep[rowNum] <- TRUE
        }
      }
      else if (numRsAlleles == 3 && numSsAlleles == 1)
      {
        # keep row if RS has 3 alleles and SS samples have one allele, and the first two RS alleles don't match
        # the SS allele, and the third RS allele only occurs once and matches the SS allele
        if ( tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[1]] 
          && tabulatedRsRowDataNames[[2]] != tabulatedSsRowDataNames[[1]] 
          && tabulatedRsRowDataNames[[3]] == tabulatedSsRowDataNames[[1]] 
          && tabulatedRsRowData[[ tabulatedRsRowDataNames[[3]] ]] == 1 )
        {
          rowsToKeep[rowNum] <- TRUE
        }
        # keep row if RS has 3 alleles and SS samples have one allele, and the first and third RS alleles don't 
        # match the SS allele, and the second RS allele only occurs once and matches the SS allele
        else if ( tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[1]] 
          && tabulatedRsRowDataNames[[3]] != tabulatedSsRowDataNames[[1]] 
          && tabulatedRsRowDataNames[[2]] == tabulatedSsRowDataNames[[1]] 
          && tabulatedRsRowData[[ tabulatedRsRowDataNames[[2]] ]] == 1 )
        {
          rowsToKeep[rowNum] <- TRUE
        }
        # keep row if RS has 3 alleles and SS samples have one allele, and the second and third RS alleles don't 
        # match the SS allele, and the first RS allele only occurs once and matches the SS allele
        else if ( tabulatedRsRowDataNames[[2]] != tabulatedSsRowDataNames[[1]] 
          && tabulatedRsRowDataNames[[3]] != tabulatedSsRowDataNames[[1]] 
          && tabulatedRsRowDataNames[[1]] == tabulatedSsRowDataNames[[1]] 
          && tabulatedRsRowData[[ tabulatedRsRowDataNames[[1]] ]] == 1 )
        {
          rowsToKeep[rowNum] <- TRUE
        }
      }
    } # end of 10-percent overlap else-statement
  } 
  
} # end for-loop


# filter out all rows except the ones which we've marked to keep
mafReport1 <- mafReport1[which(rowsToKeep), ]

write(paste0("================================================"), stdout())
write(paste0("Writing final output csv's. "), stdout())

write.csv(mafReport1, paste(path.expand(path), outputMAFcutoffReportName,sep="/"), row.names=FALSE)

write(paste0("FINISHED."), stdout())

