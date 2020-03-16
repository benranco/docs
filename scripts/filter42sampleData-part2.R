options(stringsAsFactors = FALSE, warn = 1)
options(scipen = 999) # disable scientific notation for large numbers

write("Running filter42sampleData-part2.R.", stdout())

# This script is based on (modified) the original filter46sampleData.R script from Aug 2017, which was used to filter the 46 sample data set with Exom14706 ref.


# ########################################################
# Input Parameters:

args <- commandArgs(trailingOnly = TRUE)

path <- "/home/benrancourt/Desktop/junjun/reports-42samples-corrected-allDataIncluding1000rowsfromP35"

inputMAFcutoffReportName <- "MAF_cutoff_report_filtered.csv"

#outputMAFcutoffReportName <- "MAF_cutoff_report_filtered_part2_no_overlap.csv"
outputMAFcutoffReportName <- "MAF_cutoff_report_filtered_part2_overlapOf1-b.csv"

accept10percentOverlap <- TRUE

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

rsColNames <- c("HI.4542.006.Index_9.BC11.tab", "HI.4542.007.Index_11.BC111.tab", "HI.4542.007.Index_8.BC16.tab", "HI.4542.007.Index_10.BC29.tab", "HI.4542.007.Index_20.BC615.tab", "HI.4542.007.Index_22.BC616.tab", "S001AB3_ATCACG_L002.tab", "S001AB4_CGATGT_L002.tab", "S001AB6_TTAGGC_L002.tab", "S001AA3_GATCAG_L002.tab", "S001AA8_ACTTGA_L002.tab", "S001AAF_CAGATC_L002.tab", "S0008F3_GGCTAC_L002.tab", "S0008F5_GATCAG_L002.tab", "S0008F8_TAGCTT_L002.tab")

ssColNames <- c("HI.4542.002.Index_9.HC1.tab", "HI.4542.003.Index_21.HC12.tab", "HI.4542.003.Index_23.HC14.tab", "HI.4542.003.Index_27.HC15.tab", "HI.4542.004.Index_2.HC16.tab", "HI.4542.004.Index_13.HC17.tab", "HI.4542.004.Index_6.HC19.tab", "HI.4542.004.Index_15.HC20.tab", "HI.4542.004.Index_7.HC21.tab", "HI.4542.005.Index_18.HC22.tab", "HI.4542.005.Index_14.HC23.tab", "HI.4542.005.Index_16.HC24.tab", "HI.4542.005.Index_4.HC26.tab", "HI.4542.005.Index_5.HC27.tab", "HI.4542.006.Index_12.HC28.tab", "HI.4542.006.Index_19.HC29.tab", "HI.4542.002.Index_8.HC3.tab", "HI.4542.002.Index_10.HC4.tab", "HI.4542.002.Index_11.HC5.tab", "HI.4542.002.Index_20.HC7.tab", "HI.4542.003.Index_22.HC8.tab", "HI.4542.003.Index_25.HC9.tab")

rsColumns <- mafReport1[ ,rsColNames]
ssColumns <- mafReport1[ ,ssColNames]


rowsToKeep <- logical(length=nrow(mafReport1)) # all initialized to FALSE

for (rowNum in 1:nrow(mafReport1))
{
  # The following conditions in this for-loop fulfill Jun-Jun's requirement that for each row we:
  # - check all samples from the RS and SS sample sets to see if there are two or three genotypes
  #     - if 2 genotypes, the RS group should have one, and the SS group the other, with no overlap,
  #       or an overlap with just one sample (10%)
  #     - if 3 genotypes, the RS group should have two, and the SS group the other, with at most one
  #       sample overlapping (10%)
  # - keep only those rows that fulfill the above requirements

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
      # keep row if RS has 1 allele and SS samples have 2 alleles, and not the same as RS alleles
      else if ( numRsAlleles == 1 && numSsAlleles == 2
              && tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[1]] 
              && tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[2]] )
      {
        rowsToKeep[rowNum] <- TRUE
      }
      # keep row if RS has 3 alleles (e.g. G/A is same as A/G) and SS samples have one allele, and not the same as RS alleles
      else if ( numRsAlleles == 2 && numSsAlleles == 1
              && tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[1]] 
              && tabulatedRsRowDataNames[[2]] != tabulatedSsRowDataNames[[1]] 
              && tabulatedRsRowDataNames[[3]] != tabulatedSsRowDataNames[[1]] )
      {
        rowsToKeep[rowNum] <- TRUE
      }
      # keep row if RS has 1 allele and SS samples have 3 alleles (e.g. G/A is same as A/G), and not the same as RS alleles
      else if ( numRsAlleles == 1 && numSsAlleles == 2
              && tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[1]] 
              && tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[2]] 
              && tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[3]] )
      {
        rowsToKeep[rowNum] <- TRUE
      }
      # keep row if RS has 2 alleles and SS samples have 2 alleles (e.g. G/A is same as A/G), and not the same as RS alleles
      else if ( numRsAlleles == 1 && numSsAlleles == 2
              && tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[1]] 
              && tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[2]] 
              && tabulatedRsRowDataNames[[2]] != tabulatedSsRowDataNames[[1]] 
              && tabulatedRsRowDataNames[[2]] != tabulatedSsRowDataNames[[2]] )
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
        # keep row if RS has 1 allele and SS samples have 2 alleles, and not the same as RS alleles
        if ( tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[1]] 
                && tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[2]] )
        {
          rowsToKeep[rowNum] <- TRUE
        }
        # keep row if RS has 1 allele and SS samples have two alleles, and the first SS allele is not the 
        # same as the RS allele, and the second SS allele matches the RS allele but only occurs once in 
        # the SS data set
        else if ( tabulatedRsRowDataNames[[1]] != tabulatedSsRowDataNames[[1]] 
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

  
} # end for-loop


# filter out all rows except the ones which we've marked to keep
mafReport1 <- mafReport1[which(rowsToKeep), ]

write(paste0("================================================"), stdout())
write(paste0("Writing final output csv's. "), stdout())

write.csv(mafReport1, paste(path.expand(path), outputMAFcutoffReportName,sep="/"), row.names=FALSE)

write(paste0("FINISHED."), stdout())

