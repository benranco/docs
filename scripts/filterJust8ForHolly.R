options(stringsAsFactors = FALSE, warn = 1)
options(scipen = 999) # disable scientific notation for large numbers

write("Running filterForHolly.R.", stdout())

# ########################################################
# Initial Input Parameters:

args <- commandArgs(trailingOnly = TRUE)

path <- "/home/benrancourt/Downloads/Holly"

inputReportName <- "Fold-change&FDR_pval_tol&susc.csv"

outputReportNameFirstPart <- "Fold-change&FDR_pval_tol&susc-filter"
outputReportNameRatioPart <- "-byRatio"


# ########################################################
# Function filterReport:

filterReport <- function(rowKeeper,rowNum,filterId, minRatio, goodFam1Col, goodFam2Col, badFam1Col, badFam2Col, pval1Col, pval2Col)
{
    thisRow <- report[rowNum, ]

    rowHasNAs <- any(is.na(thisRow))
    # deal with case where there is some NA data in the row
    if (rowHasNAs)
    {
      rowKeeper[rowNum] <- TRUE
    }
    else 
    {
      goodFamAverage <- (thisRow[goodFam1Col] + thisRow[goodFam2Col]) / 2
      badFamAverage <- (thisRow[badFam1Col] + thisRow[badFam2Col]) / 2
      
      # We are only interested in rows whose averages are greater than 1.8
      if (abs(goodFamAverage) >= 1.8 && abs(badFamAverage) >= 1.8)
      {
        # keep row if the bad family average is negative or 0 and the good family average is positive
        if (badFamAverage <= 0 && goodFamAverage > 0)
        {
          rowKeeper[rowNum] <- TRUE
        }
        # keep row if the absolute value of the ratio of good to bad passes the test
        else if (goodFamAverage != 0)
        {
          rowKeeper[rowNum] <- abs(goodFamAverage / badFamAverage) >= minRatio 
        }
      }
    }
    
    # return rowKeeper
    rowKeeper
} # end function filterReport


# ########################################################
# Function filterReportUsingOr:

filterReportUsingOr <- function(filterId, minRatio, goodFam1Col, goodFam2Col, badFam1Col, badFam2Col, pval1Col, pval2Col)
{
  write(paste0("Running filter ",filterId," with ratio ",minRatio," using OR."), stdout())

  keepRows <- logical(nrow(report))

  for (i in 1:nrow(report))
  {
    thisRow <- report[i, ]

    # we are only interested in rows that meet this pval condition
    if (thisRow[pval1Col] < 0.05 || thisRow[pval2Col] < 0.05)
    {
      keepRows <- filterReport(keepRows,i,filterId, minRatio, goodFam1Col, goodFam2Col, badFam1Col, badFam2Col, pval1Col, pval2Col)
    } # end pval if-condition

  } # end for-loop

  
  filteredReport <- report[which(keepRows), ] 
  write.csv(filteredReport, paste(path.expand(path), paste0(outputReportNameFirstPart, filterId, outputReportNameRatioPart, minimumRatio, "-usingOrPvals.csv"),sep="/"), row.names=FALSE)
 
} # end function filterReportUsingOr

# ########################################################
# Function filterReportUsingAnd:

filterReportUsingAnd <- function(filterId, minRatio, goodFam1Col, goodFam2Col, badFam1Col, badFam2Col, pval1Col, pval2Col)
{
  write(paste0("Running filter ",filterId," with ratio ",minRatio," using AND."), stdout())

  keepRows <- logical(nrow(report))

  for (i in 1:nrow(report))
  {
    thisRow <- report[i, ]

    # we are only interested in rows that meet this pval condition
    if (thisRow[pval1Col] < 0.05 && thisRow[pval2Col] < 0.05)
    {
      keepRows <- filterReport(keepRows,i,filterId, minRatio, goodFam1Col, goodFam2Col, badFam1Col, badFam2Col, pval1Col, pval2Col)
    } # end pval if-condition

  } # end for-loop
  

  filteredReport <- report[which(keepRows), ] 
  write.csv(filteredReport, paste(path.expand(path), paste0(outputReportNameFirstPart, filterId, outputReportNameRatioPart, minimumRatio, "-usingAndPvals.csv"),sep="/"), row.names=FALSE)

} # end function filterReportUsingAnd


# ########################################################
# Execution:


report <- read.csv(paste(path.expand(path),inputReportName,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.)



####################################
# Ninth filter

filterNum <- 9

goodFamily1Col <- 12 # col L
goodFamily2Col <- 16 # col P
badFamily1Col <- 4 # col D
badFamily2Col <- 8 # col H

pval1ColNum <- 13 # col M
pval2ColNum <- 17 # col Q


minimumRatio <- 4
filterReportUsingOr(filterNum, minimumRatio, goodFamily1Col, goodFamily2Col, badFamily1Col, badFamily2Col, pval1ColNum, pval2ColNum)





write(paste0("FINISHED."), stdout())

