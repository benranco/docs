options(stringsAsFactors = FALSE, warn = 1)
options(scipen = 999) # disable scientific notation for large numbers

args <- commandArgs(trailingOnly = TRUE)


write("Running filterForIsabel_part1.R.", stdout())

# #######################################################################

#basePath <- "/run/media/benrancourt/27a1f1ff-6590-4f38-aa2c-13f6f1e0e771/ben/SNPpipeline-isabel"
basePath <- args[1]

#pathToReports <- paste(path.expand(basePath),"reports/reports-backup/done-onlyMissingData",sep="/")
pathToReports <- paste(path.expand(basePath),args[2],sep="/")

#filledReportName <- "filled_report.csv"
filledReportName <- args[3]

#outputReportName <- "filtered_report.csv"
outputReportName <- args[4]

minNumAlternateValues <- 10

# #######################################################################

filledReport <- read.csv(paste(path.expand(pathToReports),filledReportName,sep="/"),header=TRUE)
# check if first col is not CHROM, and remove it if so (it's just a left-over old index column if it exists)
if (colnames(filledReport)[1] != "CHROM")
{
  filledReport <- filledReport[,-1]
}

# remove all rows with at least one NA
message("filtering out rows that contain at least one NA")
rowsWithNoNAs <- complete.cases(filledReport)
filledReport <- filledReport[rowsWithNoNAs, ]
#percentageSnpReport <- percentageSnpReport[rowsWithNoNAs, ]


# #######################################################################
message("filtering out rows who don't have enough differing samples")
  
#report <- as.data.frame(report)

startCol <- 4
if("COMBINED" %in% colnames(filledReport))
{
  startCol <- 5
}

rowsToRemove <- numeric(nrow(filledReport))
rowsToRemoveIndex <- 1
reportIndex <- 1

while(reportIndex <= nrow(filledReport))
{
  rowData <- as.matrix(filledReport[reportIndex,startCol:ncol(filledReport)]) 
  tabulatedRowData <- table(rowData, useNA = "no")
  tabulatedRowDataNames <- names(tabulatedRowData)

#  # remove row if there's only one type of non-NA data
#  if( length(tabulatedRowDataNames) == 1 )
#  {
#    rowsToRemove[rowsToRemoveIndex] <- reportIndex
#    rowsToRemoveIndex <- rowsToRemoveIndex + 1        
#  }

  # remove row if the number of samples that don't match the most common value is less than minNumAlternateValues
  if( sum(tabulatedRowData[-which.max(tabulatedRowData)]) < minNumAlternateValues ) {
    rowsToRemove[rowsToRemoveIndex] <- reportIndex
    rowsToRemoveIndex <- rowsToRemoveIndex + 1            
  }

  reportIndex <- reportIndex + 1
}

rowsToRemove <- rowsToRemove[!(rowsToRemove %in% c(0))] # remove all elements containing 0 from rowsToRemove
# remove the selected rows from filledReport
if(length(rowsToRemove) > 0) 
{
  filledReport <- filledReport[-rowsToRemove[1:length(rowsToRemove)], ]
}


write.csv(filledReport, paste(path.expand(pathToReports), outputReportName,sep="/"), row.names=FALSE)

