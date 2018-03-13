options(stringsAsFactors = FALSE, warn = 1)
options(scipen = 999) # disable scientific notation for large numbers

args <- commandArgs(trailingOnly = TRUE)


write("Running filterForIsabel_part3.R.", stdout())

# #######################################################################

#basePath <- "/run/media/benrancourt/27a1f1ff-6590-4f38-aa2c-13f6f1e0e771/ben/SNPpipeline-isabel"
basePath <- args[1]

#pathToReports <- paste(path.expand(basePath),"reports/reports-backup/done-exceptMissingData",sep="/")
pathToReports <- paste(path.expand(basePath),args[2],sep="/")

#depthReportName <- "filtered_report.csv_depth.csv"
depthReportName <- args[3]

#depthDetailedReportName <- "filtered_report.csv_depth_detailed.csv"
depthDetailedReportName <- args[4]

#inputReportName <- "filtered_report.csv"
inputReportName <- args[5]

#outputReportName <- "filtered_by_depth.csv"
outputReportName <- args[6]

#outputStrReplacedReportName <- "filtered_by_depth_recoded.csv"
outputStrReplacedReportName <- args[7]


# #######################################################################
# Purpose: Only keep the rows in filtered_report.csv whose coverage for each sample (column) is more than ten.


inputReport <- read.csv(paste(path.expand(pathToReports),inputReportName,sep="/"),header=TRUE)
# check if first col is not CHROM, and remove it if so (it's just a left-over old index column if it exists)
if (colnames(inputReport)[1] != "CHROM")
{
  inputReport <- inputReport[,-1]
}

depthReport <- read.csv(paste(path.expand(pathToReports),depthReportName,sep="/"),header=TRUE)
# check if first col is not CHROM, and remove it if so (it's just a left-over old index column if it exists)
if (colnames(depthReport)[1] != "CHROM")
{
  depthReport <- depthReport[,-1]
}

depthDetailedReport <- read.csv(paste(path.expand(pathToReports),depthDetailedReportName,sep="/"),header=TRUE)
# check if first col is not CHROM, and remove it if so (it's just a left-over old index column if it exists)
if (colnames(depthDetailedReport)[1] != "CHROM")
{
  depthDetailedReport <- depthDetailedReport[,-1]
}


if(nrow(depthReport) != nrow(inputReport))
{
  write(paste0("Aborting because the number of rows in ",depthReport," don't match the number of rows in ",inputReport,"."), stdout())
} else {
  startCol <- 4
  if("COMBINED" %in% colnames(depthReport))
  {
    startCol <- 5
  }

  rowsToRemove <- numeric(nrow(depthReport))
  rowsToRemoveIndex <- 1
  reportIndex <- 1
  
  while(reportIndex <= nrow(depthReport))
  {
    rowData <- as.matrix(depthReport[reportIndex,startCol:ncol(depthReport)]) 

    # remove row if at least one depth value is less than 10
    if( min(rowData) < 10 )
    {
      rowsToRemove[rowsToRemoveIndex] <- reportIndex
      rowsToRemoveIndex <- rowsToRemoveIndex + 1        
    }
    reportIndex <- reportIndex + 1
  }

  rowsToRemove <- rowsToRemove[!(rowsToRemove %in% c(0))] # remove all elements containing 0 from rowsToRemove
  # remove the selected rows from inputReport
  if(length(rowsToRemove) > 0) 
  {
    inputReport <- inputReport[-rowsToRemove[1:length(rowsToRemove)], ]
    depthReport <- depthReport[-rowsToRemove[1:length(rowsToRemove)], ]
    depthDetailedReport <- depthDetailedReport[-rowsToRemove[1:length(rowsToRemove)], ]
  }

  strReplacedReport <- inputReport
  strReplacedReport[strReplacedReport=="A/A"]<-"A"
  strReplacedReport[strReplacedReport=="C/C"]<-"C"
  strReplacedReport[strReplacedReport=="G/G"]<-"G"
  strReplacedReport[strReplacedReport=="T/T"]<-"T"

  strReplacedReport[strReplacedReport=="A/T"]<-"W"
  strReplacedReport[strReplacedReport=="T/A"]<-"W"

  strReplacedReport[strReplacedReport=="C/G"]<-"S"
  strReplacedReport[strReplacedReport=="G/C"]<-"S"

  strReplacedReport[strReplacedReport=="A/C"]<-"M"
  strReplacedReport[strReplacedReport=="C/A"]<-"M"

  strReplacedReport[strReplacedReport=="G/T"]<-"K"
  strReplacedReport[strReplacedReport=="T/G"]<-"K"

  strReplacedReport[strReplacedReport=="A/G"]<-"R"
  strReplacedReport[strReplacedReport=="G/A"]<-"R"

  strReplacedReport[strReplacedReport=="C/T"]<-"Y"
  strReplacedReport[strReplacedReport=="T/C"]<-"Y"


  write.csv(inputReport, paste(path.expand(pathToReports), outputReportName,sep="/"), row.names=FALSE)
  write.csv(strReplacedReport, paste(path.expand(pathToReports), outputStrReplacedReportName,sep="/"), row.names=FALSE)
  write.csv(depthReport, paste(path.expand(pathToReports), paste0(depthReportName,"_filtered.csv"),sep="/"), row.names=FALSE)
  write.csv(depthDetailedReport, paste(path.expand(pathToReports), paste0(depthDetailedReportName,"_filtered.csv"), sep="/"), row.names=FALSE)


}

