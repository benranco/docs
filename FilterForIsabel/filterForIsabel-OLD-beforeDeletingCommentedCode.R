options(stringsAsFactors = FALSE, warn = 1)
options(scipen = 999) # disable scientific notation for large numbers

args <- commandArgs(trailingOnly = TRUE)


#message("running report generation part 2")
#if(!require(seqinr))
#{
#  install.packages('seqinr', repos='http://cran.us.r-project.org')
#}
#
#if(!require(VariantAnnotation))
#{
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("VariantAnnotation")
#}
#
#library(VariantAnnotation)
#library(seqinr)


write("Running filterForIsabel.R.", stdout())

# #######################################################################

#MAF_CUTOFF <- 0.1

#basePath <- "/run/media/benrancourt/27a1f1ff-6590-4f38-aa2c-13f6f1e0e771/ben/SNPpipeline-isabel"
basePath <- args[1]

#pathToFastaRef <- paste(path.expand(basePath),"reference",sep="/")
pathToFastaRef <- args[2]

#fastaRefFilename <- "formatted_output.fasta"
fastaRefFilename <- args[3]

#pathToReports <- paste(path.expand(basePath),"reports/reports-backup/done-onlyMissingData",sep="/")
pathToReports <- args[4]

#filledReportName <- "filled_report.csv"
filledReportName <- args[5]

#percentageSnpReportName <- "percentage_snps.csv"

#outputReportName <- "filtered_report.csv"
outputReportName <- args[6]

#outputPercentageSnpReportName <- "MAF_cutoff_10_percentage_snps_report.csv"


# #######################################################################

filledReport <- read.csv(paste(path.expand(pathToReports),filledReportName,sep="/"),header=TRUE)
# check if first col is not CHROM, and remove it if so (it's just a left-over old index column if it exists)
if (colnames(filledReport)[1] != "CHROM")
{
  filledReport <- filledReport[,-1]
}

#percentageSnpReport <- read.csv(paste(path.expand(pathToReports),percentageSnpReportName,sep="/"),header=TRUE)
## check if first col is not CHROM, and remove it if so (it's just a left-over old index column if it exists)
#if (colnames(percentageSnpReport)[1] != "CHROM")
#{
#  percentageSnpReport <- percentageSnpReport[,-1]
#}

# remove all rows with at least one NA
message("filtering out rows that contain at least one NA")
rowsWithNoNAs <- complete.cases(filledReport)
filledReport <- filledReport[rowsWithNoNAs, ]
#percentageSnpReport <- percentageSnpReport[rowsWithNoNAs, ]


# #######################################################################
message("filtering out rows whose samples are all the same")
  
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

  # remove row if there's only one type of non-NA data
  if( length(tabulatedRowDataNames) == 1 )
  {
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

# #######################################################################


# #######################################################################







 
#rowsThatPassMAFCutoff <- logical(nrow(percentageSnpReport))
#
#for(i in 1:nrow(percentageSnpReport))
#{
#    if (percentageSnpReport[i,"MAF"] >= MAF_CUTOFF)
#    {
#        rowsThatPassMAFCutoff[i] <- TRUE
#    }
#}
#
#filledReport <- filledReport[rowsThatPassMAFCutoff, ]
#percentageSnpReport <- percentageSnpReport[rowsThatPassMAFCutoff, ]

write.csv(filledReport, paste(path.expand(pathToReports), outputReportName,sep="/"), row.names=FALSE)
#write.csv(percentageSnpReport, paste(path.expand(pathToReports), outputPercentageSnpReportName,sep="/"), row.names=FALSE)

# #######################################################################
#message("generating site mutation percentage data")
#
#fastaRef <- read.fasta(file = paste(path.expand(pathToFastaRef), fastaRefFilename,sep="/"), as.string = TRUE)
#mutationReport <- data.frame()
#
#for(sector in 1:length(names(fastaRef)))
#{
#  sectorName <- attributes(fastaRef[[sector]])$name
#  numRowsForSectorName <- length(grep(sectorName, filledReport[, "CHROM"], fixed=TRUE))
#  numCharsInSequence <- nchar(fastaRef[[sector]][1])
#
#  mutationReport[sectorName, "role"] <- attributes(fastaRef[[sector]])$Annot  
#  mutationReport[sectorName, "snp"] <- numRowsForSectorName
#  mutationReport[sectorName, "length"] <- numCharsInSequence
#  mutationReport[sectorName, "percentage SNP"] <- numRowsForSectorName / numCharsInSequence
#}
#
##mutationReport <- mutationReport[order(-mutationReport$`percentage SNP`), ]
#mutationReport$`percentage SNP` <- mutationReport$`percentage SNP` * 100
#
#write.csv(mutationReport, paste(path.expand(pathToReports), "MAF_cutoff_10_mutation_percentage.csv", sep = "/"), row.names=TRUE)

