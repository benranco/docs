options(stringsAsFactors = FALSE, warn = 1)
options(scipen = 999) # disable scientific notation for large numbers

write("Running getDepthStats.R.", stdout())


# ########################################################
# Input Parameters:

args <- commandArgs(trailingOnly = TRUE)

path <- "/home/benrancourt/Desktop/junjun/fluidigm-SNPpipelineFupan-dev"

reportsSubDir <- "reports"
inputMAFcutoffReportName <- "MAF_cutoff_report.csv"

vcfFilesSubDir <- "outputTemp/single"
cutoffFileNamePostFix <- "_cutoff"

bamFilesSubDir <- "dataTemp/single"
bamFileNamePostFix <- "_sorted.bam"

depthFilesSubDir <- "dataTemp/depthfiles"
depthFileNamePostFix <- "_depth.txt"

outputMAFcutoffReportNameDepth <- "MAF_cutoff_report_depth.csv"
outputMAFcutoffReportNameDepthDetailed <- "MAF_cutoff_report_depth_detailed.csv"

samtoolsPathAndExecutable <- "/home/benrancourt/Desktop/junjun/fluidigm-SNPpipelineFupan-dev/tools/samtools-1.3.1/samtools"

# ########################################################
# Execution code:

# Basic Process:
# - open MAF_cutoff_report
# - create two copies of it
# - for each sample (col):
#     - prepare vcf file name
#     - call samtools to create a depth file from the .bam file
#     - iterate down the col, and for each cell:
#         - find the row in the vcf file with the CHROM id and POS in the VCF file
#         - if the row exists:
#             - save the DP data (search for ";DP=165;") in the cell in copy 1
#             - save the DP;DPB;RO;AO (search for ";DP=165;", ";DPB=165;", ";AO=127;", ";RO=38;") in copy 2
#         - if the row doesn't exist:
#             - find the row in the depth file with the CHROM id and POS, and save its depth value in both 
#               copy 1 and copy 2 tables, using NA for the other three values in copy 2.
# - write output to .csv files

mafReport1 <- read.csv(paste(path.expand(path),reportsSubDir,inputMAFcutoffReportName,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.


# check if first col is not CHROM, and remove it if so (it's just a left-over old index column if it exists)
if (colnames(mafReport1)[1] != "CHROM")
{
  mafReport1 <- mafReport1[,-1]
}

mafReport2 <- mafReport1

startingCol <- 4
# check if COMBINED col exists, and adjust starting col num 
if (colnames(mafReport1)[4] == "COMBINED")
{
  startingCol <- 5
}

for (colNum in startingCol:10)
#for (colNum in startingCol:ncol(mafReport1))
{
  write(paste0("Processing column ",colNum), stdout())
  sampleName <- colnames(mafReport1)[colNum]
  if ( substr(sampleName, nchar(sampleName)-3, nchar(sampleName)) == ".tab" )
  {
    sampleName <- paste0( substr(sampleName, 1, nchar(sampleName)-4) )
  }

  vcfFileName <- paste0(sampleName, cutoffFileNamePostFix)
  vcfFilePath <- paste(path.expand(path), vcfFilesSubDir, vcfFileName,sep="/")

  # use samtools to generate a depth file for the current .bam file (for this column):
  depthOutputFileName <- paste0(sampleName, depthFileNamePostFix)
  depthOutputFilePath <- paste(path.expand(path), depthFilesSubDir, depthOutputFileName, sep="/")

  bamFileName <- paste0(sampleName, bamFileNamePostFix)
  bamFilePath <- paste(path.expand(path), bamFilesSubDir, bamFileName, sep="/")

  system2( samtoolsPathAndExecutable, c("depth", "-aa", bamFilePath), stdout=depthOutputFilePath, stderr=NULL )

  options(warn=-1) # globally suppress warnings

  for (rowNum in 1:100)  
  #for (rowNum in 1:nrow(mafReport1))
  {
    sequenceName <- mafReport1[rowNum,1]
    position <- mafReport1[rowNum,2]
    dp <- NA
    dpb <- NA
    ao <- NA
    ro <- NA
    # get data from vcf file
    #vcfLine <- system( paste0("grep \"^",sequenceName,"\\s",position,"\\s\" ", vcfFilePath), intern=TRUE)
    vcfLine <- system2( "grep", c(paste0("\"^",sequenceName,"\\s",position,"\\s\""), vcfFilePath), stdout=TRUE, stderr=NULL )
    if (length(vcfLine) == 1) # if there is one line in vcfLine
    { 
      # extract the DP, DPB, AO, and RO data from the line

      # remove everything preceding the value we want to extract
      dp <- sub(".*;DP=", "", vcfLine)
      dpb <- sub(".*;DPB=", "", vcfLine)
      ao <- sub(".*;AO=", "", vcfLine)
      ro <- sub(".*;RO=", "", vcfLine)

      # remove everything following the value we want to extract
      dp <- sub(";.*", "", dp)
      dpb <- sub(";.*", "", dpb)
      ao <- sub(";.*", "", ao)
      ro <- sub(";.*", "", ro)      
    } 
    else if (length(vcfLine) == 0) # if the data is not found in the vcf file
    {
      # get data from the depth file
      depthLine <- system2( "grep", c(paste0("\"^",sequenceName,"\\s",position,"\\s\""), depthOutputFilePath), stdout=TRUE, stderr=NULL )
      if (length(depthLine) == 1)
      {
        dp <- strsplit(depthLine, "\\s")[[1]][3]
      }
    }
    
    mafReport1[rowNum, colNum] <- dp  # basic data table
    mafReport2[rowNum, colNum] <- paste(dp,dpb,ao,ro,sep=";")  # detailed data table
  }
  
  options(warn=0) # reenable warnings
}


write(paste0("================================================"), stdout())
write(paste0("Writing final output csv's. "), stdout())

write.csv(mafReport1, paste(path.expand(path), reportsSubDir, outputMAFcutoffReportNameDepth,sep="/"), row.names=FALSE)
write.csv(mafReport2, paste(path.expand(path), reportsSubDir, outputMAFcutoffReportNameDepthDetailed,sep="/"), row.names=FALSE)

write(paste0("FINISHED."), stdout())

