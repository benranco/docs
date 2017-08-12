options(stringsAsFactors = FALSE, warn = 1)
options(scipen = 999) # disable scientific notation for large numbers

#write("Running getDepthStats.R.", stdout())


# ########################################################
# Input Parameters:

args <- commandArgs(trailingOnly = TRUE)
sampleNumber <- as.numeric(args[1])

#path <- "/home/benrancourt/Desktop/junjun/fluidigm-SNPpipelineFupan-dev"
path <- args[1]

reportsSubDir <- "reports"
inputMAFcutoffReportName <- "MAF_cutoff_report.csv"

outputMAFcutoffReportNameDepth <- "MAF_cutoff_report_depth.csv"
outputMAFcutoffReportNameDepthDetailed <- "MAF_cutoff_report_depth_detailed.csv"


vcfFilesSubDir <- "outputTemp/single"
cutoffFileNamePostFix <- "_cutoff"

bamFilesSubDir <- "dataTemp/single"
bamFileNamePostFix <- "_sorted.bam"

depthFilesSubDir <- "dataTemp/depthfiles"
depthFileNamePostFix <- "_depth.txt"

samtoolsPathAndExecutable <- paste(path.expand(path), "tools/samtools-1.3.1/samtools",sep="/")

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

if(!require(doParallel))
{
  install.packages('doParallel', repos='http://cran.us.r-project.org')
}  
library(doParallel)  


# determine the number of cores we can use for parallel processing
ncore <- as.numeric(system2( "grep", c("-c", "^processor", "/proc/cpuinfo"), stdout=TRUE, stderr=NULL ))
if (ncore > 4)
{
  ncore <- ncore - 4 # leave some processor cores free for other applications
}
write(paste0("Number of cores to use for parallel processing: ",ncore), stdout())

registerDoParallel(cores=ncore)


mafReport <- read.csv(paste(path.expand(path),reportsSubDir,inputMAFcutoffReportName,sep="/"),header=TRUE)

# check if first col is not CHROM, and remove it if so (it's just a left-over old index column if it exists)
if (colnames(mafReport)[1] != "CHROM")
{
  mafReport <- mafReport[,-1]
}

startingCol <- 4
# check if COMBINED col exists, and adjust starting col num 
if (colnames(mafReport)[4] == "COMBINED")
{
  startingCol <- 5
}

mafChromAndPos <- mafReport[ , 1:2] # just need the CHROM and POS columns
sampleNames <- colnames( mafReport[ , startingCol:ncol(mafReport)] )






options(warn=-1) # globally suppress warnings

#parallelResults <- foreach (colNum=startingCol:10) %dopar% {
parallelResults <- foreach (colNum=startingCol:ncol(mafReport)) %dopar% {
  #write(paste0("Processing column ",colNum), stdout())

  columnA <- character(nrow(mafChromAndPos))
  columnB <- columnA

  sampleName <- colnames(mafReport)[colNum]
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



  #for (rowNum in 1:10)  
  for (rowNum in 1:nrow(mafChromAndPos))
  {
    sequenceName <- mafChromAndPos[rowNum,1]
    position <- mafChromAndPos[rowNum,2]
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
    
    columnA[rowNum] <- dp  # basic data table
    columnB[rowNum] <- paste(dp,dpb,ao,ro,sep=";")  # detailed data table
  } # end inner for loop

  result <- data.frame(columnA, columnB)
  colnames(result) <- c(paste0(sampleName, "-DP"), paste0(sampleName, "-DP;DPB;AO;RO"))
  result # the return value of this iteration/parallel process.  

} # end parallel processing loop

options(warn=0) # reenable warnings



write(paste0("Combining individual depth columns into one table."), stdout())

table1 <- mafReport[ , 1:startingCol-1]
table2 <- table1

for (i in 1:length(parallelResults))
{
  table1 <- cbind(table1, parallelResults[[i]][1])
  table2 <- cbind(table2, parallelResults[[i]][2])
}



write(paste0("Combining columns completed."), stdout())

write(paste0("================================================"), stdout())
write(paste0("Writing final output csv's. "), stdout())

write.csv(table1, paste(path.expand(path), reportsSubDir, outputMAFcutoffReportNameDepth,sep="/"), row.names=FALSE)
write.csv(table2, paste(path.expand(path), reportsSubDir, outputMAFcutoffReportNameDepthDetailed,sep="/"), row.names=FALSE)

write(paste0("FINISHED."), stdout())




