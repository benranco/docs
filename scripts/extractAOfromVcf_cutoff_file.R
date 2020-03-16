options(stringsAsFactors = FALSE, warn = 1)
options(scipen = 999) # disable scientific notation for large numbers

write("Running extractAOfromVcf_cutoff_file.R.", stdout())


# ########################################################
# Input Parameters:

# output file will be saved in the same folder as the input file, so make sure you have write permissions
pathToInputFile <- "/work2/AfterJunJunsFiltering"
pathToVcfMasterFolder <- "/work2/vcffiles"


inputFileName <- "snp418.csv"
outputFileName <- "out-snp418.csv"

aoColNames <- c("R_alt","S_alt")
aoCountColNames <- c("R_ao_count","S_ao_count")
roCountColNames <- c("R_ro_count","S_ro_count")
dpColNames <- c("R_dp","S_dp")
relativeCountColNames <- c("Revative-C","Relative-c")
vcfFileNames <- c("Bigfile_RS_cutoff","Bigfile_SC_cutoff")

vcfFilesSubDirs <- c("vcf-pipe01", "vcf-pipe02", "vcf-pipe03", "vcf-pipe04", "vcf-pipe05", "vcf-pipe06", "vcf-pipe07", "vcf-pipe08", "vcf-pipe09", "vcf-pipe10")

# ########################################################
# Execution code:
#
# Basic Process:
# - open the input file
# - for each sample (in parallel):
#     - iterate down the alt col (indicated in input param aoColNames), and for each cell:
#         - iterate through the vcf file subfolders (should be one for each pipeline) and query 
#           the vcf files that correspond to the current sample (one in each subfolder)
#             - find the row in the vcf file with the CHROM id and POS in the VCF file
#             - if the row exists:
#                 - get AO value, etc, save this new data in their own vectors
# - finally, after parallel processing is complete, update the report with the new vectors
# - write output to .csv files
#
# Here's the VCF specification that we use:
# http://samtools.github.io/hts-specs/VCFv4.1.pdf
#
# Here are the definitions from the VCF file of the depth data we've extracted:
#   DP  = "Total read depth at the locus"
#   DPB = "Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype"
#   RO  = "Reference allele observation count, with partial observations recorded fractionally"
#   AO  = "Alternate allele observations, with partial observations recorded fractionally"
###########################################################

if(!require(doParallel))
{
  install.packages('doParallel', repos='http://cran.us.r-project.org')
}  
library(doParallel)  


# determine the number of cores we can use for parallel processing
ncore <- as.numeric(system2( "grep", c("-c", "^processor", "/proc/cpuinfo"), stdout=TRUE, stderr=NULL ))
if (ncore > 1)
{
  ncore <- ncore - 1 # leave some processor cores free for other applications
}
write(paste0("Number of cores to use for parallel processing: ",ncore), stdout())

registerDoParallel(cores=ncore)


report <- read.csv(paste(path.expand(pathToInputFile),inputFileName,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.


chromPosRef <- report[ , 1:3] 



options(warn=-1) # globally suppress warnings for the parallel operations

parallelResults <- foreach (whichAO=1:length(aoColNames)) %dopar% {
  #write(paste0("Processing column ",whichAO), stdout())

  altCol <- report[ , aoColNames[whichAO]]
  altCountCol <- report[ , aoCountColNames[whichAO]]
  relativeCountCol <- report[ , relativeCountColNames[whichAO]]
  refCol <- character(nrow(report))

  # this will be used to record which pipeline folder vcf file that had the data was in
  pipelineCol <- character(nrow(report))  

  for (rowNum in 1:nrow(chromPosRef))
  {
    sequenceName <- chromPosRef[rowNum,1]
    position <- chromPosRef[rowNum,2]
    ref <- chromPosRef[rowNum,3]
    aoCount <- NA
    aoVal <- NA
    roVal <- NA

    for (vcfDir in vcfFilesSubDirs) {
      vcfFilePath <- paste(path.expand(pathToVcfMasterFolder), vcfDir, vcfFileNames[whichAO],sep="/")
      # get data from vcf file
      vcfLine <- system2( "grep", c(paste0("\"^",sequenceName,"\\s",position,"\\s\""), vcfFilePath), stdout=TRUE, stderr=NULL )
      # This commented out version also makes sure the REF matches too:
      #vcfLine <- system2( "grep", c(paste0("\"^",sequenceName,"\\s",position,"\\s.*\\s",ref,"\\s\""), vcfFilePath), stdout=TRUE, stderr=NULL )

      if (length(vcfLine) == 1) # if there is one line in vcfLine
      { 
        lineContents <- strsplit(vcfLine, "\\s")[[1]]
        aoVal <- lineContents[5] # according to the VCF specification
        roVal <- lineContents[4]

        # sometimes there can be more than one ao value, separated by commas:
        if (grepl(",",aoVal,fixed=TRUE)) {        
          aoCount <- sub(".*;AO=", "", vcfLine) # remove everything preceding the value we want to extract        
          aoCount <- sub(";.*", "", aoCount) # remove everything following the value we want to extract

          # find index of the max AO count, and use that to choose the right AO val
          aoCountOpts <- as.numeric( strsplit(aoCount,",")[[1]] )
          indexOfMaxAO <- which(aoCountOpts==max(aoCountOpts, na.rm=TRUE))
          aoCount <- max(aoCountOpts, na.rm=TRUE)
          
          aoValOpts <- strsplit(aoVal,",")[[1]]
          aoVal <- aoValOpts[indexOfMaxAO]
          
          #message(paste0("RS rowNum ",rowNum,": max=",ao,", ao=",r_depthDetailed[3]) )
        }  

        pipelineCol[rowNum] <- vcfDir
        break; # we assume that if we've found data for this CHROM + POS in one vcf file we won't find it the others
      } 
      else if (length(vcfLine) > 1) {
        write(paste0("WARNING!!!: The vcf query for ",sequenceName," ",position," returned more than one line from ",vcfFileNames[whichAO]," in ",vcfDir,"."), stdout())
      }
    } # end inner for-loop that iterates through the vcf file directories and queries the vcf files

    if (!is.na(roVal)) {
      refCol[rowNum] <- roVal
    } 
    if (!is.na(aoVal)) {
      altCol[rowNum] <- aoVal
    }
    # We only queried the AO count if there was more than one AO
    if (!is.na(aoCount)) {
      altCountCol[rowNum] <- aoCount
      relativeCountCol[rowNum] <- (aoCount + report[rowNum, roCountColNames[whichAO]]) / report[rowNum, dpColNames[whichAO]]
    }

  } # end for loop that iterates through the rows of the input file

  # derive ro column name and pipeline column name from the ao column name (these columns are entirely new)
  roColName <- sub("alt","ref",aoColNames[whichAO])
  pipelineColName <- sub("alt","pipe",aoColNames[whichAO])
  if (roColName == aoColNames[whichAO]) {
    roColName <- paste0(roColName,"_ref")
  }
  if (pipelineColName == aoColNames[whichAO]) {
    pipelineColName <- paste0(pipelineColName,"_pipe")
  }

  result <- data.frame(refCol, altCol, altCountCol, relativeCountCol, pipelineCol)
  colnames(result) <- c(roColName, aoColNames[whichAO], aoCountColNames[whichAO], relativeCountColNames[whichAO], pipelineColName)
  result # the return value of this iteration/parallel process.  

} # end parallel processing loop

options(warn=0) # reenable warnings



write(paste0("Combining results into final report."), stdout())

for (i in 1:length(parallelResults)) {
  report[, aoColNames[i]] <- parallelResults[[i]][2]
  report[, aoCountColNames[i]] <- parallelResults[[i]][3]
  report[, relativeCountColNames[i]] <- parallelResults[[i]][4]

  # insert the new ref and pipe columns into the dataframe before the alt column:
  indexOfAltCol <- which(colnames(report)==aoColNames[i])
  report <- cbind(report[, 1:(indexOfAltCol-1)], parallelResults[[i]][5], parallelResults[[i]][1], report[, indexOfAltCol:ncol(report) ] )
}


write(paste0("================================================"), stdout())
write(paste0("Writing final output csv's. "), stdout())

write.csv(report, paste(path.expand(pathToInputFile), outputFileName,sep="/"), row.names=FALSE)

write(paste0("FINISHED."), stdout())




