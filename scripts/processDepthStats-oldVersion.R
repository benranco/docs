args <- commandArgs(trailingOnly = TRUE)


message("Running processDepthStats.R.")

#path <- "~/Desktop/SNPpipeline"
options(stringsAsFactors = FALSE, warn = 1)



path <- "/home/benrancourt/Desktop/junjun/33sampleRNA/reports-pipe02-done"
mafReportName <- "MAF_cutoff_report.csv"
depthDetailedReportName <- "MAF_cutoff_report_depth_detailed.csv"
outputReportName <- "coverages_test.csv"

# Here are the definitions from the VCF file of the depth data we've extracted:
#   DP  = "Total read depth at the locus"
#   DPB = "Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype"
#   RO  = "Reference allele observation count, with partial observations recorded fractionally"
#   AO  = "Alternate allele observations, with partial observations recorded fractionally"

# SNPid-pos-ref   
# 
# dp
# dpb        
# ro ref       
# ro count      
# ro frequency     = ro / ro + ao
# ao ?        
# ao count      
# ao frequency     = ao / ro + ao
# 
# coverage        
# variant1        
# v1coverage      
# v1frequency     
# variant2        
# v2coverage      
# v2frequency     
# 
# Rv1freq/Sv1freq   
# Rv2req/Sv2freq
# 
# R ao matches S ao?
#
# R:
# SNPid coverage  variant1  v1coverage v1frequency  variant2  v2coverage v2frequency
# S:
# SNPid coverage  variant1  v1coverage v1frequency  variant2  v2coverage v2frequency
# Also:
# Rv1freq/Sv1freq   Rv2req/Sv2freq
# 
# v1frequency = v1coverage / v1coverage + v2coverage


#########

mafReport <- read.csv(paste(path.expand(path),mafReportName,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

depthDetailedReport <- read.csv(paste(path.expand(path),depthDetailedReportName,sep="/"),header=TRUE, check.names=FALSE)

report <- mafReport[ , 1:3]
snpID <- paste(report$CHROM,report$POS,report$REF,sep="-")

R_variants <- mafReport[ , "Bigfile_RS.tab"]
S_variants <- mafReport[ , "Bigfile_SC.tab"]

numRows <- nrow(mafReport)

R_coverage <- numeric(numRows)
R_variant1 <- character(numRows)
R_variant2 <- character(numRows)
R_v1coverage <- numeric(numRows)
R_v2coverage <- numeric(numRows)
R_v1frequency <- numeric(numRows)
R_v2frequency <- numeric(numRows)

S_coverage <- numeric(numRows)
S_variant1 <- character(numRows)
S_variant2 <- character(numRows)
S_v1coverage <- numeric(numRows)
S_v2coverage <- numeric(numRows)
S_v1frequency <- numeric(numRows)
S_v2frequency <- numeric(numRows)

v1ratio_RtoS <- numeric(numRows)
v2ratio_RtoS <- numeric(numRows)


# globally supress warnings, because of likely "NAs introduced by coercion"
oldWarningLevel <- getOption("warn")
options(warn = -1)

for (rowNum in 1:numRows) {

  # R data:

  # extract variant 1 and 2
  r_allele <- strsplit(mafReport[rowNum, "Bigfile_RS.tab"],"/")[[1]]
  R_variant1[rowNum] <- r_allele[1]
  R_variant2[rowNum] <- r_allele[2]

  # extract ro and ao
  r_depthDetailed <- strsplit(depthDetailedReport[rowNum, "Bigfile_RS-DP;DPB;AO;RO"],";")[[1]]
  R_coverage[rowNum] <- as.numeric(r_depthDetailed[1])
  R_v1coverage[rowNum] <- as.numeric(r_depthDetailed[4])
  R_v2coverage[rowNum] <- as.numeric(r_depthDetailed[3])

  # calc v1 and v2 frequencies
  R_v1frequency[rowNum] <- round( R_v1coverage[rowNum] / (R_v1coverage[rowNum] + R_v2coverage[rowNum]), digits=3)
  R_v2frequency[rowNum] <- round( R_v2coverage[rowNum] / (R_v1coverage[rowNum] + R_v2coverage[rowNum]), digits=3)

  # S data:

  # extract variant 1 and 2
  s_allele <- strsplit(mafReport[rowNum, "Bigfile_SC.tab"],"/")[[1]]
  S_variant1[rowNum] <- s_allele[1]
  S_variant2[rowNum] <- s_allele[2]

  # extract ro and ao
  s_depthDetailed <- strsplit(depthDetailedReport[rowNum, "Bigfile_SC-DP;DPB;AO;RO"],";")[[1]]
  S_coverage[rowNum] <- as.numeric(s_depthDetailed[1])
  S_v1coverage[rowNum] <- as.numeric(s_depthDetailed[4])
  S_v2coverage[rowNum] <- as.numeric(s_depthDetailed[3])

  # calc v1 and v2 frequencies
  S_v1frequency[rowNum] <- round( S_v1coverage[rowNum] / (S_v1coverage[rowNum] + S_v2coverage[rowNum]), digits=3)
  S_v2frequency[rowNum] <- round( S_v2coverage[rowNum] / (S_v1coverage[rowNum] + S_v2coverage[rowNum]), digits=3)


  # ratios:
  v1ratio_RtoS[rowNum] <- round( R_v1frequency[rowNum] / S_v1frequency[rowNum], digits=3)
  v2ratio_RtoS[rowNum] <- round( R_v2frequency[rowNum] / S_v2frequency[rowNum], digits=3)

} # end for-loop

# globally return warning level to normal
options(warn = oldWarningLevel)

report <- cbind(report, snpID, R_variants, R_coverage, R_variant1, R_v1coverage, R_v1frequency, R_variant2, R_v2coverage, R_v2frequency, S_variants, S_coverage, S_variant1, S_v1coverage, S_v1frequency, S_variant2, S_v2coverage, S_v2frequency, v1ratio_RtoS, v2ratio_RtoS)


#########


write.csv(report, paste(path.expand(path), outputReportName, sep = "/"), row.names=FALSE)

message("...done.")
