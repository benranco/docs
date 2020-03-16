args <- commandArgs(trailingOnly = TRUE)


message("Running processDepthStats.R.")

#path <- "~/Desktop/SNPpipeline"
options(stringsAsFactors = FALSE, warn = 1)



path <- "/home/benrancourt/Desktop/junjun/33sampleRNA/reports-pipe05-done"
mafReportName <- "MAF_cutoff_report.csv"
depthDetailedReportName <- "MAF_cutoff_report_depth_detailed.csv"
outputReportName <- "33sampleRNA-pipe05-coverage_frequencies.csv"

# Here are the definitions from the VCF file of the depth data we've extracted:
#   DP  = "Total read depth at the locus"
#   DPB = "Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype"
#   RO  = "Reference allele observation count, with partial observations recorded fractionally"
#   AO  = "Alternate allele observations, with partial observations recorded fractionally"

# SNPid-pos-ref   
# 
# alt
# dp
# dpb        
# ro ref           
# ro count      
# ro frequency     = ro / ro + ao
# ao ?        
# ao count      
# ao frequency     = ao / ro + ao
# 
# alt
# dp
# dpb        
# ro ref           
# ro count      
# ro frequency     = ro / ro + ao
# ao ?        
# ao count      
# ao frequency     = ao / ro + ao
# 
# ro_ratio_RtoS    = R_ro_freq/S_ro_freq   
# ao_ratio_RtoS    = R_ao_freq/S_ao_freq
# 
# TODO: add column to check if R ao matches S ao?
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


depthDetailedReport <- read.csv(paste(path.expand(path),depthDetailedReportName,sep="/"),header=TRUE, check.names=FALSE)  # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

mafReport <- read.csv(paste(path.expand(path),mafReportName,sep="/"),header=TRUE, check.names=FALSE)  # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.


# remove all rows with indels (ie where nchar of REf is > 1) from both maf and depth reports
toRemoveIndels <- nchar(mafReport$REF) == 1
mafReport <- mafReport[ toRemoveIndels, ]
depthDetailedReport <- depthDetailedReport[ toRemoveIndels, ]



mafDataColumns <- mafReport[ , 4:ncol(mafReport)]

report <- depthDetailedReport[ , 1:3]
snpID <- paste(report$CHROM,report$POS,report$REF,sep="-")

numRows <- nrow(depthDetailedReport)

R_alt <- character(numRows)
R_dp <- numeric(numRows)
R_dpb <- numeric(numRows)
R_ro_count <- numeric(numRows)
R_ao_count <- numeric(numRows)
R_ro_frequency <- numeric(numRows)
R_ao_frequency <- numeric(numRows)

S_alt <- character(numRows)
S_dp <- numeric(numRows)
S_dpb <- numeric(numRows)
S_ro_count <- numeric(numRows)
S_ao_count <- numeric(numRows)
S_ro_frequency <- numeric(numRows)
S_ao_frequency <- numeric(numRows)

ro_ratio_RtoS <- numeric(numRows)
ao_ratio_RtoS <- numeric(numRows)


# globally supress warnings, because of likely "NAs introduced by coercion"
oldWarningLevel <- getOption("warn")
options(warn = -1)

for (rowNum in 1:numRows) {

  # R data:

  # TODO: get the ALT
  R_alt[rowNum] <- "?"

  # extract ro and ao
  r_depthDetailed <- strsplit(depthDetailedReport[rowNum, "Bigfile_RS-DP;DPB;AO;RO"],";")[[1]]
  R_dp[rowNum] <- as.numeric(r_depthDetailed[1])
  R_dpb[rowNum] <- as.numeric(r_depthDetailed[2])
  R_ro_count[rowNum] <- as.numeric(r_depthDetailed[4])
  # sometimes there can be more than one ao value, separated by commas:
  if (grepl(",",r_depthDetailed[3],fixed=TRUE)) {
    ao <- as.numeric( strsplit(r_depthDetailed[3],",")[[1]] )
    ao <- max(ao, na.rm=TRUE)
    R_ao_count[rowNum] <- max(ao, na.rm=TRUE)
    #message(paste0("RS rowNum ",rowNum,": max=",ao,", ao=",r_depthDetailed[3]) )
  }
  else {
    R_ao_count[rowNum] <- as.numeric(r_depthDetailed[3])
  }

  # calc v1 and v2 frequencies
  R_ro_frequency[rowNum] <- round( R_ro_count[rowNum] / (R_ro_count[rowNum] + R_ao_count[rowNum]), digits=3)
  R_ao_frequency[rowNum] <- round( R_ao_count[rowNum] / (R_ro_count[rowNum] + R_ao_count[rowNum]), digits=3)

  # S data:

  # TODO: get the ALT
  S_alt[rowNum] <- "?"

  # extract ro and ao
  s_depthDetailed <- strsplit(depthDetailedReport[rowNum, "Bigfile_SC-DP;DPB;AO;RO"],";")[[1]]
  S_dp[rowNum] <- as.numeric(s_depthDetailed[1])
  S_dpb[rowNum] <- as.numeric(s_depthDetailed[2])
  S_ro_count[rowNum] <- as.numeric(s_depthDetailed[4])
  # sometimes there can be more than one ao value, separated by commas:
  if (grepl(",",s_depthDetailed[3],fixed=TRUE)) {
    ao <- as.numeric( strsplit(s_depthDetailed[3],",")[[1]] )
    ao <- max(ao, na.rm=TRUE)
    S_ao_count[rowNum] <- max(ao, na.rm=TRUE)
    #message(paste0("SC rowNum ",rowNum,": max=",ao,", ao=",s_depthDetailed[3]) )
  }
  else {
    S_ao_count[rowNum] <- as.numeric(s_depthDetailed[3])
  }

  # calc v1 and v2 frequencies
  S_ro_frequency[rowNum] <- round( S_ro_count[rowNum] / (S_ro_count[rowNum] + S_ao_count[rowNum]), digits=3)
  S_ao_frequency[rowNum] <- round( S_ao_count[rowNum] / (S_ro_count[rowNum] + S_ao_count[rowNum]), digits=3)


  # ratios:
  ro_ratio_RtoS[rowNum] <- round( R_ro_frequency[rowNum] / S_ro_frequency[rowNum], digits=3)
  ao_ratio_RtoS[rowNum] <- round( R_ao_frequency[rowNum] / S_ao_frequency[rowNum], digits=3)

} # end for-loop

# globally return warning level to normal
options(warn = oldWarningLevel)

report <- cbind(report, snpID, mafDataColumns, R_alt, R_dp, R_dpb, R_ro_count, R_ao_count, R_ro_frequency, R_ao_frequency, S_alt, S_dp, S_dpb, S_ro_count, S_ao_count, S_ro_frequency, S_ao_frequency, ro_ratio_RtoS, ao_ratio_RtoS)


#########


write.csv(report, paste(path.expand(path), outputReportName, sep = "/"), row.names=FALSE)

message("...done.")
