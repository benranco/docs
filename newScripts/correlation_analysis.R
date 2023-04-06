# Some useful documentation of Correlation and P-value calculation.
#
# This seems like a decent simple introduction to Correlation and P-value: 
# https://dataschool.com/fundamentals-of-analysis/correlation-and-p-value/
#
# Here are some online R and p value calculators:
# R-value calculator:
# https://www.socscistatistics.com/tests/pearson/
# p-value calculator:
# https://www.socscistatistics.com/pvalues/pearsondistribution.aspx
# https://www.danielsoper.com/statcalc/formulas.aspx?id=44
#
#
# Calculating R: 
#
# Calculate the R coefficient for Pearson correlation calculation in R:
# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
# http://www.sthda.com/english/wiki/t-distribution-table
# 
# Using data row two above:
# x <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0)
# y <- c(1,1,1,1,1,1,1,1,1,1,1,0.430769,0.495575,0.836066,0.820755,0.390625,1,0.793103,0.659574,0.78,0.984848,0.730159,0.844828,0.492063,1,1,1,0.971831,1,0.673077,1,1,0.913043,0.697368,1,0.578947)
# 
# To calulate R (it doesn't matter which is x and which is y: 
# cor(x, y,  method = "pearson", use = "na.or.complete") 
#
#
# Calculating t for a Pearson correlation: 
#
# From: https://www.danielsoper.com/statcalc/formulas.aspx?id=44
#
# Calculate the t-value for a Pearson correlation (where n is the number of samples): 
#   t <- R / sqrt( (1-R**2)/(n-2) ) 
# same as:
#   t <- R / sqrt(1-R**2)) * sqrt(n-2)
#
#
# How to manually calculate the p-value given the t value/score (their online calculator using a two-sided test yields correct results):
# https://statkat.com/find-p-value/t-value.php
# https://statkat.com/online-calculators/p-value-given-t-value.php
# https://statkat.com/degrees-of-freedom-t-test.php
#
#
# Calculating the p-value:
#
# One way to calculate the p-value (which this script doesn't use), from:
# https://stats.stackexchange.com/questions/153937/finding-p-value-in-pearson-correlation-in-r
# Calculate the p-value (still need to figure out how to tell it what to do for NA values, see documentation): 
# cor.test(x, y,  method = "pearson")$p.value
#
# Another description of how to calculate p given t in R (I used this method): 
# https://www.statology.org/p-value-of-t-score-r/
# If using a negative t value, I get the proper results using lower.tail=TRUE.
# If using a positive t value, I get the proper results using lower.tail=FALSE.
# This will always work, whether t is positive or negative: 
#   pt(q=abs(t), df=34, lower.tail=FALSE)*2
#
# #######################################################################

options(stringsAsFactors = FALSE, warn = 1)

# Input parameters below:

path <- "."

naChar <- "-"
numMetaCols <- 4

inputRefAllelicFreqsCsv <- "LP36_pipe1_filtered-ratio-first10000.csv"
outputCorrelationAnalysis <- "LP36_pipe1_filtered-correlation-first10000.csv"

inputStatsCsv <- "LP36_pipe1_filtered-stats-first10000.csv"
outputStatsCsv <- "LP36_pipe1_filtered-stats-correlation-first10000.csv"


# Assuming this is the exact order of the sample columns: 
# C10,C11,C12,C13,C1,C3,C4,C5,C6,C7,C8,C9,LP-C2,RN1,RS10,RS1,RS2,RS3,RS4,RS5,RS6,RS7,RS8,RS9,SN1,SN2,SS10,SS1,SS2,SS3,SS4,SS5,SS6,SS7,SS8,SS9
# Assign Res samples with resistance level at 1, and Sus samples with resistance level at 0, and C samples with resistance level at 0: 
#resistanceLevels <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0)
#resistanceLevels <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0)
resistanceLevels  <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0)

# End of input parameters
# #######################################################################

doCorrelation <- function(xRow, yVec) {
  
  R <- cor(as.numeric(xRow), yVec,  method = "pearson", use = "na.or.complete")  # TODO: check on whether to set use = "complete.obs"
  
  numNAs <- sum(is.na(xRow))
  n <- ncol(xRow) - numNAs
  
  t <- R / sqrt( (1-R**2)/(n-2) ) 
  
  # the *2 is for the two-sided test, the abs(t) ensures lower.tail=FALSE will always give the correct result)
  pVal <- pt(q=abs(t), df=(n-2), lower.tail=FALSE)*2  
  
  #message(paste0("R == ",R,", t == ",t,", p == ",pVal))
  
  results <- c(R=R, p=pVal, n=n, t=t)
  results # return results
}


# #######################################################################
message(paste0("Processing --- ",Sys.time()))

refAllelicFreqs <- read.csv(paste(path.expand(path),inputRefAllelicFreqsCsv,sep="/"),header=TRUE, check.names=FALSE, na.strings=naChar)


correlData <- cbind(refAllelicFreqs[, 1:numMetaCols], R=numeric(nrow(refAllelicFreqs)), p=numeric(nrow(refAllelicFreqs)), n=numeric(nrow(refAllelicFreqs)), t=numeric(nrow(refAllelicFreqs)) ) 

dataCols <- refAllelicFreqs[,(numMetaCols+1):ncol(refAllelicFreqs)]

message(paste0("Single row --- ",Sys.time()))
correl <- doCorrelation( dataCols[1,], resistanceLevels)

message(paste0("for-loop --- ",Sys.time()))
for (i in 1:nrow(dataCols)) {
  correl <- doCorrelation( dataCols[i,], resistanceLevels)
  correlData[i,"R"] <- round(correl["R"], digits=6)
  correlData[i,"p"] <- round(correl["p"], digits=6)
  correlData[i,"n"] <- correl["n"]
  correlData[i,"t"] <- round(correl["t"], digits=6)
}


message(paste0("Finished correlation analysis --- ",Sys.time()))

message(paste0("Writing output correlation file --- ",Sys.time()))
write.csv(correlData, paste(path.expand(path), outputCorrelationAnalysis, sep = "/"), row.names=FALSE)

message(paste0("Reading stats file --- ",Sys.time()))
stats <- read.csv(paste(path.expand(path),inputStatsCsv,sep="/"),header=TRUE, check.names=FALSE, na.strings=naChar)

# used to preserve the order of stats after doing the merge
stats$index  <- 1:nrow(stats)

message(paste0("Merging stats data and correlation data --- ",Sys.time()))
newstats <- merge(stats, correlData, by=c("CHROM","POS","REF","ALT"), all=TRUE)

# set the order back to the original order of stats
newstats <- newstats[order(newstats$index), ]
newstats <- newstats[, -which(names(newstats) %in% c("index"))]


message(paste0("Writing merged data --- ",Sys.time()))
write.csv(newstats, paste(path.expand(path), outputStatsCsv, sep = "/"), row.names=FALSE)


