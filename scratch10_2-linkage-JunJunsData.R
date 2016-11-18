#args <- commandArgs(trailingOnly = TRUE)
#path <- args[1]

path <- "/home/benrancourt/Downloads/"

#MAF_CUTOFF <- args[2]
#MAF_CUTOFF <- as.double(MAF_CUTOFF)

options(stringsAsFactors = FALSE, warn = 1)

message("running report generation part 2")
if(!require(seqinr))
{
  install.packages('seqinr', repos='http://cran.us.r-project.org')
}

if(!require(VariantAnnotation))
{
  source("https://bioconductor.org/biocLite.R")
  biocLite("VariantAnnotation")
}

library(VariantAnnotation)
library(seqinr)

filenameMinusExtension <- "fromJunJun-toConvertToLinkage-R28-S91-M28K-Ben"
reportc <- read.csv(paste0(path,filenameMinusExtension,".csv"))
#reportc <- reportc[, -1]

if(!("COMBINED" %in% colnames(reportc)))
{
  s <- 4
}else
{
  s <- 5
}

# #######################################################################
message("converting the chi square report to .linkage format")

# Information on .linkage file format (for use with Le-MAP2 software):
#
# Data is tab-delimited.
#
# First 6 Columns contain "pedigree information":
# 1. Family ID (can be alphanumeric)
# 2. Individual ID (must be unique within family, can be alphanumeric)
# 3. Father ID (0 if father is not in family)
# 4. Mother ID (0 if mother is not in family)
# 5. Gender (0=unknown, 1=male, 2=female)
# 6. Affection status (0=unknown, 1=unaffected, 2=affected)
#
# Columns 7 and onward describe the phenotype data, separated by tabs. There 
# are four different types of phenotype data supported by the LINKAGE format 
# (Numbered Alleles, Binary Factors, Affection Status, Quantitative Traits), 
# and the Lep-MAP2 documentation uses Numbered Alleles.
#
# With Numbered Alleles, each genotype is represented as a pair of numbers 
# (eg: 1 2     2 2   1 1    1 2). Each number represents an allele, and 0 
# represents an unknown allele.
#
# For our purposes (input for the Lep-MAP2 software), we are setting:
# H = "1 2", A = "1 1", B = "2 2", NA = "0 0".
#
# Lep-MAP2 documentation: https://sourceforge.net/p/lepmap2/wiki/browse_pages/
#
# Official LINKAGE file format documentation available:
# http://www.jurgott.org/linkage/LinkagePC.html#__RefHeading__137_1806185151
# http://www.jurgott.org/linkage/LinkageUser.pdf

reportLinkageGenotypes <- reportc[, -1] # remove the marker id column
reportLinkageGenotypes <- cbind(parent1 = c("A"), parent2 = c("H"), reportLinkageGenotypes) # add two samples to use as parents
reportLinkageGenotypes <- t(reportLinkageGenotypes) # transpose the report (so it's columns are now rows)

reportLinkageGenotypes[reportLinkageGenotypes=="H"] <- "1 2"
reportLinkageGenotypes[reportLinkageGenotypes=="A"] <- "1 1"
reportLinkageGenotypes[reportLinkageGenotypes=="B"] <- "2 2"
reportLinkageGenotypes[is.na(reportLinkageGenotypes)] <- "0 0"
reportLinkageGenotypes[reportLinkageGenotypes=="-"] <- "0 0" # in case NA "-" has already been substituted with "-"

reportLinkage <- cbind(family = c("Try-1-S91"), id = c(paste0("S",(1:nrow(reportLinkageGenotypes))-2)), fatherId = c("P1"), motherId = c("P2"), gender = c(0), affectionStatus = c(0), reportLinkageGenotypes)
reportLinkage[1,2:5] <- c("P1","0","0","1") # change id from S-1 to P1, no parents, male
reportLinkage[2,2:5] <- c("P2","0","0","2") # change id from S0 to P2, no parents, female

write.table(reportLinkage, file= paste0(path,filenameMinusExtension,".linkage"), append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

message("report_gen part 2 complete")
