############################################################################################
# This script is a revision of convert_csv_to_linkage.R, to convert a .csv file in chi 
# square root table format into a .linkage file. The difference between this and 
# convert_csv_to_linkage.R is that the input .csv is oriented with samples being represented
# in rows, and SNPs being represented by columns. I.e. each row is going to be really long,
# which is not recommended in a .csv file. This script should hardly, if ever, need to be 
# used. This script is also not built into the LepMAP2pipeline, so if plan to use the 
# pipeline, you should just transpose your input .csv file to list the samples as columns 
# and the SNPs as rows before running the pipeline, in which case you won't need this script.
#
# Input .csv file format:
# The chi square root table consists of one row listing SNP ids, and a row for 
# each sample. The data in the chi square root table is represented as either "H", "A", or
# something ("-" or "NA" for example, you can specify which you use in the 
# missingDataIndicator input parameter below) to represent missing or irrelevant data. The 
# most frequent genotype in an SNP is represented as H and the second most frequent type 
# is represented as A. 
#
#
# Information on .linkage file format (for use with Lep-MAP2 software):
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
#
#
############################################################################################
# Input Parameters:

#args <- commandArgs(trailingOnly = TRUE)
#path <- args[1]

path <- "/home/benrancourt/Downloads/r38-753+5092"

#MAF_CUTOFF <- args[2]
#MAF_CUTOFF <- as.double(MAF_CUTOFF)

options(stringsAsFactors = FALSE, warn = 1)



filenameMinusExtension <- "Data-753+5,092"

familyName <- "r38"


reportc <- read.csv(paste(path.expand(path),paste0(filenameMinusExtension,".csv"), sep="/"), check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.
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


reportLinkageGenotypes <- reportc[, -1] # remove the sample id column
message(paste0("ncols: ",ncol(reportLinkageGenotypes)))
reportLinkageGenotypes <- rbind(parent1 = c("A"), parent2 = c("H"), reportLinkageGenotypes) # add two samples to use as parents
message(paste0("ncols: ",ncol(reportLinkageGenotypes)))
#reportLinkageGenotypes <- t(reportLinkageGenotypes) # transpose the report (so it's columns are now rows)
message(paste0("nrows: ",nrow(reportLinkageGenotypes)))

reportLinkageGenotypes[reportLinkageGenotypes=="H"] <- "1 2"
reportLinkageGenotypes[reportLinkageGenotypes=="A"] <- "1 1"
reportLinkageGenotypes[reportLinkageGenotypes=="B"] <- "2 2"
reportLinkageGenotypes[is.na(reportLinkageGenotypes)] <- "0 0"
reportLinkageGenotypes[reportLinkageGenotypes=="-"] <- "0 0" # in case NA "-" has already been substituted with "-"

reportLinkage <- cbind(family = c(familyName), id = c(paste0("S",(1:nrow(reportLinkageGenotypes))-2)), fatherId = c("P1"), motherId = c("P2"), gender = c(0), affectionStatus = c(0), reportLinkageGenotypes)
reportLinkage[1,2:5] <- c("P1","0","0","1") # change id from S-1 to P1, no parents, male
reportLinkage[2,2:5] <- c("P2","0","0","2") # change id from S0 to P2, no parents, female

write.table(reportLinkage, file= paste(path.expand(path),paste0(filenameMinusExtension,".linkage"), sep="/"), append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

message("report_gen part 2 complete")
