args <- commandArgs(trailingOnly = TRUE)


message("Running combineRSandSSfor33sampleTrinityData.R.")

#path <- "~/Desktop/SNPpipeline"
options(stringsAsFactors = FALSE, warn = 1)


######################################################
# Input parameters:

path <- "/home/benrancourt/Desktop/junjun/33sampleRNA/Trinity-SNP-set5"

#RSReportName <- "RS1-first10000.csv"
#SSReportName <- "SS1-first10000.csv"
#outputReportName <- "RS1-SS1-combined-first10000.txt"

#RSReportName <- "RS1-variants2,963,462.txt"
#SSReportName <- "SS1-variants2,837,080.txt"
#outputReportName <- "RS1-SS1-combined.txt"

#RSReportName <- "RS2-variants1,466,720.txt"
#SSReportName <- "SS2-variants1,267406.txt"
#outputReportName <- "RS2-SS2-combined.txt"

#RSReportName <- "RS3-variants2,291,276.txt"
#SSReportName <- "SS3-variants2,045,667.txt"
#outputReportName <- "RS3-SS3-combined.txt"

#RSReportName <- "RS4-variants2,099,778.txt"
#SSReportName <- "SS4-variants1,740,149.txt"
#outputReportName <- "RS4-SS4-combined.txt"

RSReportName <- "RS5-variants514,800.txt"
SSReportName <- "SS5-variants517,226.txt"
outputReportName <- "RS5-SS5-combined.txt"


# The character used as the field delimiter in the .csv file. 
# Use "\t" for tab delimited, or "," for comma delimited fields.
fieldDelimiter <- "\t" 


######################################################
# What this script does, but maybe in a slightly different way:
# - read both tables in
# - only keep rows whose "Mapping", "Reference Position" and "Reference" match between tables, and only of
#   "Variant Type" SNV.
# For each table:
# - keep the "Coverage" column as is
# - check if "Variant #1" matches "Reference"
#   - if yes (default): 
#     - copy "Allele variants" as is
#     - copy "Count of #1" as "Count of Ref", and "Count of #2" as "Count of Alt"
#   - if not, check if "Variant #2" matches "Reference", if so:
#     - set the row's value for "swapped_alleles" to "swapped"
#     - create new value for "Allele variants" by combining "Variant #2"/"Variant #1"
#     - copy "Count of #1" as "Count of Alt", and "Count of #2" as "Count of Ref"
#   - if "Variant #2" doesn't match "Reference" either:
#     - set the row's value for "swapped_alleles" to "no_reference"
#     - copy "Allele variants" as is
#     - copy "Count of #1" as "Count of Ref", and "Count of #2" as "Count of Alt"

######################################################
# Execution code:

RSReport <- read.csv(paste(path.expand(path),RSReportName,sep="/"),header=TRUE, sep=fieldDelimiter, check.names=FALSE)  # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

SSReport <- read.csv(paste(path.expand(path),SSReportName,sep="/"),header=TRUE, sep=fieldDelimiter, check.names=FALSE)  # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.


mergeCols <- c("Mapping","Reference Position","Reference","Variant type")
report <- merge(x=RSReport, y=SSReport, all=FALSE, by=mergeCols)

# After the above merge, all RS column names will be suffixed by ".x" and all SS column names will be
# suffixed by ".y"

# remove all rows whose "Variant type" != "SNV"
onlyKeepSNVVariantType <- report[, "Variant type"] == "SNV"
report <- report[ onlyKeepSNVVariantType, ]

# remove the unnecessary " mapping" from the sequence ids
report$Mapping <- gsub(" mapping", "", report$Mapping, fixed=TRUE)
colnames(report)[1] <- "seqId"

numRows <- nrow(report)

RS_coverage <- report[ , "Coverage.x"]
RS_num_variants <- report[ , "Variants.x"]
RS_allele_variants <- report[ , "Allele variants.x"]
RS_count_ref <- report[ , "Count of #1.x"]
RS_count_alt <- report[ , "Count of #2.x"]
RS_swapped_alleles <- character(numRows)

SS_coverage <- report[ , "Coverage.y"]
SS_num_variants <- report[ , "Variants.y"]
SS_allele_variants <- report[ , "Allele variants.y"]
SS_count_ref <- report[ , "Count of #1.y"]
SS_count_alt <- report[ , "Count of #2.y"]
SS_swapped_alleles <- character(numRows)

for (rowNum in 1:numRows) {

  # process RS data

  if (report[rowNum, "Variant #1.x"] != report[rowNum, "Reference"]) {
    
    if (report[rowNum, "Variant #2.x"] == report[rowNum, "Reference"]) {
      RS_swapped_alleles[rowNum] <- "swapped"
      RS_allele_variants[rowNum] <- paste(report[rowNum, "Variant #2.x"], report[rowNum, "Variant #1.x"], sep="/")
      RS_count_ref[rowNum] <- report[rowNum, "Count of #2.x"]
      RS_count_alt[rowNum] <- report[rowNum, "Count of #1.x"]

      # just in case there are more than 2 variants
      if (report[rowNum, "Variants.x"] > 2) {
        vars <- strsplit(report[rowNum, "Allele variants.x"], split="/")[[1]]
        vars <- paste(vars[3:length(vars)],collapse="/")
        RS_allele_variants[rowNum] <- paste(RS_allele_variants[rowNum], vars, sep="/")
      }
    }
    else {
      RS_swapped_alleles[rowNum] <- "no_reference"
    }
  }
  else {
    RS_swapped_alleles[rowNum] <- "-"
  }

  # process SS data

  if (report[rowNum, "Variant #1.y"] != report[rowNum, "Reference"]) {
    
    if (report[rowNum, "Variant #2.y"] == report[rowNum, "Reference"]) {
      SS_swapped_alleles[rowNum] <- "swapped"
      SS_allele_variants[rowNum] <- paste(report[rowNum, "Variant #2.y"], report[rowNum, "Variant #1.y"], sep="/")
      SS_count_ref[rowNum] <- report[rowNum, "Count of #2.y"]
      SS_count_alt[rowNum] <- report[rowNum, "Count of #1.y"]

      # just in case there are more than 2 variants
      if (report[rowNum, "Variants.y"] > 2) {
        vars <- strsplit(report[rowNum, "Allele variants.y"], split="/")[[1]]
        vars <- paste(vars[3:length(vars)],collapse="/")
        SS_allele_variants[rowNum] <- paste(SS_allele_variants[rowNum], vars, sep="/")
      }
    }
    else {
      SS_swapped_alleles[rowNum] <- "no_reference"
    }
  }
  else {
    SS_swapped_alleles[rowNum] <- "-"
  }
  
} # end for-loop


report <- cbind(report[ , 1:4], RS_coverage, RS_num_variants, RS_allele_variants, RS_count_ref, RS_count_alt, RS_swapped_alleles, SS_coverage, SS_num_variants, SS_allele_variants, SS_count_ref, SS_count_alt, SS_swapped_alleles)

write.table(report, paste(path.expand(path), outputReportName, sep = "/"), row.names=FALSE, sep=fieldDelimiter)

message("...done.")
