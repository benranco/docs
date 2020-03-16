args <- commandArgs(trailingOnly = TRUE)


message("Running combineRSandSSfor33sampleTrinityData.R.")

#path <- "~/Desktop/SNPpipeline"
options(stringsAsFactors = FALSE, warn = 1)


######################################################
# Input parameters:

path <- "/home/benrancourt/Desktop/junjun/XPGWAS/Rust-GWAS"

#report1Name <- "avcr2-AA3_variants-375,673.csv"
#report2Name <- "vcr2_ mapping variants-453,143.csv"
#outputReportName <- "XPGWAS-avcr2-vcr2-combined.csv"

report1Name <- "avc2-AA3_ vari646,258-ref24,535.csv"
report2Name <- "vcr2-AA5-vari758,163-ref24,535.csv"
outputReportName <- "Rust-GWAS-avc2AA3-vcr2AA5-combined.csv"


# The character used as the field delimiter in the .csv file. 
# Use "\t" for tab delimited, or "," for comma delimited fields.
fieldDelimiter <- "," 


######################################################
# What this script does, but maybe in a slightly different way:
# - read both tables in
# - only keep rows whose "Mapping", "Reference Position" and "Reference" match between tables, and only of
#   "Variant Type" SNV.
# - create new id from those
# For each XPGWAS table:
# - if number of variants > 2, remove row
# - in both tables, check to make sure variant1 == REF, and if not, switch the variants
#   - if necessary, switch all other columns related to variant 1 and 2
# - keep the Coverage column as is
# - if, after switching to variant1 == REF, if variant2 doesn't match between tables, remove row
#
######################################################
# Execution code:

report1 <- read.csv(paste(path.expand(path),report1Name,sep="/"),header=TRUE, sep=fieldDelimiter, check.names=FALSE)  # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

report2 <- read.csv(paste(path.expand(path),report2Name,sep="/"),header=TRUE, sep=fieldDelimiter, check.names=FALSE)  # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.



mergeCols <- c("Mapping","Reference Position","Reference","Variant type")
report <- merge(x=report1, y=report2, all=FALSE, by=mergeCols)

# After the above merge, all report1 column names will be suffixed by ".x" and all report2 column names will be
# suffixed by ".y"

# remove all rows whose "Variant type" != "SNV"
onlyKeepSNVVariantType <- report[, "Variant type"] == "SNV"
report <- report[ onlyKeepSNVVariantType, ]

# remove the unnecessary " mapping" from the sequence ids
report$Mapping <- gsub(" mapping", "", report$Mapping, fixed=TRUE)
#colnames(report)[1] <- "seqId"


# for each dataset, remove all rows whose number of variants is > 2
onlyKeepRowsWithTwoOrLessVariants <- report[, "Variants.x"] <= 2
report <- report[ onlyKeepRowsWithTwoOrLessVariants, ]
onlyKeepRowsWithTwoOrLessVariants <- report[, "Variants.y"] <= 2
report <- report[ onlyKeepRowsWithTwoOrLessVariants, ]

numRows <- nrow(report)

r1_coverage <- report[ , "Coverage.x"]
r1_num_variants <- report[ , "Variants.x"]
r1_variant_ref <- report[ , "Variant #1.x"]
r1_variant_alt <- report[ , "Variant #2.x"]
r1_frequency_ref <- report[ , "Frequency of #1.x"]
r1_frequency_alt <- report[ , "Frequency of #2.x"]
r1_count_ref <- report[ , "Count of #1.x"]
r1_count_alt <- report[ , "Count of #2.x"]
r1_swapped_variants <- character(numRows)

r2_coverage <- report[ , "Coverage.y"]
r2_num_variants <- report[ , "Variants.y"]
r2_variant_ref <- report[ , "Variant #1.y"]
r2_variant_alt <- report[ , "Variant #2.y"]
r2_frequency_ref <- report[ , "Frequency of #1.y"]
r2_frequency_alt <- report[ , "Frequency of #2.y"]
r2_count_ref <- report[ , "Count of #1.y"]
r2_count_alt <- report[ , "Count of #2.y"]
r2_swapped_variants <- character(numRows)


# For each XPGWAS table:
# - if number of variants > 2, remove row
# - in both tables, check to make sure variant1 == REF, and if not, switch the variants
#   - if necessary, switch all other columns related to variant 1 and 2
# - keep the Coverage column as is
# - if, after switching to variant1 == REF, if variant2 doesn't match between tables, remove row

otherRowsToRemove <- logical(numRows)

for (rowNum in 1:numRows) {

  # process r1 data

  if (report[rowNum, "Variant #1.x"] != report[rowNum, "Reference"]) {
    
    if (report[rowNum, "Variant #2.x"] == report[rowNum, "Reference"]) {
      r1_swapped_variants[rowNum] <- "swapped"
      r1_variant_ref[rowNum] <- report[rowNum, "Variant #2.x"]
      r1_variant_alt[rowNum] <- report[rowNum, "Variant #1.x"]
      r1_frequency_ref[rowNum] <- report[rowNum, "Frequency of #2.x"]
      r1_frequency_alt[rowNum] <- report[rowNum, "Frequency of #1.x"]
      r1_count_ref[rowNum] <- report[rowNum, "Count of #2.x"]
      r1_count_alt[rowNum] <- report[rowNum, "Count of #1.x"]
      
    }
    else {
      otherRowsToRemove[rowNum] <- TRUE
    }
  }
  else {
    r1_swapped_variants[rowNum] <- "-"
  }

  # process r2 data

  if (report[rowNum, "Variant #1.y"] != report[rowNum, "Reference"]) {
    
    if (report[rowNum, "Variant #2.y"] == report[rowNum, "Reference"]) {
      r2_swapped_variants[rowNum] <- "swapped"
      r2_variant_ref[rowNum] <- report[rowNum, "Variant #2.y"]
      r2_variant_alt[rowNum] <- report[rowNum, "Variant #1.y"]
      r2_frequency_ref[rowNum] <- report[rowNum, "Frequency of #2.y"]
      r2_frequency_alt[rowNum] <- report[rowNum, "Frequency of #1.y"]
      r2_count_ref[rowNum] <- report[rowNum, "Count of #2.y"]
      r2_count_alt[rowNum] <- report[rowNum, "Count of #1.y"]

    }
    else {
      otherRowsToRemove[rowNum] <- TRUE
    }
  }
  else {
    r2_swapped_variants[rowNum] <- "-"
  }
  
  # if, after switching to variant1 == REF, variant2 doesn't match between tables, remove row
  if (r1_variant_alt[rowNum] != r2_variant_alt[rowNum]) {
    otherRowsToRemove[rowNum] <- TRUE
  }

} # end for-loop

id <- paste(report[,"Mapping"], report[,"Reference Position"], report[,"Reference"], sep="-")

report <- cbind(id, report[ , 1:3], r1_coverage, r1_num_variants, r1_variant_ref, r1_variant_alt, r1_count_ref, r1_count_alt, r1_frequency_ref, r1_frequency_alt, r1_swapped_variants, r2_coverage, r2_num_variants, r2_variant_ref, r2_variant_alt, r2_count_ref, r2_count_alt, r2_frequency_ref, r2_frequency_alt, r2_swapped_variants)

report <- report[!otherRowsToRemove, ]

write.table(report, paste(path.expand(path), outputReportName, sep = "/"), row.names=FALSE, sep=fieldDelimiter)

message("...done.")
