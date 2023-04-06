args <- commandArgs(trailingOnly = TRUE)
path <- "." #args[1]

#MAF_CUTOFF <- args[2]
MAF_CUTOFF <- 0.05 #as.double(MAF_CUTOFF)

#GENERATE_CHI_SQ_REPORT <- args[3]
GENERATE_CHI_SQ_REPORT <- 0 #as.integer(GENERATE_CHI_SQ_REPORT)

#GENERATE_PROBABILITY_REPORT <- args[4]
GENERATE_PROBABILITY_REPORT <- 0 #as.integer(GENERATE_PROBABILITY_REPORT)

#HAPLOID_OR_DIPLOID <- args[5] # 1 == haploid, 2 == diploid
HAPLOID_OR_DIPLOID <- 2 #as.integer(HAPLOID_OR_DIPLOID)

HAS_INDELS <- TRUE  
#if (as.integer(args[6]) == 1) {
  #HAS_INDELS <- FALSE
#}

# (1) run the full pipeline
# (2) just process the data and do not generate reports (ie. just run the first half of the pipeline), 
# (3) just generate reports based on data that has already been processed by the first half of the pipeline (ie. just run the second half of the pipeline assuming the first half has already been run).
# (4) just generate reports beginning AFTER the filled_report.csv, assuming it has already been generated.
WHAT_TO_RUN <- 4 #args[7]


report <- data.frame()
options(stringsAsFactors = FALSE, warn = 1)

message(paste0("running report generation part 2 --- ",Sys.time()))


if(!require(seqinr)) {
  install.packages('seqinr', repos='http://cran.us.r-project.org')
}
library(seqinr)

if(GENERATE_PROBABILITY_REPORT) {
  # Used in report_gen_p2.R
  # This install proceedure is tested to work in R 4.2.1, specifically: BiocManager::install(version = "3.15"
  # It may require keyboard input in response to a question (select "a" for all).
  # If there are errors (you can ignore warnings) during installation, it might be because some required 
  # linux libraries need to be installed first. You can figure that out be reading through the error messages.
  # For example, I needed to install libxml2-devel and openssl-devel via dnf.  
  if(!require(VariantAnnotation)) {
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos='http://cran.us.r-project.org')
    }
    BiocManager::install(version = "3.15")
    BiocManager::install("VariantAnnotation")
  }
  library(VariantAnnotation)
}


# #######################################################################
if(WHAT_TO_RUN == 4) {
  # assume filled_report already exists
  report <- read.csv(paste(path.expand(path),"reports","filled_report.csv",sep="/"),header=TRUE, check.names=FALSE)
} else {
  # concat the files in reporttemp to make the filled_report.
  message(paste0("reading and merging split data --- ",Sys.time()))
  for(single1 in list.files(paste0(path, "/reporttemp"))) {
  
    if(grepl("filled.Rds", single1, fixed=TRUE)) {
      #message(single1)
      sr <- readRDS(paste0(path, "/reporttemp/", single1))
      if(nrow(report) == 0 && ncol(report) == 0) {
        report <- sr
      } else {
        report <- rbind(report, sr)
      }
    }
  }

  write.csv(report, paste(paste(path.expand(path), "reports", sep = "/"), "filled_report.csv", sep = "/"), row.names=FALSE)
}

# #######################################################################
if(!("COMBINED" %in% colnames(report))) {
  startCol <- 4
  s <- 4 # we have both s and startCol because some legacy code used s as a variable name, and searching and replacing it would be tricky. Haven't done that yet.
} else {
  startCol <- 5
  s <- 5
}

# #######################################################################
# Function: isRowNonBiallelic - determines if the given row data has more than two alleles.
# input:  a list of all the unique data items in the report. Each data item is formatted:
#         "allele/allele" or possibly (for legacy purposes) "allele/". An allele is usually
#         only one character long, but in the case where indels were not removed in the first
#         part of the pipeline, an allele can have two or more characters.
# returns: TRUE if the row has more than two alleles, FALSE if it only has two (or one) allele
isRowNonBiallelic <- function(namesOfGenotypesInRow) {

  alleleIndex <- 1
  individualAlleles <- character(length(namesOfGenotypesInRow) * 2) # because, for example "A/A" has two alleles, "A/G" has two. If a type is of format "A/", this function will pretend it is in format "A/A" and count the single A twice. This is okay to do, because the count is only used internally to this function.
  for ( type in namesOfGenotypesInRow ) {
    allelesInType <- strsplit(type, "/")[[1]]
    individualAlleles[alleleIndex] <- allelesInType[1]
    alleleIndex <- alleleIndex + 1
    individualAlleles[alleleIndex] <- ifelse (length(allelesInType) > 1, allelesInType[2], allelesInType[1]) # Because for haploid data, it is in format "A/" instead of "A/A".
    alleleIndex <- alleleIndex + 1
  }
  # if there are 3 or more unique alleles in the row, return TRUE, else return FALSE
  return ( ifelse( length(names(table(individualAlleles, useNA = "no"))) >= 3, TRUE, FALSE ) )
}


# #######################################################################
message(paste0("editing for cutoffs --- ",Sys.time()))
  
report <- as.data.frame(report)

if((ncol(report) > 24 && "COMBINED" %in% colnames(report)) || (ncol(report) > 23 && !("COMBINED" %in% colnames(report)))) {
  doFullEditingForCutoffs <- TRUE
} else {
  doFullEditingForCutoffs <- FALSE
  message("not enough samples for full cutoff editing, just removing non-biallelic rows")
}

# CUSTOM: # keep only those rows whose REF value has only one character in it (i.e. it's not an indel):
report <- report[ (nchar(report$REF) == 1), ]
HAS_INDELS <- FALSE

numrows <- nrow(report)
numColsInReport <- ncol(report)
rowsToRemove <- NULL
rowsToRemove <- logical(numrows) # initializes all to FALSE


numDataCols <- numColsInReport - (startCol - 1)

for (curRow in 1:nrow(report)) {
  rowData <- as.matrix(report[curRow,startCol:numColsInReport]) 
  tabulatedRowData <- table(rowData, useNA = "no")
  tabulatedRowDataNames <- names(tabulatedRowData)

  if( length(tabulatedRowDataNames) <= 1 & doFullEditingForCutoffs)  {
  # remove row if there's only one (or 0) type of non-NA data
    rowsToRemove[curRow] <- TRUE
  } else if( sum(tabulatedRowData) < numDataCols - 4 & doFullEditingForCutoffs) {
  # remove row if it has more than 4 NAs # CUSTOM
    rowsToRemove[curRow] <- TRUE
  } else if( isRowNonBiallelic(tabulatedRowDataNames) ) {
  # remove row if it has 3 or more unique alleles (eg. a row with A/A,A/C,C/C is kept, but 
  # a row with A/A,A/C,C/T is removed). We only keep biallelic data.
    rowsToRemove[curRow] <- TRUE
  }

}

report <- report[!rowsToRemove, ]

write.csv(report, paste(paste(path.expand(path), "reports", sep = "/"), "edited_report.csv", sep = "/"), row.names=FALSE)

# #######################################################################
# Function: sortDFByCol
# Sorts a data.frame by specified column(s). This was copied from: 
#   http://www.markvanderloo.eu/yaRb/2014/08/15/sort-data-frame/
# See also: 
#   https://www.r-bloggers.com/2021/02/it-has-always-been-wrong-to-call-order-on-a-data-frame/
sortDFByCol <- function(x, decreasing=FALSE, by=1, ... ) {
  f <- function(...) order(...,decreasing=decreasing)
  i <- do.call(f,x[by])
  x[i,,drop=FALSE]
}

# #######################################################################
# Function: sortDFByRow
# Sorts a data.frame by specified row.
sortDFByRow <- function(df, decreasing=FALSE, by=1) {
  dt <- as.data.frame(t(df))
  dt <- sortDFByCol(dt, decreasing=decreasing, by=by)
  as.data.frame(t(dt))
}

# #######################################################################
# Function: tallyAllelesInFactoredRow 
# input: takes as input a row of SNP data that has already been factored into totals of each genotype. 
# returns: a data.frame of one row containing the totals of each unique allele, sorted  in decreasing order.
tallyAllelesInFactoredRow <- function(factoredRow) {

  alleleSums <- data.frame("temp_col"=0) # data frame initialised with one row, and a default column "temp_col" with value of 0 (which won't be used)
  alleleSums$temp_col <- NULL # remove the temp_col, which was needed to initialize the dataframe to one row
  
  # count how many times each allele occurs    
  for ( type in names(factoredRow) ) {
    allelesInType <- strsplit(type, "/")[[1]]
    firstAllele <- allelesInType[1]
    if (firstAllele %in% names(alleleSums)) {
      alleleSums[firstAllele] <- alleleSums[firstAllele][[1]] + factoredRow[type][[1]]
    } else {
      alleleSums[firstAllele] <- factoredRow[type][[1]]
    }

    if (HAPLOID_OR_DIPLOID == 2) {
      # 1 == haploid, 2 == diploid. If it's haploid, we follow the format in the .tab file of "A/",
      # whereas if it's diploid we follow the format in the .tab file of "A/A". 
 
      secondAllele <- ifelse (length(allelesInType) > 1, allelesInType[2], allelesInType[1]) # There used to be e.g. "A/" as a shorthand for "A/A" even in the diploid data, so count it twice if it still happens to be that way
    
      if (secondAllele %in% names(alleleSums)) {
        alleleSums[secondAllele] <- alleleSums[secondAllele][[1]] + factoredRow[type][[1]]
      } else {
        alleleSums[secondAllele] <- factoredRow[type][[1]]
      }
    }
  } # end for-loop
  alleleSums <- sortDFByRow(alleleSums, decreasing = TRUE)
  alleleSums
}

# #######################################################################
message(paste0("finding snp percentage per site --- ",Sys.time()))

snpp <- report[, c(1:(startCol-1))] # the metadata columns

numColsInReport <- ncol(report)
snpp <- as.data.frame(snpp)
numRows <- nrow(snpp)

snpp$A <- integer(numRows)
snpp$C <- integer(numRows)
snpp$T <- integer(numRows)
snpp$G <- integer(numRows)
if (HAS_INDELS) {
  snpp$indel1 <- integer(numRows)
  snpp$indel2 <- integer(numRows)
  snpp$other_indels <- integer(numRows)
  snpp$indel1_val <- character(numRows)
  snpp$indel2_val <- character(numRows)
}
snpp$empty <- integer(numRows)
snpp$max <- integer(numRows)
snpp$second_max <- integer(numRows)
snpp$sum <- integer(numRows)
snpp$MAF <- numeric(numRows)

for(curRow in 1:nrow(report)) {

  rowAsMatrix <- as.matrix(report[curRow, s:numColsInReport])
  factoredGenotypes <- table(rowAsMatrix, useNA = "no")
  alleleSums <- tallyAllelesInFactoredRow(factoredGenotypes)
  alleleSums <- sortDFByRow(alleleSums, decreasing = TRUE)

  indel1 <- NA
  indel2 <- NA

  # this for-loop requires alleleSums to be sorted in decreasing order
  for ( allele in names(alleleSums)) {
    if (allele %in% c("A","C","T","G")) {
      snpp[curRow,allele] <- snpp[curRow,allele][[1]] + alleleSums[allele][[1]]
    } else if (nchar(allele) > 1 && HAS_INDELS == TRUE) {
      if (is.na(indel1)) {
        indel1 <- allele
        snpp[curRow,"indel1"] <- snpp[curRow,"indel1"][[1]] + alleleSums[allele][[1]]
        snpp[curRow,"indel1_val"] <- allele
      } else if (is.na(indel2)) {
        indel2 <- allele
        snpp[curRow,"indel2"] <- snpp[curRow,"indel2"][[1]] + alleleSums[allele][[1]]
        snpp[curRow,"indel2_val"] <- allele
      } else {
        snpp[curRow,"other_indels"] <- snpp[curRow,"other_indels"][[1]] + alleleSums[allele][[1]]
      }
    }
  }

  snpp[curRow, "empty"] <- sum(is.na(rowAsMatrix))

  if (HAS_INDELS) {  
    snpp[curRow, "max"] <- max(snpp[curRow, c("A", "C", "T", "G", "indel1", "indel2")])
    snpp[curRow, "second_max"] <- sortDFByRow(snpp[curRow, c("A", "C", "T", "G", "indel1", "indel2")], decreasing=TRUE)[2]
    snpp[curRow, "sum"] <- sum(snpp[curRow, c("A", "C", "T", "G", "indel1", "indel2")])
  } else {
    snpp[curRow, "max"] <- max(snpp[curRow, c("A", "C", "T", "G")])
    snpp[curRow, "second_max"] <- sortDFByRow(snpp[curRow, c("A", "C", "T", "G")], decreasing=TRUE)[2]
    snpp[curRow, "sum"] <- sum(snpp[curRow, c("A", "C", "T", "G")])
  }

  snpp[curRow, "MAF"] <- snpp[curRow, "second_max"] / (snpp[curRow, "second_max"] + snpp[curRow, "max"]) 
  #snpp[curRow, "chi"] <- ((snpp[curRow,"max"] - snpp[curRow,"sum"]%/%2)^2)/snpp[curRow, "sum"]%/%2  + ((snpp[curRow,"second_max"] - snpp[curRow,"sum"]%/%2)^2)/snpp[curRow, "sum"]%/%2 
}

write.csv(snpp, paste(paste(path.expand(path), "reports", sep = "/"), "percentage_snps.csv", sep = "/"), row.names=FALSE)

# #######################################################################
message(paste0("generating MAF_cutoff_report --- ",Sys.time()))

snpp <- snpp[order(-snpp$MAF), ]
snpp <- snpp[snpp$MAF >= MAF_CUTOFF,]
# CUSTOM:
#report <- report[rownames(snpp), ]
mergeBy <- c("CHROM","POS","REF")
snpp <- snpp[ , mergeBy]
message("filtering MAF report by filtered percentage_snps")
report <- merge(report, snpp, by=mergeBy, all=FALSE)

write.csv(report, paste(paste(path.expand(path), "reports", sep = "/"), "MAF_cutoff_report.csv", sep = "/"), row.names=FALSE)

# CUSTOM: 
message(paste0("snpp rownames: ",rownames(snpp)[1:100]))

# #######################################################################
message(paste0("generating site mutation percentage data --- ",Sys.time()))

fastaRef <- read.fasta(file = paste(path.expand(path), "reference/formatted_output.fasta", sep = "/"), as.string = TRUE)
mutationReport <- data.frame()

for(sector in 1:length(names(fastaRef))) {
  sectorName <- attributes(fastaRef[[sector]])$name
  sectorNameRegEx <- paste0("^",sectorName,"$")
  numRowsForSectorName <- length(grep(sectorNameRegEx, report[, "CHROM"], fixed=FALSE))
  numCharsInSequence <- nchar(fastaRef[[sector]][1])

  mutationReport[sectorName, "role"] <- attributes(fastaRef[[sector]])$Annot  
  mutationReport[sectorName, "snp"] <- numRowsForSectorName
  mutationReport[sectorName, "length"] <- numCharsInSequence
  mutationReport[sectorName, "percentage SNP"] <- numRowsForSectorName / numCharsInSequence
}

mutationReport <- mutationReport[order(-mutationReport$`percentage SNP`), ]
mutationReport$`percentage SNP` <- mutationReport$`percentage SNP` * 100

write.csv(mutationReport, paste(paste(path.expand(path), "reports", sep = "/"), "mutation_percentage.csv", sep = "/"), row.names=TRUE)


# #######################################################################
# Function: replaceSnpsWithChiCodes
# input: 
#   reportc - the report (in data.frame format) whose SNP values are to be 
#             replaced with "H", "A" or "B".
#   doExtraFiltering - default FALSE. TRUE indicates do extra row filtering,
#             which means removing rows with more than 20% NA, or if the most 
#             frequent value occurs more than 95% of the time, or if the 
#             second most frequent allele occurs less than 5% of the time.
# return: the newly updated report.
replaceSnpsWithChiCodes <- function(reportc, doExtraFiltering=FALSE) {

    # New logic:
    # remove rows with NA in MORE than 20% of sites
    # calculate occurences of each allele (eg. A/A, A/G, G/G = 3 A, 3 G)
    #   if the most frequent allele has more than 95% of occurences, remove row
    #   if the second most frequent allele has less than 5% of occurences, remove row (TODO: Or do we sum all the minor characters to get 5%, if there's more than one minor character? - I doubt this, because H,A,B encoding only allows for three characters.)
    #   set the heterozygous genotype (eg. G/T) to H, the most frequent homo type to A (eg. G/G), and the second most frequent (if it exists, eg. T/T) to B (if there's a tie in homo type frequency then set the REF/REF one to A)
    #   if there are no heterozygous genotypes, set the most frequent homo type to H and the second most to A (if there's a tie in homo type frequency then set the REF/REF one to H)
    #   set all cells that were converted to H, A or B to NA

    # TODO: The above logic is fine for biallelic sites, but what about triallelic sites? They would have three
    # possible heterozygous genotypes (REF/ALT1, REF/ALT2, ALT1/ALT2); if I set the most frequently occuring one to
    # H, what do I set the other two to? They would also have three possible homozygous genotypes (REF/REF, ALT1/ALT1, ALT2/ALT2); if I set the most frequently occuring one to A and the second most to B, what do I set the other to?

    chiRowsToRemove <- logical(nrow(reportc)) # initializes all to FALSE
    # new columns that will be added to the chi sq report
    H_Type <- character(nrow(reportc))
    A_Type <- character(nrow(reportc))
    B_Type <- character(nrow(reportc)) 
    includeB_Type <- FALSE # B type col will only be added to the report if there is a B type
    H_Total <- numeric(nrow(reportc))
    A_Total <- numeric(nrow(reportc))
    B_Total <- numeric(nrow(reportc))

    numDataCols <- ncol(reportc) - (s-1)

    for (curRow in 1:nrow(reportc)) {
      datap <- reportc[curRow, s:ncol(reportc)]

      numNAsInRow <- sum(is.na(datap)) 
      if( doExtraFiltering && numNAsInRow/numDataCols > 0.2 ) {
        # Remove the row if it has more than 20% NA's
        chiRowsToRemove[curRow] <- TRUE
      } else {
        datapTabled <- as.data.frame(table(as.matrix(datap), useNA = "no"))
        factoredData <- sortDFByCol(datapTabled, decreasing = TRUE, by=2)
        factoredRow <- factoredData[,2]
        names(factoredRow) <- factoredData[,1]
        totalAlleleCount <- sum(factoredRow)
        if (HAPLOID_OR_DIPLOID == 2) {
          # 1 == haploid, 2 == diploid. If it's haploid, we follow the format in the .tab file of "A/",
          # whereas if it's diploid we follow the format in the .tab file of "A/A".  
          totalAlleleCount <- totalAlleleCount * 2 # because, for example "A/A" has two alleles, "A/G" has two...
        }
        alleleSums <- tallyAllelesInFactoredRow(factoredRow)
        alleleSums <- sortDFByRow(alleleSums, decreasing = TRUE)
        
        if ( doExtraFiltering && alleleSums[1]/totalAlleleCount > 0.95 ) {
          # Remove the row if the most frequent allele occurs more than 95% of the time
          chiRowsToRemove[curRow] <- TRUE
        } else if ( doExtraFiltering && alleleSums[2]/totalAlleleCount < 0.05 ) {
          # Remove the row if the second most frequent allele occurs less than 5% of the time
          chiRowsToRemove[curRow] <- TRUE
        } else {
          # Process the row that we keep
          HType <- NA
          AType <- NA
          BType <- NA

          # For the code below, the possibility of length(names(alleleSums)) == 1 isn't a problem 
          # because in that case the row would be removed from the chi sq table (above) and we 
          # wouldn't be in this else clause.

          # haploid data is formatted like "G/", diploid like "G/G"
          # 1 == haploid, 2 == diploid.
          if (HAPLOID_OR_DIPLOID == 2) { 
            # Figure out which way around the heterozygous genotype should be (eg. "G/T" or "T/G") 
            # by taking the one with the most occurences (factoredRow is sorted by number of occurences)
            # (most likely only one of them will have any occurences at all)
            HTypeOpt1 <- paste0(names(alleleSums)[1],"/",names(alleleSums)[2])
            HTypeOpt2 <- paste0(names(alleleSums)[2],"/",names(alleleSums)[1])
            HTypeOpt1Frequency <- 0
            HTypeOpt2Frequency <- 0
            for (type in names(factoredRow)) {
              if (type == HTypeOpt1) {
                if (is.na(HType)) { 
                  HType <- HTypeOpt1 
                }
                HTypeOpt1Frequency <- factoredRow[type]
              } else if (type == HTypeOpt2) {
                if (is.na(HType)) { 
                  HType <- HTypeOpt2
                }
                HTypeOpt2Frequency <- factoredRow[type]
              }
            }
            # if HTypeOpt1Frequency and HTypeOpt2Frequency are tied, pick the one that starts with the REF
            if ( HTypeOpt1Frequency > 0 & HTypeOpt1Frequency == HTypeOpt2Frequency ) {
              if ( names(alleleSums)[1] == reportc[curRow, "REF"] ) {
                HType <- HTypeOpt1
              } else if ( names(alleleSums)[2] == reportc[curRow, "REF"] ) {
                HType <- HTypeOpt2
              }
            }
          } # end of dealing with HType for Diploid data 

          firstAlleleIndex <- 1
          secondAlleleIndex <- 2
          # if the first two alleles are tied in terms of frequency, choose the one that matches the REF
          # as the first one.    
          if ( alleleSums[1] == alleleSums[2] & names(alleleSums)[2] == reportc[curRow, "REF"] ) {
            firstAlleleIndex <- 2
            secondAlleleIndex <- 1
          }
          
          # if there was no heterozygous site containing the two most frequent alleles,
          # such as would be the case for Haploid data:
          if (is.na(HType)) {
            # haploid data is formatted like "G/", diploid like "G/G"
            # 1 == haploid, 2 == diploid.
            if (HAPLOID_OR_DIPLOID == 1) {
              HType <- paste0(names(alleleSums)[firstAlleleIndex],"/")
              AType <- paste0(names(alleleSums)[secondAlleleIndex],"/")
            } else {
              HType <- paste0(names(alleleSums)[firstAlleleIndex],"/",names(alleleSums)[firstAlleleIndex])
              AType <- paste0(names(alleleSums)[secondAlleleIndex],"/",names(alleleSums)[secondAlleleIndex])
            }
          } else {
            AType <- paste0(names(alleleSums)[firstAlleleIndex],"/",names(alleleSums)[firstAlleleIndex])
            BType <- paste0(names(alleleSums)[secondAlleleIndex],"/",names(alleleSums)[secondAlleleIndex])
          }
          
          reportc[curRow, s:ncol(reportc)] <- gsub(paste0("^",HType), "H", as.matrix(reportc[curRow, s:ncol(reportc)]), fixed = FALSE)
          reportc[curRow, s:ncol(reportc)] <- gsub(paste0("^",AType), "A", as.matrix(reportc[curRow, s:ncol(reportc)]), fixed = FALSE)
          H_Type[curRow] <- HType
          A_Type[curRow] <- AType
          H_Total[curRow] <- factoredRow[HType]
          A_Total[curRow] <- factoredRow[AType]

          if (!is.na(BType)) {      
            reportc[curRow, s:ncol(reportc)] <- gsub(paste0("^",BType), "B", as.matrix(reportc[curRow, s:ncol(reportc)]), fixed = FALSE)
            B_Type[curRow] <- BType
            B_Total[curRow] <- factoredRow[BType]
            if (includeB_Type == FALSE) {            
              includeB_Type <- TRUE
            }
          }

          # deal with the potential problem if the data is haploid but the data was formatted in both haploid 
          # and diploid format (e.g. both "T/" and "T/T", and HType is "T/", then the cases of "T/T" would be 
          # updated to "HT"). This shouldn't be the case with the updated pipeline, but just in case.
          if (HAPLOID_OR_DIPLOID == 1) {
            reportc[curRow, s:ncol(reportc)] <- gsub("^H[ACGTacgt]*$", "H", as.matrix(reportc[curRow, s:ncol(reportc)]), fixed = FALSE)
            reportc[curRow, s:ncol(reportc)] <- gsub("^A[ACGTacgt]*$", "A", as.matrix(reportc[curRow, s:ncol(reportc)]), fixed = FALSE)
          }
        

          # replace all elements consisting of any character(s) followed by a "/" followed by an optional character(s) with NA (i.e. all elements that haven't already been replaced with H, A or B):
          reportc[curRow, s:ncol(reportc)] <- gsub("^.*/.*$", NA, as.matrix(reportc[curRow, s:ncol(reportc)]), fixed = FALSE)
        } # finish processing row that we keep
      } 

    } # end for-loop
    
    firstDataColIndexInChiReport <- s
    # add the new columns indicating which the H, A (and B) types are to the report 
    if (includeB_Type == TRUE) {
      reportc <- cbind(reportc[ ,1:(s-1)], H_Type, A_Type, B_Type, H_Total, A_Total, B_Total, reportc[ ,s:ncol(reportc)])
      firstDataColIndexInChiReport <- s + 3
    } else {
      reportc <- cbind(reportc[ ,1:(s-1)], H_Type, A_Type, H_Total, A_Total, reportc[ ,s:ncol(reportc)])
      firstDataColIndexInChiReport <- s + 2
    }

    # remove all rows marked for removal (and return the final result to the function caller)
    reportc[!chiRowsToRemove, ]  
}

# #######################################################################
# Function: writeLinkageReport
writeLinkageReport <- function(chiInputReport, outputFileName) {

    # #######################################################################
    #message("converting the chi square report to .linkage format")

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

    reportLinkageGenotypes <- chiInputReport[ , s:ncol(chiInputReport)]
    reportLinkageGenotypes <- cbind(parent1 = c("A"), parent2 = c("H"), reportLinkageGenotypes) # add two samples to use as parents
    reportLinkageGenotypes <- t(reportLinkageGenotypes) # transpose the report (so it's columns are now rows)

    reportLinkageGenotypes[reportLinkageGenotypes=="H"] <- "1 2"
    reportLinkageGenotypes[reportLinkageGenotypes=="A"] <- "1 1"
    reportLinkageGenotypes[reportLinkageGenotypes=="B"] <- "2 2"
    reportLinkageGenotypes[is.na(reportLinkageGenotypes)] <- "0 0"
    reportLinkageGenotypes[reportLinkageGenotypes=="-"] <- "0 0" # in case NA "-" has already been substituted with "-"
    reportLinkageGenotypes[reportLinkageGenotypes==""] <- "0 0" # in case NA has been substituted with ""
    reportLinkageGenotypes[reportLinkageGenotypes==" "] <- "0 0" # in case NA has been substituted with " "
    reportLinkageGenotypes[reportLinkageGenotypes=="N/A"] <- "0 0" 
    reportLinkageGenotypes[reportLinkageGenotypes=="na"] <- "0 0" 
    reportLinkageGenotypes[reportLinkageGenotypes=="n/a"] <- "0 0" 

    reportLinkage <- cbind(family = c("chi"), id = c(paste0("S",(1:nrow(reportLinkageGenotypes))-2)), fatherId = c("P1"), motherId = c("P2"), gender = c(0), affectionStatus = c(0), reportLinkageGenotypes)
    reportLinkage[1,2:5] <- c("P1","0","0","1") # change id from S-1 to P1, no parents, male
    reportLinkage[2,2:5] <- c("P2","0","0","2") # change id from S0 to P2, no parents, female

    write.table(reportLinkage, file= paste(paste(path.expand(path), "reports", sep = "/"), outputFileName, sep = "/"), append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

} # end of function writeLinkageReport


# #######################################################################
if (HAS_INDELS) {
  message(paste0("generating indels report from MAF report --- ",Sys.time()))

  # keep only those rows whose REF value has more than one character in it (i.e. it's an indel):
  indelReport <- report[ (nchar(report$REF) > 1), ]
  indelReport <- replaceSnpsWithChiCodes(indelReport, doExtraFiltering=FALSE)

  write.csv(indelReport, paste(paste(path.expand(path), "reports", sep = "/"), "MAF_cutoff_report_indels.csv", sep = "/"), row.names=FALSE, na="-")
}

# #######################################################################
# The chi sq report, and the .linkage file derived from it, will only be 
# generated if GENERATE_CHI_SQ_REPORT == 1
if (GENERATE_CHI_SQ_REPORT == 1) {
    message(paste0("replacing alleles with characters for chi square test, with extra filtering --- ",Sys.time()))
    reportcfiltered <- replaceSnpsWithChiCodes(report, doExtraFiltering=TRUE)
    write.csv(reportcfiltered, paste(paste(path.expand(path), "reports", sep = "/"), "MAF_cutoff_report_chi_extraFiltered.csv", sep = "/"), row.names=FALSE, na="-")    
    writeLinkageReport(reportcfiltered, "MAF_cutoff_report_chi_extraFiltered.linkage")

    message(paste0("replacing alleles with characters for chi square test --- ",Sys.time()))
    reportc <- replaceSnpsWithChiCodes(report, doExtraFiltering=FALSE)
    write.csv(reportc, paste(paste(path.expand(path), "reports", sep = "/"), "MAF_cutoff_report_chi.csv", sep = "/"), row.names=FALSE, na="-")
    writeLinkageReport(reportcfiltered, "MAF_cutoff_report_chi.linkage")

} # end of the optional block to generate the chi sq report, and the .linkage file derived from it.


# #######################################################################
# The probability report will only be generated if GENERATE_PROBABILITY_REPORT == 1
if (GENERATE_PROBABILITY_REPORT == 1) {

    message(paste0("generate probability values --- ",Sys.time()))

    reportd <- report 

    startingCol <- ifelse(("COMBINED" %in% colnames(reportd)), s-1, s) # this makes sure we include the COMBINED col in our for loop if it exists

    for(x in startingCol:ncol(reportd)) {
    
      print(paste0("------------col ",x)) 
      if(colnames(reportd)[x] == "COMBINED") {
        for(cutf in list.files(paste0(path, "/outputTemp/pooled"))) {
          if(grepl("cutoff", cutf, fixed=TRUE)) {
            fil <- paste0(path, "/outputTemp/pooled/", cutf)
          }
        }
      } else {
        fil <- paste0(path, "/outputTemp/single/", substr(colnames(reportd)[x], 0, nchar(colnames(reportd)[x]) - 4), "_cutoff")
      }
      # For documentation on the VariantAnnotation packaged used for scanVcf and related operations, see:
      #     http://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html
      # For specification of the vcf file format, see:
      #     https://github.com/samtools/hts-specs (has canonical specifications of VCF format)
      #     http://www.htslib.org/doc/vcf.html (seems to briefly document the VCF format)
      # Of particular relevance to the code below, see the descriptions of GT and GL in GENO, section 1.4.2, 
      # page 5 of:    
      #     http://samtools.github.io/hts-specs/VCFv4.1.pdf
      #     In the formula given in the documentation for GL: F(j/k) = (k*(k+1)/2)+j,  j and k are integers
      #     representing REF or ALT or ALT2 (ALT only for triallelic, which we don't have). 
      #     (0==REF, 1==ALT, 2==ALT2). So if the genotype is REF/ALT, use j=0 and k=1. Then, 
      #     the result of the formula will be either 0, 1, or 2 (or also 3, 4, 5 if triallelic data), which we
      #     then use to choose the right index for the GL probability. However, since indices in R start 
      #     at 1, we will treat 0 as GL position 1, 1 as GL position 2, 2 as GL position 3. 
      ap <- scanVcf(fil) # import vcf file (see page 35 of http://bioconductor.org/packages/release/bioc/manuals/VariantAnnotation/man/VariantAnnotation.pdf)
      rr <- ap$`*:*-*`$rowRanges
      for(y in 1:nrow(reportd)) {
      
        ind <- NA
        if(!is.na(reportd[y, x])) {
          # Haploid data, and the character matches the REF.
          # In the case where the format is "A/" or "A" instead of "A/A", it is because the data is haploid.
          # If the character matches the REF, set the likelihood value to 1.
          if(nchar(reportd[y,x]) <= 2 & strsplit(reportd[y,x], "")[[1]][1] == reportd[y, "REF"] ) {
            reportd[y,x] <- 1
          } else if(nchar(reportd[y,x]) <= 2 & strsplit(reportd[y,x], "")[[1]][1] != reportd[y, "REF"] ) {
          # Haploid data, and the character doesn't match the REF.
            # see comment below on the equivalent assignment of ind in the next else-if clause (in the
            # three character case):
            ind <- intersect(grep(paste0(":", reportd[y, "POS"], "_"), names(rr), value=FALSE, fixed=TRUE), grep(reportd[y, "CHROM"], names(rr), value=FALSE, fixed=TRUE))
            
            if(!is.null(ind) && length(ind) > 0 && !is.na(ind)) {
              # Check the Genotype (GT) to determine which GL position to use.
              # Since this is haploid data, and we know already that it doesn't
              # match the REF, we expect the GT to be 1 (==ALT allele).
              if( as.list(ap$`*:*-*`$GENO$GT[ind])[[1]][1] == "1" ) { # "1" means GL position 2
                reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][2]
              } else { # Just in case the GT is "0" (==REF allele), it means GL position 1
                reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][1]
              }
            } else { # if ind is NA or null or has no value
              reportd[y,x] <- NA
            }
          } else if(nchar(reportd[y,x]) == 3) {
          # three characters, eg. "A/A" or "A/T".
            # Set "ind" to the list of the indices of those row range names which contain both this row's "CHROM"
            # value and ":<this row's POS value>_". For our purposes, we only care about the first one in the
            # list if there are more than one. It'll be used to look up the Genotype Likelihood values for the
            # current cell.
            # See the descriptions of GT and GL in GENO, section 1.4.2, page 5 of:    
            #     http://samtools.github.io/hts-specs/VCFv4.1.pdf 
            # as well as my previous comments in the code regarding that document, above this for-loop.
            # (Using the fixed=TRUE parameter in grep prevents it from evaluating the first string as a regular
            # expression - important, because it would cause problems for special characters otherwise).
            # (Using the value=FALSE parameter in grep makes it return the integer indices of matching elements,
            # rather than the elements themselves. It's the default, but it doesn't hurt to be specific). 
            ind <- intersect(grep(paste0(":", reportd[y, "POS"], "_"), names(rr), value=FALSE, fixed=TRUE), grep(reportd[y, "CHROM"], names(rr), value=FALSE, fixed=TRUE))

            if(!is.null(ind) && length(ind) > 0 && !is.na(ind)) {
              # if the 1st and 3rd match each other:
              if(strsplit(reportd[y,x], "")[[1]][1] == strsplit(reportd[y,x], "")[[1]][3] ) {
                # if they match the REF:
                if(strsplit(reportd[y,x], "")[[1]][1] == reportd[y, "REF"] ) {
                  # (GL position 1)
                  reportd[y,x] <- 1
                  # TODO: Not yet sure if I just want to set it to 1 automatically here (as above), but if I 
                  # don't, use this code instead:
                  #reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][1]
                } else {
                # if they DON'T match the REF.
                  # Check the Genotype (GT) to determine which GL position to use
                  if( as.list(ap$`*:*-*`$GENO$GT[ind])[[1]][1] == "1/1" ) { # "1/1" means GL position 3
                    reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][3]
                  } else { # "2/2" means GL position 6
                    reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][6]
                  }
                }
              } else {
              # if the 1st and 3rd DON'T match each other.
                # if the first character matches the REF:
                if(strsplit(reportd[y,x], "")[[1]][1] == reportd[y, "REF"] ) {
                  # Check the Genotype (GT) todetermine which GL position to use
                  if( as.list(ap$`*:*-*`$GENO$GT[ind])[[1]][1] == "0/1" ) { # "0/1" means GL position 2
                    reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][2]
                  } else { # "0/2" means GL position 4
                    reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][4]
                  }
                } else {
                # if the first character DOESN'T match match the REF:
                  # it must be GL position 5 (only possible for triallelic sites)
                  reportd[y, x] <- 10 ^ as.list(ap$`*:*-*`$GENO$GL[ind])[[1]][5]
                }
              }
            } else { # if ind is NA or null or has no value

              # if the 1st and 3rd character match each other and match the REF
              if( strsplit(reportd[y,x], "")[[1]][1] == reportd[y, "REF"] & strsplit(reportd[y,x], "")[[1]][1] == strsplit(reportd[y,x], "")[[1]][3]) {
                reportd[y,x] <- 1
              } else {
                reportd[y, x] <- NA
              }
            }
            # end of 3 character cases
          } else {
            reportd[y, x] <- NA
          }

          # Old Comments:
          # TODO: If we raise 10 ^ (much greater than -320) I think the number is too large, resulting in 0.


        } # end of if(!is.na(reportd[y, x]))
      } # end for loop (iterating down the rows)
    } # end for loop (iterating across the colums)

    write.csv(reportd, paste(paste(path.expand(path), "reports", sep = "/"), "probability.csv", sep = "/"), row.names=FALSE)

} # end of the optional block to generate the probability report.

message(paste0("report_gen part 2 complete --- ",Sys.time()))
