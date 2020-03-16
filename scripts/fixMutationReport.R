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




path <- "/home/benrancourt/Desktop/junjun/reports-42-samples-corrected-without1000rowsfromP35"

report <- read.csv(paste(path.expand(path),".","MAF_cutoff_report.csv",sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.


# #######################################################################
message("generating site mutation percentage data")

fastaRef <- read.fasta(file = paste(path.expand(path), "reference/formatted_output.fasta", sep = "/"), as.string = TRUE)
mutationReport <- data.frame()

for(sector in 1:length(names(fastaRef)))
{
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

write.csv(mutationReport, paste(paste(path.expand(path), ".", sep = "/"), "mutation_percentage-NEW.csv", sep = "/"), row.names=TRUE)


