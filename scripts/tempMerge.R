args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
reportName <- args[2]
report <- data.frame()
single <- FALSE

message("temp Merge.")

#path <- "~/Desktop/SNPpipeline"
options(stringsAsFactors = FALSE, warn = 1)

message(path)





#########

twoColReport <- read.csv(paste(path.expand(path),"reports-BAD",reportName,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

sFilesReport <- read.csv(paste(path.expand(path),"reports-justSfiles",reportName,sep="/"),header=TRUE, check.names=FALSE)

report <- merge(twoColReport, sFilesReport, by = c("CHROM", "POS", "REF"), all = TRUE)

#########


write.csv(report, paste(path.expand(path), "reports", reportName, sep = "/"), row.names=FALSE)

message("temp merge completed")
