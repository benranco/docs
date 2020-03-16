options(stringsAsFactors = FALSE, warn = 1)

write("Running combineFilledReportsByColumn.R.", stdout())


# ########################################################
# Input parameters:

path <- "/work"

inputReportNames <- c("filled_report-01-backup.csv", 
                      "filled_report-02-backup.csv", 
                      "filled_report-03-backup.csv", 
                      "filled_report-04-backup.csv", 
                      "filled_report-05-backup.csv", 
                      "filled_report-06-backup.csv", 
                      "filled_report-07-backup.csv", 
                      "filled_report-08-backup.csv", 
                      "filled_report-09-backup.csv", 
                      "filled_report-10-backup.csv" )

mergeColNames <- c("CHROM","POS","REF")

outputFileName <- "filled_report-combined.csv"

# ########################################################

# using check.names=FALSE for these tables in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.
write(paste0("Loading report 1"), stdout())
master <- read.csv(paste(path.expand(path),inputReportNames[1],sep="/"),header=TRUE, check.names=FALSE)

for (i in 2:length(inputReportNames)) {
  write(paste0("Merging report ",i), stdout())
  nextReport <- read.csv(paste(path.expand(path),inputReportNames[i],sep="/"),header=TRUE, check.names=FALSE)
  master <- merge(x=master, y=nextReport, all=TRUE, by=mergeColNames)
}

write.table(master, file= paste(path.expand(path),outputFileName,sep="/"), append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())



