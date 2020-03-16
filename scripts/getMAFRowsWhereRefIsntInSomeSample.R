args <- commandArgs(trailingOnly = TRUE)


message("Running getMAFRowsWhereRefIsNotInSomeSample.R.")

#path <- "~/Desktop/SNPpipeline"
options(stringsAsFactors = FALSE, warn = 1)



path <- "/home/benrancourt/Desktop/junjun/33sampleRNA/reports-pipe02-done/again"
mafReportName <- "MAF_cutoff_report.csv"
outputReportName <- "MAF_RowsWithRefNotInSC.csv"


#########

mafReport <- read.csv(paste(path.expand(path),mafReportName,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

flaggedRows <- logical(nrow(mafReport)) # initiates all to FALSE


for (rowNum in 1:nrow(mafReport)) {

  ref <- mafReport[rowNum,"REF"]

  rs <- mafReport[rowNum,"Bigfile_RS.tab"]
  sc <- mafReport[rowNum,"Bigfile_SC.tab"]

#  if (!is.na(rs)) {
#    rs <- strsplit(rs,"/")[[1]]
#    if ( !(ref == rs[1] || ref == rs[2]) ) {
#      flaggedRows[rowNum] <- TRUE
#    }
#  }

  if (!is.na(sc)) {
    sc <- strsplit(sc,"/")[[1]]
    if ( !(ref == sc[1] || ref == sc[2]) ) {
      flaggedRows[rowNum] <- TRUE
    }
  }

} # end for-loop

report <- mafReport[flaggedRows, ]


#########


write.csv(report, paste(path.expand(path), outputReportName, sep = "/"), row.names=FALSE)

message("...done.")
