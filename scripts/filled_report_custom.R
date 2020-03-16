args <- commandArgs(trailingOnly = TRUE)
path <- args[1]

filledReport <- read.csv(paste(path.expand(path),"reports","filled_report.csv",sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

single1 <- "report_p35_filled.Rds"

# #######################################################################
message("reading and merging split data")
#for(single1 in list.files(paste0(path, "/reporttemp")))
#{
  if(grepl("filled.Rds", single1, fixed=TRUE))
  {
    #message(single1)
    sr <- readRDS(paste0(path, "/reporttemp/", single1))
    if(nrow(report) == 0 && ncol(report) == 0)
    {
      report <- sr
    }
    else
    {
      report <- rbind(report, sr)
    }
  }
#}

report <- rbind(filledReport, report)

write.csv(report, paste(paste(path.expand(path), "reports", sep = "/"), "filled_report_NEW.csv", sep = "/"), row.names=FALSE)

