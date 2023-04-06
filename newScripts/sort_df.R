options(stringsAsFactors = FALSE, warn = 1)
path <- "."

old <- read.csv(paste(path.expand(path),"",paste0("MAF_cutoff_report-oldMerge.csv"),sep="/"),header=TRUE, check.names=FALSE)
new <- read.csv(paste(path.expand(path),"",paste0("MAF_cutoff_report-oldMerge.csv"),sep="/"),header=TRUE, check.names=FALSE)

# see examples in the R documentation for the order function.
old <- old[ do.call(order, old[,1:3]), ]
new <- new[ do.call(order, new[,1:3]), ]

write.csv(old, "old.csv", row.names=FALSE)
write.csv(new, "new.csv", row.names=FALSE)

# Also, see this for an example of how to merge two dataframes while preserving the order of the original: 
#   https://stackoverflow.com/questions/17878048/merge-two-data-frames-while-keeping-the-original-row-order
#
# Example, merging two dataframes called stats and correlData:
#
## used to preserve the order of stats after doing the merge
# stats$index  <- 1:nrow(stats)
#
# message(paste0("Merging stats data and correlation data --- ",Sys.time()))
# newstats <- merge(stats, correlData, by=c("CHROM","POS","REF","ALT"), all=TRUE)
#
## set the order back to the original order of stats
# newstats <- newstats[order(newstats$index), ]
# newstats <- newstats[, -which(names(newstats) %in% c("index"))]


