
MAF_CUTOFF <- 0.05
path <- "."

snpp <- read.csv(paste(path.expand(path),"","percentage_snps-10000.csv",sep="/"),header=TRUE, check.names=FALSE)

# #######################################################################
message(paste0(" --- ",Sys.time()))

snpp <- snpp[order(-snpp$MAF), ]
snpp <- snpp[snpp$MAF >= MAF_CUTOFF,]

write.csv(snpp, paste(paste(path.expand(path), "", sep = "/"), "percentage_snps_filtered.csv", sep = "/"), row.names=FALSE)
message(paste0("DONE --- ",Sys.time()))

