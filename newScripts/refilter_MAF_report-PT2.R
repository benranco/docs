options(stringsAsFactors = FALSE, warn = 1)

MAF_CUTOFF <- 0.05
path <- "."

reportPrefix <- "pipeAll_"

snpp <- read.csv(paste(path.expand(path),"",paste0(reportPrefix,"percentage_snps.csv"),sep="/"),header=TRUE, check.names=FALSE)


# #######################################################################
message(paste0("filtering ",reportPrefix,"percentage_snps by MAF_CUTOFF >= ",MAF_CUTOFF," --- ",Sys.time()))

snpp <- snpp[order(-snpp$MAF), ]
snpp <- snpp[snpp$MAF >= MAF_CUTOFF,]

write.csv(snpp, paste(paste(path.expand(path), "", sep = "/"), paste0(reportPrefix,"percentage_snps_filtered.csv"), sep = "/"), row.names=FALSE)

# #######################################################################
message(paste0("generating ",reportPrefix,"MAF_cutoff_report_filtered --- ",Sys.time()))

maf <- read.csv(paste(path.expand(path),"",paste0(reportPrefix,"MAF_cutoff_report.csv.noIndels.max4NAs"),sep="/"),header=TRUE, check.names=FALSE)

mergeBy <- c("CHROM","POS","REF")
snpp <- snpp[ , mergeBy]

message("filtering MAF report by filtered percentage_snps")
maf <- merge(maf, snpp, by=mergeBy, all=FALSE)
write.csv(maf, paste(paste(path.expand(path), "", sep = "/"), paste0(reportPrefix,"MAF_cutoff_report_filtered.csv"), sep = "/"), row.names=FALSE)


message(paste0("DONE --- ",Sys.time()))


