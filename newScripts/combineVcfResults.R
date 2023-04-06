options(stringsAsFactors = FALSE, warn = 1)

MAF_CUTOFF <- 0.05

path <- "."

mafFilename <- "../pipe1_MAF_reports/pipe1_MAF_cutoff_report_filtered.csv"

samplesDir <- "../pipe1csv-V1"
sampleFilesPattern <- "_filtered_cutoff\\.csv$"


# Pseudocode:
#
# read MAF_cutoff_report.csv
# set MAF report "CHROM","POS","REF" as mafMaster
# for each vcf report
# 	read report
#	inner-join report with mafMaster to filter out rows in report which don't correspond to mafMaster
# 	replace all 0/1 etc with REF/ALT etc
# 	if first report
# 		set these columns of the report as masterForDeepBSA: c("CHROM","POS","REF","ALT","RO","AO")
# 		set the report as masterGeneral
# 	else 
# 		merge with masterForDeepBSA (outer-join) by "CHROM","POS","REF","ALT", adding cols "sample_RO","sample_AO"
# 		merge with masterGeneral (outer-join) by "CHROM","POS","REF", adding cols "sample_ALT", "sample_GT", "sample_DP", "sample_RO", "sample_AO"
# done
#
# To potentially do later:
# split the master report into separate reports according to use:
# 	genotype report: "CHROM","POS","REF", "sample_GT"...
#	depth stats report: "CHROM","POS","REF", "sample_DP", "sample_RO", "sample_AO"...
#
# If operating this script on different pipelines: finally, concat final reports from all three pipelines


# #######################################################################
message(paste0("Running combineVcfResults.R --- ",Sys.time()))

message(paste0("Main directory is ",path.expand(path)," --- ",Sys.time()))
message(paste0("Maf report is ",paste(path.expand(path), mafFilename, sep="/")," --- ",Sys.time()))
message(paste0("Samples director is ",paste(path.expand(path), samplesDir, sep="/")," --- ",Sys.time()))
message(paste0("sampleFilesPattern is ",sampleFilesPattern," --- ",Sys.time()))

maf <- read.csv(paste(path.expand(path), mafFilename, sep="/"), header=TRUE, check.names=FALSE)

intersectColsMin <- c("CHROM","POS","REF")
intersectColsMax <- c("CHROM","POS","REF","ALT")
mafMaster <- maf[,intersectColsMin]
masterForDeepBSA <- NULL
masterGeneral <- NULL

isFirst <- TRUE

for(fn in list.files(paste(path.expand(path), samplesDir, sep="/"), pattern=sampleFilesPattern)) {
  message(paste0("Processing ",fn," --- ",Sys.time()))
  sample <- read.csv(paste(path.expand(path), samplesDir, fn, sep="/"), header=TRUE, check.names=FALSE)
  # This is important for the following merges to work as intended: filter out rows in sample which don't correspond to mafMaster 
  sample <- merge(mafMaster, sample, by=intersectColsMin, all=FALSE)
  
  genotypeSplit <- strsplit(sample$GT, split="/", fixed=TRUE)
  genotypeSplit <- as.data.frame(genotypeSplit)
  names(genotypeSplit) <- NULL
  genotypeSplit <- t(genotypeSplit) # transpose to put it back into two vectors, one containing the first alleles and the other containing the second alleles
  
  # replace all 0's with REF and replace all 1's with ALT
  allele1 <- ifelse(genotypeSplit[,1]=="0", sample$REF, sample$ALT)
  allele2 <- ifelse(genotypeSplit[,2]=="0", sample$REF, sample$ALT)  
  sample$GT <- paste(allele1,allele2,sep="/")
 
  sampleForDeepBSA <- sample[, c("CHROM","POS","REF","ALT","RO","AO")]

  # remove the common ending portion of the filename to get the sampleName
  sampleName <- sub(sampleFilesPattern, "", fn)
  
  colnames(sample)[which(names(sample) == "ALT")] <- paste0(sampleName,"_","ALT")
  colnames(sample)[which(names(sample) == "GT")] <- paste0(sampleName,"_","GT")
  colnames(sample)[which(names(sample) == "DP")] <- paste0(sampleName,"_","DP")
  colnames(sample)[which(names(sample) == "RO")] <- paste0(sampleName,"_","RO")
  colnames(sample)[which(names(sample) == "AO")] <- paste0(sampleName,"_","AO")
  
  colnames(sampleForDeepBSA)[which(names(sampleForDeepBSA) == "RO")] <- paste0(sampleName,"_","RO")
  colnames(sampleForDeepBSA)[which(names(sampleForDeepBSA) == "AO")] <- paste0(sampleName,"_","AO")

  if (isFirst) {
    isFirst <- FALSE
    masterForDeepBSA <- sampleForDeepBSA
    masterGeneral <- sample
  } else {
    masterForDeepBSA <- merge(masterForDeepBSA, sampleForDeepBSA, by=intersectColsMax, all=TRUE)
    masterGeneral <- merge(masterGeneral, sample, by=intersectColsMin, all=TRUE)
  }
  
}

write.csv(masterForDeepBSA, paste(path.expand(path), "MAF_cutoff_depth_forDeepBSA.csv", sep = "/"), row.names=FALSE)
write.csv(masterGeneral, paste(path.expand(path), "MAF_cutoff_all_data.csv", sep = "/"), row.names=FALSE)


message(paste0("DONE --- ",Sys.time()))


