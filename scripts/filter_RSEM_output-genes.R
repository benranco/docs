options(stringsAsFactors = FALSE, warn = 1)

write("Running combineCsvFilesByColumn.R.", stdout())

# columns by which to merge tables. Set to c(1) for workingwith RSEM.genes.results and c(1,2) for RSEM.isoforms.results
mergeCols <- c(1)

fileNameGeneric <- "RSEM.genes.results.small10000"
#fileNameGeneric <- "RSEM.isoforms.results"
inPath <-  "."
outPath <- "./RSEM_filtered"


# The character used to delimit fields in the input files (e.g. "," for comma, "\t" for tab) 
inputFileFieldDelimiter <- "\t"




# ########################################################
# No need to edit below this line.


inFileList <- system2("find", args=c(inPath, paste0(" -mindepth 2 -maxdepth 2 -name ",fileNameGeneric," -type f")), stdout=TRUE)

master <- NULL
chosenGeneIds <- NULL

# merge files
for (i in 1:length(inFileList)) {
#for (i in 1:1) {
  write(paste0("Loading ",inFileList[i]), stdout())
  nextReport <- read.csv(inFileList[i],header=TRUE, sep=inputFileFieldDelimiter, check.names=FALSE) # using check.names=FALSE for these tables in case the column names have dashes (-) in them. This will prevent them from being converted to periods. (However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them).
  
  nextReportFiltered <- nextReport[ nextReport$TPM >= 1, ]

  splitFilePath <- strsplit(inFileList[i], "/", fixed=TRUE)[[1]]
  prefix <- splitFilePath[length(splitFilePath)-1]
  
  write.table(nextReportFiltered, file= paste0(outPath,"/",prefix,".",fileNameGeneric,".filtered.tab"), append=FALSE, quote=FALSE, sep=inputFileFieldDelimiter, row.names=FALSE, col.names=TRUE, na="")
  
  idColNames <- names(nextReport)[mergeCols]
  newNames <- paste(prefix,names(nextReport),sep=".")
  newNames[mergeCols] <- idColNames
  names(nextReport) <- newNames
  
  if (is.null(chosenGeneIds)) {
    chosenGeneIds <- nextReportFiltered[,idColNames[1]]
  }
  else {
    chosenGeneIds <- c(chosenGeneIds, nextReportFiltered[,idColNames[1]])
  }
  
  if (is.null(master)) {
    master <- nextReport
  }
  else {
    master <- merge(x=master, y=nextReport, all=TRUE, by=idColNames)
  }
  
  
}

uniqChosenGeneIds <- unique(chosenGeneIds)

#master <- master[ master[,1] %in% uniqChosenGeneIds, ]
master <- master[ match(uniqChosenGeneIds, master[,1]), ]  # assumes only one occurence of each gene id in master 


# add three more column as 'total_expected_count', 'total_TPM', and 'total_FPKM'

expCount <- master[, grep("expected_count", names(master), fixed=TRUE)]
tpm <- master[, grep("TPM", names(master), fixed=TRUE)]
fpkm <- master[, grep("FPKM", names(master), fixed=TRUE)]

total_expected_count <- rowSums(expCount)
total_TPM <- rowSums(tpm)
total_FPKM <- rowSums(fpkm)

mergeColNames <- names(nextReport)[mergeCols]
master <- cbind(master[, mergeCols], total_expected_count, total_TPM, total_FPKM, master[, -mergeCols])
names(master)[mergeCols] <- mergeColNames  # this is only necessary if there was only one mergeCol in the above cbind 

write.table(master, file= paste0(outPath,"/","ALLsamples.",fileNameGeneric,".filtered.tab"), append=FALSE, quote=FALSE, sep=inputFileFieldDelimiter, row.names=FALSE, col.names=TRUE, na="")

write.table(master[,1], file= paste0(outPath,"/","ALLsamples.",fileNameGeneric,".selected.tab"), append=FALSE, quote=FALSE, sep=inputFileFieldDelimiter, row.names=FALSE, col.names=FALSE, na="")

write(uniqChosenGeneIds, file= paste0(outPath,"/","ALLsamples.",fileNameGeneric,".chosen.tab"), append=FALSE, ncolumns=1, sep=inputFileFieldDelimiter)


