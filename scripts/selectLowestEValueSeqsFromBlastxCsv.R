options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)

#write("Running selectSeqsFromFasta.R.", stdout())
write("Running selectLowestEValueSeqsFromBlastxCsv.R.", stdout())

# The path to the folder containing the files. All file paths should be relative to this path, 
# even if they are not in the folder itself.
path <- "."


# Csv file, Col1: sequence ids, Col2: sequence lengths:
#seqFile <- "Hemlock30-blastx-tabular-100lines.csv"
seqFile <- "Hemlock30_longest-90,730-blastx-6outfmt.MASTER.csv"

#outputFile <- "Hemlock30-blastx-tabular-100lines-topPicks.csv"
outputFile <- "Hemlock30_longest-90,730-blastx.MASTER-topPicks.csv"


##############################################################
# Excecution Code. There should be no need to edit anything 
# below this line.



write(paste0("Input seqFile: ",seqFile), stdout())

write(paste0("Reading the seqFile."), stdout())
sequences <- read.csv(paste(path.expand(path),seqFile,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.


#seqIds <- strsplit(sequences[,1], "_seq", fixed=TRUE)
#seqIds <- as.data.frame(do.call(rbind, seqIds))[,1]
#
#sequences <- cbind(seqIds, sequences)


uniqSeqIds <- unique(sequences[,1])
len <- length(uniqSeqIds)

selected <- data.frame(character(len),character(len),character(len),numeric(len),numeric(len),character(len))
names(selected) <- names(sequences)

write(paste0("Finding the row per sequence with the lowest e-value."), stdout())

for (i in 1:len) {
  group <- sequences[sequences[,1]==uniqSeqIds[i],]
  # selecting the top pick based on smallest e-value, assuming the 5th column contains the e-values
  selected[i,] <- group[which.min(group[,5]), ]  
}

write.csv(selected, paste(path.expand(path),outputFile,sep="/"), row.names=FALSE)

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())


