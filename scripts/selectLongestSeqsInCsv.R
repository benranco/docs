options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)

write("Running selectSeqsFromFasta.R.", stdout())

# The path to the folder containing the files. All file paths should be relative to this path, 
# even if they are not in the folder itself.
path <- "."


# Csv file, Col1: sequence ids, Col2: sequence lengths:
seqFile <- "Ht5313+5396-Trinity85,183-length.csv"

outputFile <- "Ht5313+5396-Trinity85,183-longest.csv"



##############################################################
# Excecution Code. There should be no need to edit anything 
# below this line.



write(paste0("Input seqFile: ",seqFile), stdout())

write(paste0("Reading the seqFile."), stdout())
sequences <- read.csv(paste(path.expand(path),seqFile,sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.


genes <- strsplit(sequences[,1], "_seq", fixed=TRUE)
genes <- as.data.frame(do.call(rbind, genes))[,1]

sequences <- cbind(genes, sequences)

uniqGenes <- unique(genes)
len <- length(uniqGenes)

selected <- data.frame(character(len),character(len),numeric(len))
names(selected) <- c("Gene","Seq","Length")

write(paste0("Finding longest sequences."), stdout())

for (i in 1:len) {
  group <- sequences[sequences[,1]==uniqGenes[i],]
  selected[i,] <- group[which.max(group[,3]), ]
}

write.csv(selected, paste(path.expand(path),outputFile,sep="/"), row.names=FALSE)

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())


