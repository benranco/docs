options(stringsAsFactors = FALSE, warn = 1)

write("Running combineCsvFilesByColumn.R.", stdout())


# ########################################################
# Input parameters:

# Use "." if the files are in the same folder as this script, otherwise enter the folder path 
# in which to find the files. E.g. on Windows: "C:\path\to\folder"
path <- "/home/benrancourt/Desktop/junjun/33sampleRNA/reports-pipe02-done/again"

# List the names of the files you want to combine, in the order you want to combine them.
# Import: the file with the master list of id's needs to be first.
inputFileNames <- c("MAF_RowsWithRefNotInSC.csv", "depth-rows_with_SC_ro_NOTzero.csv" )

# The character used to delimit fields in the input files (e.g. "," for comma, "\t" for tab) 
inputFileFieldDelimiter <- ","

# The index or name of the column to merge by. If a name is used, the column name must be 
# identical across all files.
mergeCol <- c("CHROM","POS","REF")

# The name of the output file.
#outputFileName <- "merged_files_RSdata_roNOTzero.csv"
outputFileName <- "merged_files_SCdata_roNOTzero.csv"

##########
# These parameters are for optional use:


# ########################################################
# No need to edit below this line.

write(paste0("Loading report 1"), stdout())
master <- read.csv(paste(path.expand(path),inputFileNames[1],sep="/"),header=TRUE, sep=inputFileFieldDelimiter, check.names=FALSE)
# using check.names=FALSE for these tables in case the column names have dashes (-) in them. This will prevent them from being converted to periods. (However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them).


# merge files
for (i in 2:length(inputFileNames)) {
  write(paste0("Merging report ",i), stdout())
  nextReport <- read.csv(paste(path.expand(path),inputFileNames[i],sep="/"),header=TRUE, sep=inputFileFieldDelimiter, check.names=FALSE)
  master <- merge(x=master, y=nextReport, all=FALSE, by=mergeCol)
}

write.csv(master, paste(path.expand(path), outputFileName, sep = "/"), row.names=FALSE)

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())



