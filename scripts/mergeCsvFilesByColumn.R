options(stringsAsFactors = FALSE, warn = 1)

write("Running combineCsvFilesByColumn.R.", stdout())


# ########################################################
# Input parameters:

# Use "." if the files are in the same folder as this script, otherwise enter the folder path 
# in which to find the files. E.g. on Windows: "C:\path\to\folder"
path <- "."

# List the names of the files you want to combine, in the order you want to combine them.
# Import: the file with the master list of id's needs to be first.
inputFileNames <- c("file1.csv", "file2.csv", "file3.csv" )

# The character used to delimit fields in the input files (e.g. "," for comma, "\t" for tab) 
inputFileFieldDelimiter <- "\t"

# The index or name of the column to merge by. If a name is used, the column name must be 
# identical across all files.
mergeCol <- 1

# The name of the output file.
outputFileName <- "merged_files.csv"

##########
# These parameters are for optional use:

# Indicates whether to create a new gene name column from extracting the gene names out of 
# pre-existing columns. TRUE if you want to, FALSE if you don't. If TRUE, it will also sort
# the output file by gene name in decreasing order.
createNewGeneNameColumn <- TRUE

# if the above createNewGeneNameColumn parameter is TRUE, these two columns are used to 
# extract the gene name. 
firstColumnToLookForGeneName <- "sprot_Top_BLASTX_hit"
secondColumnToLookForGeneName <- "sprot_Top_BLASTP_hit"

# the character you'd like to use for the cases where there is no gene name.
noGeneNameIndicator <- "."

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

# optionally, create a new gene column based on data in pre-existing columns 
if (createNewGeneNameColumn == TRUE) {
  gene <- character(nrow(master))

  # extract gene names from specific columns and save them in a new column
  for (i in 1:nrow(master)) {
    gname <- master[i, firstColumnToLookForGeneName]
    if(is.na(gname) || nchar(gname) == 0) {
      gname <- master[i, secondColumnToLookForGeneName]
    }

    if(!is.na(gname) && nchar(gname) > 0) {
      gname <- strsplit(gname,"Full=")[[1]][2]
      gname <- strsplit(gname,";")[[1]][1]
      if(is.na(gname) || nchar(gname) == 0) {
        gname <- noGeneNameIndicator
      }
    }
    else {
      gname <- noGeneNameIndicator
    }

    gene[i] <- gname
  } # end for-loop

  # sort it by the new gene column, in decreasing order so empty values are at the bottom
  master <- cbind(gene, master)
  master <- master[order(master$gene, decreasing=TRUE), ]

}

write.table(master, file= paste(path.expand(path),outputFileName,sep="/"), append=FALSE, quote=FALSE, sep=inputFileFieldDelimiter, row.names=FALSE, col.names=TRUE)

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())



