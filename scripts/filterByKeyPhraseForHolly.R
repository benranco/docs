options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)


write("Running filterByKeyPhraseForHolly.R.", stdout())

##########################################################
# What this script does:
#
# - accept a list of key phrases as input, and a data input file in .csv format
# - keep rows containing at least one of the key phrase and discard the other rows
#     - if the key phrase is "cell wall", ignore all occurences of: 
#         "^cellular_component^cell wall" 
#         "^cellular_component^plant-type cell wall"
#
##########################################################
# Input Parameters. Edit these as required:

# The path to the folder containing the file
path <- "/home/benrancourt/Downloads"

# Input file name (it must be a .csv file)
inputFile <- "ExamplefileforBen.csv"

# The character used as the field delimiter in the .csv file. 
# Use "\t" for tab delimited, or "," for comma delimited fields.
fieldDelimiter <- "\t" # 

# A list of key phrases. Eg:
# keyPhrases <- c("phrase one", "phrase two", "phrase three", "etc")
keyPhrases <- c("cell wall", "DFIRtrinTRINITY", "molecular_function")



##############################################################
# Excecution Code. There should be no need to edit anything 
# below this line.

cellWall <- "cell wall"
cellWallException1 <- "^cellular_component^cell wall"
cellWallException2 <- "^cellular_component^plant-type cell wall"

input <- readLines( paste(path.expand(path),inputFile,sep="/") )

linesToKeep <- logical(length=length(input)) # initializes all elements to FALSE
linesToKeep[1] <- TRUE # assuming the first line is a header

for(lineNum in 2:length(input)) {
  line <- input[lineNum]

  for(phrase in keyPhrases) {
    if (phrase == cellWall) {
      line <- gsub(cellWallException1, "", line, fixed=TRUE)
      line <- gsub(cellWallException2, "", line, fixed=TRUE)  
    }    
    contains <- grepl(phrase,line,fixed=TRUE)
    if(contains) {
      linesToKeep[lineNum] <- TRUE
      break
    }    
  } # end inner for-loop

} # end outer for-loop


input <- input[linesToKeep]

outputFileName <- gsub(".csv","-filtered.csv",inputFile,fixed=TRUE)

writeLines(input, paste(path.expand(path),outputFileName,sep="/") )

write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())



