


options(stringsAsFactors = FALSE, warn = 1)

write("Running convertCsvToCircosInput.R.", stdout())


# ########################################################
# This script extracts rows that relate to each other from a specific .csv file,
# and outputs two .csv files as a result.
#
# Firstly, it extracts any rows that belong to the same LG (column 10) and
# have the same accession (E-value) (column 5), whose consensus values (column 12) 
# are within 5 of each other. All other rows are discarded. 
#
# Secondly, it simplifies the output still further by creating a second .csv file
# with only one row per LG/accession value row-groups. The contents of this second
# file are:
# Col1: accession (E-value)
# Col2: gene, which is a concatenation of the values in column 1 of the first file, 
#   separated by the pipe "|" character.
# Col3: consensus difference, which is the difference between the consensus values of
#   adjacent rows, concatenated by a " | " if more that two rows are in the group.
# Col4: LG
# Col5: consensus mean of the rows in the group
#  
#
# Input Parameters:

path <- "/home/benrancourt/Downloads/"

inputFileName <- "Mapped9546-Scaffolds-Pila-v15_--8619e-100.csv"
#inputFileName <- "Mapped9546-Scaffolds-Pila-v15_--9546-scaffold-posi.csv"

firstOutputFileName <- "Mapped9546-Scaffolds-Pila-v15_--8619e-100-filtered.csv"
#firstOutputFileName <- "Mapped9546-Scaffolds-Pila-v15_--9546-scaffold-posi-filtered.csv"

secondOutputFileName <- "Mapped9546-Scaffolds-Pila-v15_--8619e-100-simplified.csv"
#secondOutputFileName <- "Mapped9546-Scaffolds-Pila-v15_--9546-scaffold-posi-simplified.csv"

# ########################################################

input <- read.csv(paste0(path,inputFileName),header=TRUE)

# col5 = accesion (E-value), col10 = LG, col12 = consensus

LGs <- sort(unique(input[ , 10]))
#LGsInGroup2 <- sort(unique(input[ , 2]))


allFilteredData <- NULL
allSimplifiedData <- NULL

for (i in 1:length(LGs))
{
  write(paste0("LG ",i), stdout())
  dataInLG <- input[ input[,10]==LGs[i], ]
  tabledAccessionData <- table(dataInLG[,5], useNA="no")
  # get only the accession values that occur more than once in the current LG
  duplicateAccesionNamesInLG <- sort(names(tabledAccessionData[tabledAccessionData>=2]))
  
  for (j in 1:length(duplicateAccesionNamesInLG))
  {
    # get the rows in the current LG who's col5 accession value equals the current accession value
    dataOfAccession <- dataInLG[ dataInLG[,5]==duplicateAccesionNamesInLG[j], ]
    # order the rows by the value in col12 (consensus value)
    dataOfAccession <- dataOfAccession[order(dataOfAccession[,12]), ]

    keepRows <- logical(nrow(dataOfAccession)) # initializes all elements to FALSE
    # keep only those rows whose consensus value in col12 is within 5 of 
    # the consensus in either of the adjacent rows
    if (nrow(dataOfAccession) > 1)
    {
      for (k in 1:nrow(dataOfAccession))
      {
        if (k == 1) {
          if ( dataOfAccession[k+1, 12] - dataOfAccession[k, 12] <= 5 )
            keepRows[k] <- TRUE
        } else if (k == nrow(dataOfAccession)) {
          if ( dataOfAccession[k, 12] - dataOfAccession[k-1, 12] <= 5 ) 
            keepRows[k] <- TRUE
        } else {
          if ( (dataOfAccession[k+1, 12] - dataOfAccession[k, 12] <= 5) | 
               (dataOfAccession[k, 12] - dataOfAccession[k-1, 12] <= 5) )
            keepRows[k] <- TRUE        
        }
          
      }
    }

    dataOfAccession <- dataOfAccession[keepRows, ]
    
    if (nrow(dataOfAccession) > 0)
    {
      # all filtered data
      if (is.null(allFilteredData))
      {
        allFilteredData <- dataOfAccession
      }
      else
      {
        allFilteredData <- rbind(allFilteredData, dataOfAccession)
      }

      # all simplified filtered data

      distance <- character()
      for (k in 2:nrow(dataOfAccession))
      {
        distance <- append(distance, toString(dataOfAccession[k,12] - dataOfAccession[k-1,12]))
      }

      simplifiedData <- cbind(dataOfAccession[1,5], 
                              paste(dataOfAccession[,1], collapse="|"), 
                              paste(distance, collapse=" | "), 
                              dataOfAccession[1,10], 
                              mean(dataOfAccession[,12]))

      if (is.null(allSimplifiedData))
      {
        allSimplifiedData <- simplifiedData
      }
      else
      {
        allSimplifiedData <- rbind(allSimplifiedData, simplifiedData)
      }
    }

  } # duplicate accession names for loop
  
} # end LG for loop

colnames(allSimplifiedData) <- c("accession (E-value)","gene","consensus difference","LG","consensus mean")

write.table(allFilteredData, file= paste0(path,firstOutputFileName), append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

write.table(allSimplifiedData, file= paste0(path,secondOutputFileName), append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)


write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())



