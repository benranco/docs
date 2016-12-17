#args <- commandArgs(trailingOnly = TRUE)
#path <- args[1]

path <- "/home/benrancourt/Downloads/LepMAP3Tests/LepMAP3-r45-60592-p0.01/"

#MAF_CUTOFF <- args[2]
#MAF_CUTOFF <- as.double(MAF_CUTOFF)

options(stringsAsFactors = FALSE, warn = 1)

message("running remove_duplicate_markers.R. Output is a .csv sorted by linkage_group, position, marker_id, error_estimate.")

if(!require(stringr))
{
  install.packages('stringr', repos='http://cran.us.r-project.org')
}

library(stringr)

# #######################################################################
# x1. Read -map_js.txt file
# x2. Create a vector of unique Linkage Group (LG) ids from that file
# x3. Read chrX.SA.txt files (one for each chromosome/linkage group)
# 4. For each LG, create a dataframe of markers in each LG:
#        - Col.1: LG
#        - Col.2: male_position from SA file
#        - Col.3: marker id (in col.1 of csv data) corresponding to that line number (next)
#        - Col.4: marker number from SA file (which is the line number in both map_js.txt and .csv file, excluding the header line from the count. i.e. the first data row would be line 1.)
#        - Col.5: error estimate from SA file
#        - Col.6: duplicate OR phases from SA file
#        - Col.7-x: the data columns from csv file
#        - sort by order of columns given above
# 5. For each unique male_position (Col.2 above):
#        - find the ones that aren't listed as duplicates, 
#        - for each of these, find the ones listed as duplicates that share their genotype (first part of marker id) and their error estimate
#        - then find in that group (including the non-duplicate) the one with the most complete data and save its rownum in a vector which we'll use at the end to take only those rows that we're not deleting.
# 6. Delete all the rows we aren't keeping.
# 7. Delete the duplicate/phases column.
# #######################################################################

# TODO: define what is used to represent "no data" in the csv (NA or "-", etc), and update the code down below accordingly

csvFilenameMinusExtension <- "r45-60592-p0.01"
sepChromFilenameMinusExtension <- "r45-60592-p0.01-map"
joinSinglesFilenameMinusExtension <- "r45-60592-p0.01-map_js"

csvData <- read.csv(paste0(path,csvFilenameMinusExtension,".csv"))
numDataColsInCsvData <- ncol(csvData)-1
colnames(csvData)[1] <- "marker_id"
# we can set the marker_number sequentially like this because the marker_number in the
# Order Markers output file refers to the row number in both the .csv and the map file.
csvData$marker_number <- seq.int(nrow(csvData))

mapData <- read.table(paste0(path,sepChromFilenameMinusExtension,".txt"))
linkageGroups <- sort(unique(mapData)[ , 1])

# this will contain the combined and filtered data from all linkage groups in one table
finalCombinedData <- NULL

#linkageGroups <- 1 # TODO: temporary for testing
for (lg in linkageGroups) 
{
  if (lg == 0) 
  {
    next; # lg 0 is not a linkage group, skip to the next for-loop iteration
  }
  message(paste0("Processing linkage group ",lg))
  # The file output from the LepMAP3 OrderMarkers2 module is organized into space/tab
  # delimited columns, which the below read.table command reads in (it excludes the first few lines because they're metadata and don't match the column pattern):
  orderMarkersFile <- paste0(path,joinSinglesFilenameMinusExtension,"-chr",lg,".SA.txt")
  orderMarkersData <- read.table(orderMarkersFile, col.names = c("marker_number","position","female_position","left_paren","phased_data","right_paren"))
  
  # get rid of unnecessary columns:
  orderMarkersData <- subset(orderMarkersData , select=-c(female_position,left_paren,phased_data,right_paren))
  
  combinedData <- merge(orderMarkersData, csvData, by="marker_number", all=FALSE, sort=FALSE)
  # this sorting order of the combinedData is essential for the successful execution 
  # of the below while-loop. If a different order is desired for presentation (eg. the 
  # marker_id before the position), this can be done at the very end of this 
  # for-loop (which iterates through the linkage groups), after the inner while-loop 
  # has been exited.
  combinedData <- combinedData[order(combinedData$position, combinedData$marker_id), ]
#  print(combinedData[1:10,1:3])
#  write.csv(combinedData, paste0(path,csvFilenameMinusExtension,"postLepMAP3-lg",lg,"-noneRemoved.csv"), row.names=FALSE)

  # we will use these temporary columns below
  combinedData$record_number <- seq.int(nrow(combinedData)) 
  
  # library(stringr) required for str_split_fixed:
  combinedData$genotype <- str_split_fixed(combinedData$marker_id,"-",2)[ ,1]

  index <- 1
  rowsToKeep <- numeric(nrow(combinedData)) # initializes all elements to 0
  rowsToKeepIndex <- 1

  # This while-loop is used to determine which rows to keep and which rows (duplicates)
  # to get rid of. The indexes of the rows to keep are stored in rowsToKeep.
  # NOTE: this while-loop only works in this iterative fashion because combinedData is 
  # already sorted by the columns position, marker_id (from which genotype is derived),
  # and error_estimate.
  while (index <= nrow(combinedData))
  {
    curPos <- combinedData[index, "position"][[1]]
    curGenotype <- combinedData[index, "genotype"][[1]]

    # for our purposes a rowGroup is defined as all rows that match the current position and genotype
    rowGroup <- combinedData[ combinedData$position == curPos & combinedData$genotype == curGenotype , ]
    if (nrow(rowGroup) == 1)
    {
      rowsToKeep[rowsToKeepIndex] <- rowGroup[1, "record_number"]
      rowsToKeepIndex <- rowsToKeepIndex + 1
    }
    else 
    {
      # we're going to keep only the one with the most complete data
      nodataCounts <- NULL
      nodataCounts <- numeric(nrow(rowGroup))
      nodataIndex = 1
      # TODO: define what is used to represent "no data" in the csv (NA or "-", etc), and update the code accordingly
      while (nodataIndex <= nrow(rowGroup))
      {
        nodataCounts[nodataIndex] <- length(grep("^-$", rowGroup[nodataIndex, 4:numDataColsInCsvData+3]))
        nodataIndex <- nodataIndex + 1
      }
      rowsToKeep[rowsToKeepIndex] <- rowGroup[which.min(nodataCounts), "record_number"]
      rowsToKeepIndex <- rowsToKeepIndex + 1
    }

    index <- index + nrow(rowGroup)
  } # end while-loop

  # remove all elements containing 0 from rowsToKeep
  rowsToKeep <- rowsToKeep[!(rowsToKeep %in% c(0))] 
  # remove all rows we don't want to keep
  combinedData <- combinedData[rowsToKeep, ]
  # remove unnecessary and temporary columns
  combinedData <- subset(combinedData , select=-c(record_number, genotype))

  # re-arrange the order of the columns to make more visual sense (column 3 is the marker_id, col1 is marker_number, col2 is position)
  combinedData <- combinedData[ , c(2,3,1,4:ncol(combinedData))]

  # Modify the sorted order of the data slightly, so that it is ordered by marker_id
  # before error_estimate. (It needed to be sorted by error_estimate before 
  # marker_id for the successful execution of the above inner while-loop).
  combinedData <- combinedData[order(combinedData$position, combinedData$marker_id), ]

  # add a linkage_group column to the beginning of the dataframe
  combinedData <- cbind(linkage_group = lg, combinedData)

  # add the current linkage group to the table containing all linkage groups
  if (nrow(combinedData > 0))
  {
    if (is.null(finalCombinedData))
    {
      finalCombinedData <- combinedData
    }
    else 
    {
      finalCombinedData <- rbind(finalCombinedData, combinedData)
    }
  }
} # end for-loop (iterates through the linkage groups)

write.csv(finalCombinedData, paste0(path,csvFilenameMinusExtension,"-postLepMAP3-allChr.csv"), row.names=FALSE)

message("remove_duplicate_markers complete")
