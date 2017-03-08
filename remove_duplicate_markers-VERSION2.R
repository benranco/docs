#args <- commandArgs(trailingOnly = TRUE)
#path <- args[1]

path <- "/home/benrancourt/Downloads/LepMAP2-VERSION2-r45-60592-p0-test1/"

#MAF_CUTOFF <- args[2]
#MAF_CUTOFF <- as.double(MAF_CUTOFF)

options(stringsAsFactors = FALSE, warn = 1)

message("Running remove_duplicate_markers-VERSION2.R.")

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
#        - for each of these, find the ones listed as duplicates that share their gene (first part of marker id) and their error estimate
#        - then find in that group (including the non-duplicate) the one with the most complete data and save its rownum in a vector which we'll use at the end to take only those rows that we're not deleting.
# 6. Delete all the rows we aren't keeping.
# 7. Delete the duplicate/phases column.
# #######################################################################



# TODO: define what is used to represent "no data" in the csv (NA or "-", etc), and update the code down below accordingly

csvFilenameMinusExtension <- "r45-60592-p0.01"
mapFilenameMinusExtension <- "r45-60592-p0.01-map_js"

# This will be used as a header for the map file that is generated to go along with the ouput data file of this routine. It must begin with a "#", to indicate it is a comment.
outputMapFileHeader <- paste0("#Linkage groups of tabulated output of remove_duplicate_markers-VERSION2.R performed on the output of LepMAP2 OrderMarkers on ",csvFilenameMinusExtension," data, ordered by linkage group, marker id.")

csvData <- read.csv(paste0(path,csvFilenameMinusExtension,".csv"))
numDataColsInCsvData <- ncol(csvData)-1
colnames(csvData)[1] <- "marker_id"
# we can set the marker_number sequentially like this because the marker_number in the
# Order Markers output file refers to the row number in both the .csv and the map file.
csvData$marker_number <- seq.int(nrow(csvData))

mapData <- read.table(paste0(path,mapFilenameMinusExtension,".txt"))
linkageGroups <- sort(unique(mapData)[ , 1])

# this will contain the combined and filtered data from all linkage groups in one table
finalUniqueGenes <- NULL

#linkageGroups <- 1 # TODO: temporary for testing
for (lg in linkageGroups) 
{
  if (lg == 0) 
  {
    next; # lg 0 is not a linkage group, skip to the next for-loop iteration
  }
  message(paste0("Processing linkage group ",lg))
  # The file output from the LepMAP2 OrderMarkers module is organized into space/tab
  # delimited columns, which the below read.table command reads in (it excludes the first few lines because they're metadata and don't match the column pattern):
  orderMarkersFile <- paste0(path,mapFilenameMinusExtension,"-chr",lg,".SA.txt")
  orderMarkersData <- read.table(orderMarkersFile, col.names = c("marker_number","position","female_position","left_paren","error_estimate","right_paren","duplicate_or_phases"))
  
  # get rid of unnecessary columns:
  orderMarkersData <- subset(orderMarkersData , select=-c(female_position,left_paren,right_paren))
  
# #######################################################################
#(- create a new csv that removes all duplicates with error estimates that aren't the lowest)
#x- order all remaining data (including duplicates) by lg, gene, marker/snp id, error estimate, position
#??- remove all those whose error estimate is not the lowest??
#- for each gene:
#- find the highest and lowest position
#- calculate the distance between the highest and lowest
#- calculate the median position
#- calculate the median distance between adjacent positions (standard error?)
#- find the position closest to the median position
#   - if two equi-distant positions on either side of median
#      - pick the position with the most snps/markers at it, or if tied
#      - or pick the one with the most 0 error estimates, or if tied
#      - or pick the one with the snp that has the lowest error estimate, or if tied
#      - or pick the one with the snp that has the most complete data, or if tied
#      - or pick the one with the nearest/largest group of markers at an adjacent position, or if tied
#      - just pick the lower position
#- choose the snp/marker in that position with the (lowest error and) most complete data to represent the gene
#- remove all other markers of that gene
#- save in the table: lg, marker position, gene, num markers, dist betw highest and lowest pos, median dist between adjacent positions (std error?), median position, marker id, marker number, marker error est, csv data
#- remove all genes that occur in multiple LGs (create table: gene id, snp id, lg id, position, error estimate)
#      - after establishing that there's only one occurence of each gene per LG, search the gene column for any gene that appears more than once
#      - remove those rows to their own table
# #######################################################################

  combinedData <- merge(orderMarkersData, csvData, by="marker_number", all=FALSE, sort=FALSE)
  # this sorting order of the combinedData is essential for the successful execution 
  # of the below while-loop. If a different order is desired for presentation (eg. the 
  # marker_id before the error_estimate), this can be done at the very end of this 
  # for-loop (which iterates through the linkage groups), after the inner while-loop 
  # has been exited.
  combinedData <- combinedData[order(combinedData$marker_id, combinedData$error_estimate, combinedData$position), ]
#  print(combinedData[1:10,1:5])
#  write.csv(combinedData, paste0(path,csvFilenameMinusExtension,"-lg",lg,".csv"), row.names=FALSE)
  
  # library(stringr) required for str_split_fixed:
  # TODO: set this in an input parameter
  combinedData$gene <- str_split_fixed(combinedData$marker_id,"_",2)[ ,1]

  # -----

  uniqueGenes <- unique(combinedData$gene)
  # create dataframe uniqueGenes(gene, total_markers, highest_pos, lowest_pos, median_pos, total_distance,  median_distance, marker_id)
  numUniqueGenes <- length(uniqueGenes)
  uniqueGenes <- data.frame(gene = uniqueGenes, total_markers = numeric(numUniqueGenes), highest_pos = numeric(numUniqueGenes), lowest_pos = numeric(numUniqueGenes), median_pos = numeric(numUniqueGenes), total_distance = numeric(numUniqueGenes), median_distance = numeric(numUniqueGenes), marker_id = numeric(numUniqueGenes))


  for (i in 1:nrow(uniqueGenes))
  {
    snpsInGene <- subset(combinedData, gene==uniqueGenes[i, "gene"])

    if (nrow(snpsInGene) == 1) 
    {
      uniqueGenes[i, "total_markers"] <- 1
      uniqueGenes[i, "highest_pos"] <- snpsInGene[1, "position"]
      uniqueGenes[i, "lowest_pos"] <- snpsInGene[1, "position"]
      uniqueGenes[i, "median_pos"] <- snpsInGene[1, "position"]
      uniqueGenes[i, "total_distance"] <- 0
      uniqueGenes[i, "median_distance"] <- 0
      uniqueGenes[i, "marker_id"] <- snpsInGene[1, "marker_id"]
    }
    else
    {
      # ordering snpsInGene this way is necessary for calculating distancesBetweenPositions
      snpsInGene <- snpsInGene[order(snpsInGene$position), ] 

      uniqueGenes[i, "total_markers"] <- nrow(snpsInGene)
      uniqueGenes[i, "highest_pos"] <- max(snpsInGene$position)
      uniqueGenes[i, "lowest_pos"] <- min(snpsInGene$position)
      uniqueGenes[i, "median_pos"] <- median(snpsInGene$position)
      uniqueGenes[i, "total_distance"] <- uniqueGenes[i, "highest_pos"] - uniqueGenes[i, "lowest_pos"]
      
      # calculate the median distance between positions
      distancesBetweenPositions <- numeric(nrow(snpsInGene)-1)
      for (k in 2:nrow(snpsInGene))
      {
        distancesBetweenPositions[k-1] <- snpsInGene[k, "position"] - snpsInGene[k-1, "position"]
      }
      uniqueGenes[i, "median_distance"] <- median(distancesBetweenPositions)
      
      # find the snp(s) whose position match or are nearest to the median position:

      snpsOfMedianPosition <- subset(snpsInGene, position==uniqueGenes[i, "median_pos"])

      # if there are no snps whose position match the median position
      if ( is.null(snpsOfMedianPosition) | nrow(snpsOfMedianPosition) == 0 )
      {
        # find the snps directly above and below the median point, and get their positions.
        # snpsInGene is already ordered by position, so we can do this:
        indexOfSnpBelowMedian <- nrow(snpsInGene) %/% 2
        indexOfSnpAboveMedian <- indexOfSnpBelowMedian + 1
        positionBelowMedian <- snpsInGene[indexOfSnpBelowMedian, "position"]
        positionAboveMedian <- snpsInGene[indexOfSnpAboveMedian, "position"]

        snpsOfPositionBelowMedian <- subset(snpsInGene, position==positionBelowMedian)
        snpsOfPositionAboveMedian <- subset(snpsInGene, position==positionAboveMedian)
        
        # if one group is significantly larger than the other, use it
        if (nrow(snpsOfPositionBelowMedian) / nrow(snpsOfPositionAboveMedian) >= 1.5)
        {
          snpsOfMedianPosition <- snpsOfPositionBelowMedian
        }
        else if (nrow(snpsOfPositionAboveMedian) / nrow(snpsOfPositionBelowMedian) >= 1.5)
        {
          snpsOfMedianPosition <- snpsOfPositionAboveMedian
        }
        # else, use the following logic to pick a group
        else 
        {
          # remove any snp from each group whose error estimate is not equal to the lowest
          snpsOfPositionBelowMedian <- snpsOfPositionBelowMedian[order(snpsOfPositionBelowMedian$error_estimate), ]
          snpsOfPositionBelowMedian <- subset(snpsOfPositionBelowMedian, error_estimate==snpsOfPositionBelowMedian[1, "error_estimate"])

          snpsOfPositionAboveMedian <- snpsOfPositionAboveMedian[order(snpsOfPositionAboveMedian$error_estimate), ]
          snpsOfPositionAboveMedian <- subset(snpsOfPositionAboveMedian, error_estimate==snpsOfPositionAboveMedian[1, "error_estimate"])

          # if one group is significantly larger than the other, use it
          if (nrow(snpsOfPositionBelowMedian) / nrow(snpsOfPositionAboveMedian) >= 1.5)
          {
            snpsOfMedianPosition <- snpsOfPositionBelowMedian
          }
          else if (nrow(snpsOfPositionAboveMedian) / nrow(snpsOfPositionBelowMedian) >= 1.5)
          {
            snpsOfMedianPosition <- snpsOfPositionAboveMedian
          }
          # else use the one with the lowest error_estimate
          else if ( snpsOfPositionBelowMedian[1, "error_estimate"] < snpsOfPositionAboveMedian[1, "error_estimate"] )
          {
            snpsOfMedianPosition <- snpsOfPositionBelowMedian
          }
          else if ( snpsOfPositionAboveMedian[1, "error_estimate"] < snpsOfPositionBelowMedian[1, "error_estimate"] )
          {
            snpsOfMedianPosition <- snpsOfPositionBelowMedian
          }
          # else pick the one that's nearest the median, or if they're tied use the lower one
          else
          {
            distFromUpper <- snpsOfPositionAboveMedian[1, "position"] - uniqueGenes[i, "median_pos"]
            distFromLower <- uniqueGenes[i, "median_pos"] - snpsOfPositionBelowMedian[1, "position"]
            if (distFromLower <= distFromUpper)
            {
              snpsOfMedianPosition <- snpsOfPositionBelowMedian
            }
            else
            {
              snpsOfMedianPosition <- snpsOfPositionAboveMedian
            }
          }         
        }
      } # end of all logic to pick group of snps to represent the median position

      # choose the snp/marker in snpsOfMedianPosition with the (lowest error and) most complete data to represent the gene:

      # remove any snp from each group whose error estimate is not equal to the lowest
      snpsOfMedianPosition <- snpsOfMedianPosition[order(snpsOfMedianPosition$error_estimate), ]
      snpsOfMedianPosition <- subset(snpsOfMedianPosition, error_estimate==snpsOfMedianPosition[1, "error_estimate"])

      if (nrow(snpsOfMedianPosition) == 1)
      {
        uniqueGenes[i, "marker_id"] <- snpsOfMedianPosition[1, "marker_id"]
      }
      else
      {
        # we're going to use the one with the most complete data (least amount of missing data)
        naCounts <- NULL
        naCounts <- numeric(nrow(snpsOfMedianPosition))
        naCountsIndex = 1
        # TODO: define what is used to represent "no data" in the csv (NA or "-", etc), and update the code accordingly
        while (naCountsIndex <= nrow(snpsOfMedianPosition))
        {
          naCounts[naCountsIndex] <- length(grep("^-$", snpsOfMedianPosition[naCountsIndex, 6:numDataColsInCsvData+5]))
          naCountsIndex <- naCountsIndex + 1
        }
        uniqueGenes[i, "marker_id"]  <- snpsOfMedianPosition[which.min(naCounts), "marker_id"]
      }
      
    } # end of logic to populate the one row used to represent the current gene
  } # end inner for-loop


  # join uniqueGenes and combinedData by "marker_id"
  uniqueGenes <- merge(uniqueGenes, combinedData, by=c("gene","marker_id"), all=FALSE, sort=FALSE)

  # -----

  # remove unnecessary columns
  uniqueGenes <- subset(uniqueGenes , select=-c(duplicate_or_phases)) 

  # re-arrange the order of the columns so the marker_id column comes after all the summary columns. New order is: "gene","total_markers","highest_pos","lowest_pos","median_pos","total_distance","median_distance","marker_id","marker_number","position","error_estimate",11:ncol(uniqueGenes)
  uniqueGenes <- uniqueGenes[ , c(1,3,4,5,6,7,8,2,9,10,11,12:ncol(uniqueGenes))]


  # add a linkage_group column to the beginning of the dataframe
  uniqueGenes <- cbind(linkage_group = lg, uniqueGenes)

  # add the current linkage group to the table containing all linkage groups
  if (nrow(uniqueGenes > 0))
  {
    if (is.null(finalUniqueGenes))
    {
      finalUniqueGenes <- uniqueGenes
    }
    else 
    {
      finalUniqueGenes <- rbind(finalUniqueGenes, uniqueGenes)
    }
  }
} # end outer for-loop (iterates through the linkage groups)


# after outer for-loop, remove all rows whose "gene" field occurs more than once in the table, and save their data in a new table
message("Separating genes that occur in more than one linkage group into their own table.")
geneOccurences <- table(finalUniqueGenes$gene)
genesThatOccureMoreThanOnce <- names(geneOccurences[geneOccurences>1])
genesInMultipleLGs <- finalUniqueGenes[finalUniqueGenes$gene %in% genesThatOccureMoreThanOnce, ]
finalUniqueGenes <- finalUniqueGenes[!(finalUniqueGenes$gene %in% genesThatOccureMoreThanOnce), ]

# order the genesInMultipleLGs data by gene first, then LG
genesInMultipleLGs <- genesInMultipleLGs[order(genesInMultipleLGs$gene,genesInMultipleLGs$linkage_group), ]


# output a stripped-down csv file that just has the data from the representative
# marker for each gene
representativeMarkers <- finalUniqueGenes[ , c(9,13:ncol(finalUniqueGenes))]
write.csv(representativeMarkers, paste0(path,csvFilenameMinusExtension,"-postLepMAP2-allChr-uniqueGenes-justTheMarkers.csv"), row.names=FALSE)

# build a map file that preserves the linkage group/gene relationship:
write(outputMapFileHeader, paste0(path,csvFilenameMinusExtension,"-postLepMAP2-uniqueGenes-map_js.txt"), ncolumns=1,append=FALSE)
write(finalUniqueGenes$linkage_group, paste0(path,csvFilenameMinusExtension,"-postLepMAP2-uniqueGenes-map_js.txt"), ncolumns=1,append=TRUE)


write.csv(genesInMultipleLGs, paste0(path,csvFilenameMinusExtension,"-postLepMAP2-allChr-genesInMoreThanOneLG.csv"), row.names=FALSE)

write.csv(finalUniqueGenes, paste0(path,csvFilenameMinusExtension,"-postLepMAP2-allChr-uniqueGenes.csv"), row.names=FALSE)


# To output a version of the output where our LG numbers match Pinus consensus LG numbers, we need to change:
# LG1->2
# LG2->1
# LG3->3
# LG4->8
# LG5->10
# LG6->4
# LG7->12
# LG8->6
# LG9->7
# LG10->11
# LG11->5
# LG12->9a
# LG13->9b

finalUniqueGenesWithUpdateLGs <- finalUniqueGenes

rowsPerOriginalLG <- list()
originalLGs <- unique(finalUniqueGenesWithUpdateLGs$linkage_group)
# record the row numbers of each LG
for (i in originalLGs) # this assumes that the original LGs are all numeric
{
  rowsPerOriginalLG[[ i ]] <- which(finalUniqueGenesWithUpdateLGs$linkage_group %in% i)
}

#update the LG's to match the Pinus consensus LG numbers
for (i in originalLGs) # this assumes that the original LGs are all numeric
{
  if (i == 1) {
    finalUniqueGenesWithUpdateLGs[rowsPerOriginalLG[[i]], "linkage_group"] <- 2
  } else if (i == 2) {
    finalUniqueGenesWithUpdateLGs[rowsPerOriginalLG[[i]], "linkage_group"] <- 1
  } else if (i == 3) {
    finalUniqueGenesWithUpdateLGs[rowsPerOriginalLG[[i]], "linkage_group"] <- 3
  } else if (i == 4) {
    finalUniqueGenesWithUpdateLGs[rowsPerOriginalLG[[i]], "linkage_group"] <- 8
  } else if (i == 5) {
    finalUniqueGenesWithUpdateLGs[rowsPerOriginalLG[[i]], "linkage_group"] <- 10
  } else if (i == 6) {
    finalUniqueGenesWithUpdateLGs[rowsPerOriginalLG[[i]], "linkage_group"] <- 4
  } else if (i == 7) {
    finalUniqueGenesWithUpdateLGs[rowsPerOriginalLG[[i]], "linkage_group"] <- 12
  } else if (i == 8) {
    finalUniqueGenesWithUpdateLGs[rowsPerOriginalLG[[i]], "linkage_group"] <- 6
  } else if (i == 9) {
    finalUniqueGenesWithUpdateLGs[rowsPerOriginalLG[[i]], "linkage_group"] <- 7
  } else if (i == 10) {
    finalUniqueGenesWithUpdateLGs[rowsPerOriginalLG[[i]], "linkage_group"] <- 11
  } else if (i == 11) {
    finalUniqueGenesWithUpdateLGs[rowsPerOriginalLG[[i]], "linkage_group"] <- 5
  } else if (i == 12) {
    if (length(originalLGs) < 13) {
      finalUniqueGenesWithUpdateLGs[rowsPerOriginalLG[[i]], "linkage_group"] <- 9
    } else {
      finalUniqueGenesWithUpdateLGs[rowsPerOriginalLG[[i]], "linkage_group"] <- "9a"
    }
  } else if (i == 13) {
    finalUniqueGenesWithUpdateLGs[rowsPerOriginalLG[[i]], "linkage_group"] <- "9b"
  }  

}

write.csv(finalUniqueGenesWithUpdateLGs, paste0(path,csvFilenameMinusExtension,"-postLepMAP2-allChr-uniqueGenes-PinusConsensusLGs.csv"), row.names=FALSE)

# build a map file that preserves the linkage group/gene relationship:
write(outputMapFileHeader, paste0(path,csvFilenameMinusExtension,"-postLepMAP2-uniqueGenes-PinusConsensusLGs-map_js.txt"), ncolumns=1,append=FALSE)
write(finalUniqueGenesWithUpdateLGs$linkage_group, paste0(path,csvFilenameMinusExtension,"-postLepMAP2-uniqueGenes-PinusConsensusLGs-map_js.txt"), ncolumns=1,append=TRUE)



message("remove_duplicate_markers-VERSION2 complete.")
