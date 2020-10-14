# ########################################################
# This R script converts the specified input .csv file into two data files which will 
# be read by Circos.
#
# The .csv file (specified in the path and inputFileName input parameters below) must 
# have four columns organized:
#   Group1 LG, Group1 position, Group2 LG, Group2 position.
#
# The two output data files are named:
#   mykaryotype.txt, mydata.txt
# 
# ########################################################

options(stringsAsFactors = FALSE, warn = 1)
options(scipen = 999) # disable scientific notation for large numbers

write("Running convertCsvToCircosInput.R.", stdout())


# ########################################################
# Input Parameters:

args <- commandArgs(trailingOnly = TRUE)

# The directory that your input .csv file is in.
#path <- "/home/benrancourt/Downloads/Circos"
path <- args[1]

# Input file format:
# We assume the input csv file colums are organized: 
# Group1 LG, Group1 position, Group2 LG, Group2 position.
#
# These are all numeric values, with no non-numeric characters (so no "LG"
# in the LG ids, for example).
#inputFileName <- "Lim2,110-Pinast-for-CIRCOS-Feb-7-2017-preCircos.csv"
inputFileName <- args[2]

# This is the directory in which to save the mydata.txt and mykaryotype.txt circos
# input files.  
pathToSaveCircosOutputFiles <- args[3]

# This will be used as an input file by circos, and it's name and location should
# be specified in your circos.conf file.
#dataOutputFileName <- "mydata.txt"
dataOutputFileName <- args[4]

# This will be used as an input file by circos, and it's name and location should
# be specified in your circos.conf file.
#karyotypeOutputFileName <- "mykaryotype.txt"
karyotypeOutputFileName <- args[5]

# Group1 identifiers are intended for the data set that will be 
# displayed on the left-hand side of circle graph.

# prefixes to use for the labels on the chart
# Pita_ for loblolly pine, Pipi_ for pinster (maritime pine), and Pifl_ for limber pine.
#chromosomeLabelPrefixGroup1 <- "Pifl-"
#chromosomeLabelPrefixGroup2 <- "Pipi-"
chromosomeLabelPrefixGroup1 <- args[6]
chromosomeLabelPrefixGroup2 <- args[7]

# prefixes for circos to use internally for chromosome identifiers. 
# This shouldn't be changed unless you also change the circos.conf file.
#chromosomePrefixGroup1 <- "lg-"
#chromosomePrefixGroup2 <- "LG"
chromosomePrefixGroup1 <- args[8]
chromosomePrefixGroup2 <- args[9]

# IMPORTANT - multiplier required to get rid of all decimal places in the positions.
# I believe circos likes to use integers for positional information when 
# creating links, so we use a multiplier to convert the positions to
# integers without losing granularity.
# The chromosomes_units parameter in the circos.conf file (or our circos-noTicks.conf 
# and circos-withTicks.conf) should also correspond to this, and the multiplier 
# parameter in the ticks.conf (or our ticks-noTicks.conf and ticks-withTicks.conf) file 
# should be set to what is required to scale the numbers back down to a range you want 
# to use for the labels, if using labeled ticks.
# You shouldn't need to change this unless you have values with more than 7 decimal places.
#positionMultiplier <- 10000000
positionMultiplier <- strtoi(args[10])

# ########################################################
# Execution code:

input <- read.csv(paste(path.expand(path),inputFileName,sep="/"),header=TRUE)

# Convert NA's to 0? No, don't do this, instead remove all rows/links that contain NA values
#input[is.na(input)] <- 0
# Remove all rows/links that contain NA values:
input <- na.omit(input)

# We assume the input csv file colums are organized: Group1 LG, Group1 position, Group2 LG, Group2 position

LGsInGroup2 <- sort(unique(input[ , 3]))
LGsInGroup1 <- sort(unique(input[ , 1]))

write("LG's in group 1: ", stdout())
write(LGsInGroup1, stdout())
write("LG's in group 2: ", stdout())
write(LGsInGroup2, stdout())

if (length(LGsInGroup2) != length(LGsInGroup1))
{
  write("Warning: there is a different number of LG's in each species. Continuing anyway.", stdout())
}

# organize the data for the main data input file

data <- cbind(
  paste0(chromosomePrefixGroup2, input[,3]),
  input[,4] * positionMultiplier,
  input[,4] * positionMultiplier,
  paste0(chromosomePrefixGroup1, input[,1]),
  input[,2] * positionMultiplier,
  input[,2] * positionMultiplier
)

# create the data for the karyotype input file

karyotype <- data.frame(Chr=character(),Dash=character(),Name=character(),Label=character(),MinPos=numeric(),MaxPos=numeric(),Color=character())

# for each LG in Group2 ascending numeric order
# find max pos
# find min pos
# create row
# rbind row to karyotype df

for (i in 1:length(LGsInGroup2))
{
  maxPos <- max( input[ input[,3]==LGsInGroup2[i], 4 ] ) * positionMultiplier
  minPos <- min( input[ input[,3]==LGsInGroup2[i], 4 ] ) * positionMultiplier
  
  if (minPos == maxPos && minPos > 0) {
    minPos = 0
  }
  
  row <- data.frame("chr", "-", paste0(chromosomePrefixGroup2, LGsInGroup2[i]), paste0(chromosomeLabelPrefixGroup2, LGsInGroup2[i]), minPos, maxPos, paste0("chr",i))
  names(row) <- c("Chr","Dash","Name","Label","MinPos","MaxPos","Color")

  karyotype <- rbind(karyotype, row)
}

# for each LG in Group1 descending numeric order
# find max pos
# find min pos
# create row
# rbind row to karyotype df

for (i in length(LGsInGroup1):1)
{
  maxPos <- max( input[ input[,1]==LGsInGroup1[i], 2 ] ) * positionMultiplier
  minPos <- min( input[ input[,1]==LGsInGroup1[i], 2 ] ) * positionMultiplier

  if (minPos == maxPos && minPos > 0) {
    minPos = 0
  }
  
  row <- data.frame("chr", "-", paste0(chromosomePrefixGroup1, LGsInGroup1[i]), paste0(chromosomeLabelPrefixGroup1, LGsInGroup1[i]), minPos, maxPos, paste0("chr",i))
  names(row) <- c("Chr","Dash","Name","Label","MinPos","MaxPos","Color")

  karyotype <- rbind(karyotype, row)
}

# write output

write.table(karyotype, file= paste(path.expand(pathToSaveCircosOutputFiles),karyotypeOutputFileName,sep="/"), append=FALSE, quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)

write.table(data, file= paste(path.expand(pathToSaveCircosOutputFiles),dataOutputFileName,sep="/"), append=FALSE, quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)


write(paste0("Finished converting input .csv to circos input files mydata.txt and mykaryotype.txt."), stdout())

write(paste0("================================================"), stdout())



