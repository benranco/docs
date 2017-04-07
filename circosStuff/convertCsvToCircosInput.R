# ########################################################
# To install Circos and the Circos-Tools package (on CentOS), download the current version of each from the Circos download page and unzip them. Move the "tools" folder and all its contents from the unzipped circos-tools location into the "circos" folder of the unzipped circos location.
#
# To install perl modules on CentOS for Circos, from command-line:
# perl -v # to determine if your perl version is at least 5.8 or newer
# navigate to whichever directory you extracted circos to, and navigate to its bin directory  (eg. /home/benrancourt/Desktop/circos/circos-0.69-3/bin), then:
# ./circos modules
# copy the list of modules to a text file
# sudo yum install perl-CPAN
# sudo yum install gd-devel # this includes /usr/bin/gdlib-config, which is needed by some CPAN modules
# sudo cpan Module::Build
# sudo perl -MCPAN -e shell # opens the CPAN shell
# to install the missing modules from the list, for each module, type:
# install modulename
# (eg: install Math::Bezier)
# install module Statistics::Descriptive, which is needed by circos' tableviewer tool (distributed in the tools package) but which might not be included in the list of modules needed by circos:
# install Statistics::Descriptive 
# type q to quit
# make sure you now have all the necessary modules:
# ./circos modules
# if some modules are still missing, go back into the CPAN shell and try installing them again, then read the ouput of the install operation and look for errors, such as modules that are required by this module which for some reason weren't automatically installed, or required non-perl linux libraries that need to be installed via yum. Then try installing those required modules, or resolving those dependency issues and retry installing the module that didn't work.
# ########################################################

options(stringsAsFactors = FALSE, warn = 1)
options(scipen = 999) # disable scientific notation for large numbers

write("Running convertCsvToCircosInput.R.", stdout())


# ########################################################
# Input Parameters:

path <- "/home/benrancourt/Downloads/Circos/"

# Input file format:
# We assume the input csv file colums are organized: 
# Group1 LG, Group1 position, Group2 LG, Group2 position.
#
# These are all numeric values, with no non-numeric characters (so no "LG"
# in the LG ids, for example).
#inputFileName <- "Lim2,110-Pinast-for-CIRCOS-Feb-7-2017-preCircos.csv"
#inputFileName <- "LPmerged-9511-PitaEcht4492-for-CIRCOS-Feb8-2017-preCircos.csv"
inputFileName <- "To-Ben-4-Circos-3617-e-100.csv"
#inputFileName <- "To-Ben-4-Circos-3093-e0.csv"

# This will be used as an input file by circos, and it's name and location should
# be specified in your circos.conf file.
#dataOutputFileName <- "Lim2,110-Pinast-for-CIRCOS-Feb-7-2017-preCircos-data.txt"
#dataOutputFileName <- "LPmerged-9511-PitaEcht4492-for-CIRCOS-Feb8-2017-preCircos-data.txt"
dataOutputFileName <- "To-Ben-4-Circos-3617-e-100-data.txt"
#dataOutputFileName <- "To-Ben-4-Circos-3093-e0-data.txt"

# This will be used as an input file by circos, and it's name and location should
# be specified in your circos.conf file.
#karyotypeOutputFileName <- "Lim2,110-Pinast-for-CIRCOS-Feb-7-2017-preCircos-karyotype.txt"
#karyotypeOutputFileName <- "LPmerged-9511-PitaEcht4492-for-CIRCOS-Feb8-2017-preCircos-karyotype.txt"
karyotypeOutputFileName <- "To-Ben-4-Circos-3617-e-100-karyotype.txt"
#karyotypeOutputFileName <- "To-Ben-4-Circos-3093-e0-karyotype.txt"

# Group1 identifiers are intended for the data set that will be 
# displayed on the left-hand side of circle graph.

# prefixes to use for the labels on the chart
# Pita_ for loblolly pine, Pipi_ for pinster (maritime pine), and Pifl_ for limber pine.
#chromosomeLabelPrefixGroup1 <- "Pifl-"
#chromosomeLabelPrefixGroup2 <- "Pipi-"
#chromosomeLabelPrefixGroup1 <- "Pifl-"
#chromosomeLabelPrefixGroup2 <- "Pita-"
#chromosomeLabelPrefixGroup1 <- "Pifl-"
#chromosomeLabelPrefixGroup2 <- "Pigl-"
chromosomeLabelPrefixGroup1 <- "Pifl-"
chromosomeLabelPrefixGroup2 <- "Pigl-"

# prefixes for circos to use internally for chromosome identifiers. 
# This shouldn't be changed unless you also change the circos.conf file.
chromosomePrefixGroup1 <- "lg-"
chromosomePrefixGroup2 <- "LG"

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
positionMultiplier <- 10000000


# ########################################################

input <- read.csv(paste0(path,inputFileName),header=TRUE)

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

  row <- data.frame("chr", "-", paste0(chromosomePrefixGroup1, LGsInGroup1[i]), paste0(chromosomeLabelPrefixGroup1, LGsInGroup1[i]), minPos, maxPos, paste0("chr",i))
  names(row) <- c("Chr","Dash","Name","Label","MinPos","MaxPos","Color")

  karyotype <- rbind(karyotype, row)
}

# write output

write.table(karyotype, file= paste0(path,karyotypeOutputFileName), append=FALSE, quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)

write.table(data, file= paste0(path,dataOutputFileName), append=FALSE, quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)



write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())


