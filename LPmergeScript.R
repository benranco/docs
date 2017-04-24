##########################################################################################
# This R script uses the LPmerge R package to merge linkage maps (in .csv format) from  
# different populations. The output is a consensus map, or several possible consensus maps 
# if you give more than one input value for the maxInterval parameter (set below, in the  
# Input Parameters section).
#
# Documentation for LPmerge can be found at:
# https://cran.r-project.org/web/packages/LPmerge/
# https://academic.oup.com/bioinformatics/article/30/11/1623/284175/LPmerge-an-R-package-for-merging-genetic-maps-by
# Example use-case based tutorial available at:
# http://potatobreeding.cals.wisc.edu/software
#
# The LPmerge package must be installed before you can run this script. See below for 
# instructions on how to install it on CentOS Linux, if it hasn't already been installed.
#
#
##########################################################################################
# To run this script:
#
#
# - (1) First, you need to make sure your input .csv files are all in one directory, and are
# each formatted properly. Each input file represents one population, and must contain
# three columns, in this order:
#    gene, position, linkage group
# The names of the columns are not important, but the format of the content in each column 
# is: 
#
# gene: The gene names must be of the same format/naming convention (e.g. if one file names a 
# gene *M224448 and the other file names the same gene *224448, they will not be recognized 
# as the same, because their names aren't identical). If the different files don't follow 
# the same gene naming convention, you will need to edit them first.
#
# position: The position data must be numeric.
#
# linkage group: The linkage group names must follow the same naming conventions between files (eg. if 
# linkage group 9 is called LG9 in one file and 9 in another, they will not match). You 
# also cannot have more than one linkage group name per linkage groug (eg. 9a, 9b, 9c will
# not all be treated as linkage group 9, but as three separate linkage groups).
#
#
# - (2) Second, open the script in a text editor, edit the input parameters (see the Input 
# Parameters section below), and save the file. It is good practice to make a copy of the 
# file first and edit the copy in case you mess something up. If you make a copy with a 
# new file name, use the new file name in the commands below instead of LPmergeScript.R.
#
#
# - (3) Third, and finally, execute the command to run the script. You can either do this from 
# within the R console, or by using the commandline tool Rscript. 
# To use the commandline tool Rscript, from the commandline, make sure you are in the 
# same folder as the script, and type:
#    Rscript LPmergeScript.R
# To run the script from within the R console, from the commandline, make sure you are in 
# the same folder as the script, and type:
#    R
# This opens the R console. Then type:
#    source("LPmergeScript.R")
#
#
##########################################################################################
# To install the R package LPmerge (on CentOS):
# 
# From command prompt:
#    sudo yum install glpk
#    sudo yum install glpk-devel
#    sudo yum install glpk-utils
#    sudo R
#    install.packages('Rglpk', repos='http://cran.us.r-project.org')
#    install.packages('LPmerge', repos='http://cran.us.r-project.org')
#    library(LPmerge)
#
##########################################################################################

##########################################################################################
# Input Parameters:

# The directory path of the folder that the input .csv files are in. This is also the 
# folder that the output files will be saved to.
path <- "/home/benrancourt/Downloads/LPmerge"

# This is a list of the file names of the input .csv files, without the ".csv" postfix. The list 
# must be surrounded by the characters "c( )", and the names of the files must be in 
# quotes "", and delimited by a comma , . For example:
# map-names <- c("filename1withoutPostfix","filename2withoutPostfix","filename3withoutPostfix")
map.names <- c("LG9-r45-g7727-for-2ndLPmerge-run-nonNumericLG","LG9-r38-g5006-for-2ndLPmerge-run-nonNumericLG")

# The base name of the output files, without the ".csv" postfix. This name will have 
# "-maxIntervalX.csv" appended to it for each output file, where X is the number of the 
# current max interval.
finalOutputCsvNameBase <- "LPmergeOutput-r45gene7727-r38gene5006-LG9-try3repeated"

# This is an integer list of different max intervals to try when making the concensus map.
# The list must be surrounded by the characters "c( )", and delimited by commas. If the 
# list is a sequential series from numbers x to z, it can be represented in shorthand 
# in the following way: x:z
# For examples, these are two ways of representing the same list of 1 to 4:
# maxInterval <- c(1:4)
# maxInterval <- c(1,2,3,4)
# The LPmerge documentation gives this description for max interval:
# max.interval: A whole number (K in formula below) specifying the
#          maximum interval size between bins to include in the
#          objective function.  An array of numbers can be passed to
#          test different values (one consensus map is produced for each
#          value in the array).
maxInterval <- c(1:10)

# This should just be set to NULL unless you want to give a different weight/importance
# to each input map. If you do want to give a different weight to each map, you can do so 
# with a list of integers, one for each map in the same order as the map.names list, 
# representing the relative weight of each. For example, if you have three input maps, you
# could do:
# weightings <- c(100,200,50)
weightings <- NULL
#weightings <- c(621,244,156)





##########################################################################################
# Execution code:

options(stringsAsFactors = FALSE, warn = 1)

write("Running LPmergeScript.R.", stdout())

library(methods) # the methods package is loaded by default in R, but not in Rscript, so explicitly load it here.

if(!require(LPmerge))
{
  install.packages('LPmerge', repos='http://cran.us.r-project.org')
}

library(LPmerge)


write(paste0("================================================"), stdout())
write("Running LPmergeScript.R on: ", stdout())
write(paste(map.names," ",sep=", "), stdout())
write("Max intervals used will be: ", stdout())
write(paste(maxInterval," ",sep=", "))

#write("We will be excluding chromosome 9 from this run and processing it separately, as the r38 data requires that it be processed in a few different ways.", stdout())

inputCsvs <- list()
finalOutputCsvs <- list()

i <- 1
# read in the input csv's
write("Reading input csv's.", stdout())
for (i in 1:length(map.names)) 
{
  filename <- paste0(map.names[i],".csv")
  input <- read.csv(paste(path.expand(path),filename,sep="/"),header=TRUE)
  # This assumes all input csv's first three columns of data are in the order specified
  names(input)[1] <- "gene"
  names(input)[2] <- "position"
  names(input)[3] <- "chromosome"
  inputCsvs[[i]] <- input
}
names(inputCsvs) <- map.names

# The chromosome id's in the first csv will be used to find the equivalent chromosomes in the other csv's. This assumes that the chromosome ids in each csv are named the same.
chromosomes <- unique(inputCsvs[[1]][ , 3]) 

chromIndex <- 1
# iterate over the list of chromosomes, and run LPmerge for each chromosome
#for (chromIndex in c(12))
#for (chromIndex in c(10)) # the chromosome at index 10 is chromosome 9
#for (chromIndex in c(1:9,11,12)) # the chromosome at index 10 is chromosome 9, which will be dealt with separately in this case because of the r38 data.
for (chromIndex in 1:length(chromosomes))
{
  write(paste0("================================================"), stdout())
  write(paste0("Preparing to run LPmerge for chromosome ",chromosomes[chromIndex],". "), stdout())
  result <- NULL
  inputMapsForChrom <- list()
  i <- 1
  # iterate through the input data inputCsvs and extract their data for the current chromosome
  for (i in 1:length(inputCsvs))
  {
    write(paste0("Extracting data for chromosome ",chromosomes[chromIndex]," from ",map.names[i],". "), stdout())
    curCsv <- inputCsvs[[i]]
    inputMapsForChrom[[i]] <- curCsv[which(curCsv$chromosome==chromosomes[chromIndex]),c(1,2)]
  } # end inner for-loop (extracting data from input csv's)

  names(inputMapsForChrom) <- map.names

  write(paste0("The lengths of the input linkage maps for chromosome ",chromosomes[chromIndex]," are:"), stdout())
link.map.lengths <- unlist(lapply(inputMapsForChrom,function(x){max(x$position)}))
  write(paste0(names(link.map.lengths)), stdout())
  write(link.map.lengths, stdout())
  write(paste0("The mean of the lengths of the input linkage maps for chromosome ",chromosomes[chromIndex]," is:"), stdout())
  write(mean(link.map.lengths), stdout())

  if (is.null(weightings))
  {
    write(paste0("The input linkage maps for chromosome ",chromosomes[chromIndex]," will not be weighted."), stdout())
  }
  else
  {
    write(paste0("The input linkage maps for chromosome ",chromosomes[chromIndex]," are weighted as:"), stdout())
    write(weightings, stdout())
  }

  write(paste0("Running LPmerge for chromosome ",chromosomes[chromIndex],". "), stdout())
  write(paste0("--------------------------"), stdout())
  result <- LPmerge(inputMapsForChrom, max.interval = maxInterval, weights = weightings)
  write(paste0("--------------------------"), stdout())

  write(paste0("Adding output of LPmerge for chromosome ",chromosomes[chromIndex]," to final output tables. "), stdout())
  for (intervalIndex in 1:length(maxInterval))
  {
    # add the results for the current chromosome to the tables containing results for all chromosomes
    if (nrow(result[[intervalIndex]] > 0))
    {
      if (length(finalOutputCsvs) < intervalIndex)
      {
        finalOutputCsvs[[intervalIndex]] <- cbind(chromosome = chromosomes[chromIndex], result[[intervalIndex]])
      }
      else 
      {
        finalOutputCsvs[[intervalIndex]] <- rbind(finalOutputCsvs[[intervalIndex]], cbind(chromosome = chromosomes[chromIndex], result[[intervalIndex]]) )
      }
    }
  } # end inner for-loop (adding output to final output tables)

} # end outer for-loop (iterating through chromosomes)

write(paste0("================================================"), stdout())
write(paste0("Writing final output csv's. "), stdout())
for (intervalIndex in 1:length(maxInterval))
{
  outputFilename <- paste0(finalOutputCsvNameBase,"-maxInterval",maxInterval[intervalIndex],".csv")
  write.csv(finalOutputCsvs[[intervalIndex]], paste(path.expand(path),outputFilename,sep="/"), row.names=FALSE)
}

write(paste0("FINISHED."), stdout())


