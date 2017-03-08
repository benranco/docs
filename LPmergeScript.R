# ########################################################
# To install the R package LPmerge (on CentOS):
# 
# From command promt:
# sudo yum install glpk
# sudo yum install  glpk-devel
# sudo yum install glpk-utils
# sudo R
# install.packages('Rglpk', repos='http://cran.us.r-project.org')
# install.packages('LPmerge', repos='http://cran.us.r-project.org')
# library(LPmerge)
# ########################################################

options(stringsAsFactors = FALSE, warn = 1)

write("Running LPmergeScript.R.", stdout())

if(!require(LPmerge))
{
  install.packages('LPmerge', repos='http://cran.us.r-project.org')
}

library(LPmerge)

# ########################################################
# Input Parameters:

path <- "/home/benrancourt/Downloads/LPmerge/"

map.names <- c("LG9-r45-g7727-for-2ndLPmerge-run","LG9-r38-g5006-for-2ndLPmerge-run")

finalOutputCsvNameBase <- "LPmergeOutput-r45gene7727-r38gene5006-LG9-try3"
maxInterval <- c(1:10)
weightings <- NULL
#weightings <- c(621,244,156)

# ########################################################

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
  input <- read.csv(paste0(path,filename),header=TRUE)
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
  write.csv(finalOutputCsvs[[intervalIndex]], paste0(path,finalOutputCsvNameBase,"-maxInterval",maxInterval[intervalIndex],".csv"), row.names=FALSE)
}

write(paste0("FINISHED."), stdout())


