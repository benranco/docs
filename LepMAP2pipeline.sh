#!/bin/bash

# TODO:
# - overall script info, requirements, how to run, the iterative process, why tabulate... must be done separately, because of out of memory error in OrderMarkers.
# - check code to determine how many instances to run in parallel
# - choose to run OrderMarkers only on specific LGs
# - create a repository containing LepMAP2 software and any other requirements

##########################################################################################
# This script runs the LepMAP2 software on an input .csv file containing SNP data in the 
# form of a chi square root table. 
#
# Software downloads and background information on LepMAP2 can be found at:
# https://sourceforge.net/p/lepmap2/wiki/browse_pages/
#
#The purpose of this pipeline is to divide the SNPs into
# linkage groups/chromosomes. This is handled by LepMAP2 modules. However, our
# manipulation of the data is not finished with the completion of this pipeline. 
#
# Following this pipeline, you must complete the process by manually running the R script
# "tabulate_OrderMarkers_results.R", using the files generated by the LepMAP2 pipeline as
# input. The reason you must run this manually is because occasionally during the LepMAP2 
# OrderMarkers module, which is run in parallel on different linkage groups, an 
# out-of-memory error may occur on one or more linkage group, which will require you to 
# re-run the OrderMarkers module on that linkage group. So before running 
# "tabulate_OrderMarkers_results.R", check your OrderMarkers output files (their names
# will end in ".SA.txt") to make sure none of them are empty. If one is empty, re-run 
# this pipeline by changing the input parameters to specify that you only want to run the
# OrderMarkers module, and that you only want to run it on the particular linkage group
# that needs to be redone.
#
# When you are ready to run "tabulate_OrderMarkers_results.R", see the documentation 
# that script for instruction on how to run it, and for a description of what it does.
#
# Once you have completed the entire process, including "tabulate_OrderMarkers_results.R",
# you will want to examine your linkage groups. Sometimes you might want to refine your
# linkage groups, in which case you may choose to run a second iteration of this entire
# process using the files generated by "tabulate_OrderMarkers_results.R" as the input 
# for this pipeline. In this case, since it is the second time through this pipeline, you
# will want to skip some of the LepMAP2 modules. See "Choosing which LepMAP2 modules to 
# run" in the Input Parameters section for more information on which input modules to run.
#
##########################################################################################
# To run LepMAP2pipeline:
#
# - (1) First, make sure your input .csv file is in the following format.
#
# Input .csv file format:
# The chi square root table consists of one column listing SNP ids, and a column for 
# each sample. The data in the chi square root table is represented as either "H", "A", or
# something ("-" or "NA" for example, you can specify which you use in the 
# missingDataIndicator input parameter below) to represent missing or irrelevant data. The 
# most frequent genotype in an SNP is represented as H and the second most frequent type 
# is represented as A. 
#
# The LepMAP2 software does not actually accept data in this .csv format, so the first 
# step in this pipeline is to convert the .csv file into a .linkage file which LepMAP2 can
# read. This is done by automatically calling the R script "convert_csv_to_linkage.R". 
# The .linkage file is named after the .csv and is saved in the same folder as the input 
# .csv file.
#
# In fact, all files generated during this pipeline are saved to the same folder as the 
# input .csv file, which you specify below in the input parameters.
#
# - (2) Second, edit the input parameters (see the Input Parameters section below), and 
# save the file. It is good practice to make a copy of the file first and edit the copy
# instead of the original. If you make a copy with a new file name, use the new file 
# name in the commands below instead of the original script name.
#
# - (3) Third, run the script. To do this, from a command-line Terminal, make sure you are
# in the same folder as the script, and then type:
#   ./LepMAP2pipeline.sh
# (if you've renamed the script, use your filename instead of "LepMAP2pipeline.sh")
#
##########################################################################################
# Input Parameters:

# lepMap2Bin: The absolute file path of the "bin" directory in your LepMAP2 installation
# folder. This is where all the LepMAP2 modules are kept.
lepMap2Bin="/home/benrancourt/Desktop/LepMAP2/binary/bin"


# dataPath: The directory path of the folder that any input files are in. This is also the 
# folder that the output files will be saved to. Since lots of output files are generated,
# I recommend using a dedicated directory that doesn't have anything else in it.
dataPath="/home/benrancourt/Downloads/r45-LepMAP2-final-copyForTesting/r45-inThreeStages/test2"

# inputCsvFilenameWithoutExtension: The base name of the input .csv file, without the
# ".csv" postfix. This name will be used to read the .csv input file, but it will also be
# used as a basis for all the subsequent output files that will be generated by this pipeline.
inputCsvFilenameWithoutExtension="r45-60592-p0.01"

# familyName: This is probably the least important input parameter of this script. This is 
# used as  a name for your input dataset when converting your input .csv file to .linkage 
# formatto be used as input for the LepMAP2 modules. The name you pick is important only
# to you.
familyName="test"

# missingDataIndicator: This parameter relates to the data in the original .csv data 
# file that you supplied. Some of the data in your original input .csv will likely be
# missing, in which case it might be represented as NA or "-" or maybe even be empty,
# or something like that to indicate that it is missing.
# Indicate what your .csv file uses to represent missing data.
missingDataIndicator="NA"

# Choosing which LepMAP2 modules to run:
# You may wish to skip some of these LepMAP2 modules for any particular pass of the LepMAP2
# process. For example, if you're running the LepMAP2 pipeline on the output of a previous
# run of the LepMAP2 pipeline and using that output as input for this new run, you will
# probably want to skip the Filtering, SeparateChromosomes and JoinSingles Modules, since
# the markers have already been organized into linkage groups, and you want to make sure 
# you use the same groups as before. If you do skip these modules and go straight to the
# OrderMarkers module, make sure you have a map file (generated in tabulate_OrderMarkers_results.R
# along with the final .csv output file) to accompany your .csv input file. 
#
# Some possible scenarios in which you might want to skip certainmodules:
# - skip the Filtering module because you've determined Filtering isn't needed
# - skip Filtering, SeparateChromosomes and JoinSingles because you've already run them and you just want to 
#   rerun OrderMarkers to see if you get better results (the OrderMarkers results are slightly different each 
#   time)
# - skip Filtering, SeparateChromosomes and JoinSingles because you've finished a complete run of the LepMAP2 
#   pipeline, including the OrderMarkers module, and you've run the tabulate_OrderMarkers_results.R script on 
#   the OrderMarkers output, and you now want to do a second pass of the OrderMarkers module using the files 
#   generated by tabulate_OrderMarkers_results.R as the input. In this case you should also make sure you use 
#   the new map file generated by tabulate_OrderMarkers_results.R in addition to the new .csv data file.
#
# To indicate whether to skip or do any of the following modules, use:
# 1==do the module, 0==skip the module
doFiltering=0
doSeparateChromosomes=1
doJoinSingles=0
doOrderMarkers=0

# sizeLimitSeparateChr: This is used as an input parameter for the SeparateChromosomes
# module, which forms the Linkage Groups. Any Linkage Group with fewer than sizeLimit
# number of markers will be removed.
sizeLimitSeparateChr=10

# lodLimitSeparateChr: This is used as the LOD score limit input parameter for the
# SeparateChromosomes module.
# See LepMAP2 online documentation for a minimal description of how it uses pair-wise LOD
# scores between markers to assign markers into linkage groups.
lodLimitSeparateChr=10

# lodLimitJoinSingles: This is used as the LOD score limit input parameter for the
# JoinSingles module.
lodLimitJoinSingles=6


# End of input parameters.
##########################################################################################









# ############################################################
# These parameters don't need to be edited from use to use:

# make sure to end with a "/"
dataLocationWithFinalSlash=$dataPath"/"
lepMap2Bin=$lepMap2Bin"/"

inputFileNameMainPart=$inputCsvFilenameWithoutExtension
inputFileNameDotSuffix=".linkage"


logFileNameSuffix="-lepMAP2-log.txt"

# output file naming conventions
filteredOutputSuffix="-filtered.linkage"

mainInputFileSuffix=$inputFileNameDotSuffix
separateChromosomesOutputSuffix="-map.txt"
joinSinglesOutputSuffix="-map_js.txt"

orderMarkersOutputSuffixPt1="-map_js-chr"
orderMarkersOutputSuffixPt2=".SA.txt"
orderMarkersLogFileNameSuffixPt1="-orderMarkers-chr"
orderMarkersLogFileNameSuffixPt2="-lepMAP2-log.txt"

if [[ $doFiltering -eq 1 ]]
then
  mainInputFileSuffix=$filteredOutputSuffix
  separateChromosomesOutputSuffix="-filtered-map.txt"
  joinSinglesOutputSuffix="-filtered-map_js.txt"

  orderMarkersOutputSuffixPt1="-filtered-map_js-chr"
  orderMarkersLogFileNameSuffixPt1="-filtered-orderMarkers-chr"
fi

dl=$dataLocationWithFinalSlash

mainLog=$dl$inputFileNameMainPart$logFileNameSuffix


# ############################################################
# build the commands using above parameters

#default dataTolerance is 0.01
filterCommand="java -cp "$lepMap2Bin" Filtering data="$dl$inputFileNameMainPart$inputFileNameDotSuffix

separateChromosomesCommand="java -cp "$lepMap2Bin" SeparateChromosomes data="$dl$inputFileNameMainPart$mainInputFileSuffix" lodLimit="$lodLimitSeparateChr" sizeLimit="$sizeLimitSeparateChr

joinSinglesCommand="java -cp "$lepMap2Bin" JoinSingles "$dl$inputFileNameMainPart$separateChromosomesOutputSuffix" data="$dl$inputFileNameMainPart$mainInputFileSuffix" lodLimit="$lodLimitJoinSingles

orderMarkersCommandChooseLG="java -cp "$lepMap2Bin" OrderMarkers data="$dl$inputFileNameMainPart$mainInputFileSuffix" map="$dl$inputFileNameMainPart$joinSinglesOutputSuffix" alpha=0.1 polishWindow=100 filterWindow=10  sexAveraged=1 chromosome="

#-----


# ############################################################
function printMainLogFileDivider {
  echo " " | tee -a $mainLog  
  echo "===========================================================================" >> $mainLog
  echo "===========================================================================" >> $mainLog
  date | tee -a $mainLog
}


# ############################################################
function getNumLinkageGroups {
  # parameter $1 should be the name of the map file we will be checking to determine number of linkage groups

  # number of linkage groups = the number of unique lines (minus the header) in the file
  numLinkageGroups=`sort -u $dl$1 | wc -l` 
  numLinkageGroups=$((numLinkageGroups-1)) # - 1 because of the header line
}


# ############################################################
function countMarkersInLinkageGroups {
  ### Count the number of markers in each Linkage Group

  mapFileName=$1 # $1 retrieves the first input parameter to this function

  printMainLogFileDivider
  getNumLinkageGroups $mapFileName

  echo "Number of Linkage Groups in "$mapFileName": "$numLinkageGroups" (might be only "$((numLinkageGroups-1))" if some markers are in 0=no linkage group)" | tee -a $mainLog  

  for i in `seq 0 $numLinkageGroups`
  do 
    count=`grep -c "^$i$" $dataLocationWithFinalSlash$mapFileName`
    echo "Number of markers in LG "$i": "$count | tee -a $mainLog  
  done

}



# ############################################################
# execute the commands:

echo "Converting $inputCsvFilenameWithoutExtension.csv to .linkage format." | tee $mainLog
echo " " >> $mainLog

Rscript convert_csv_to_linkage.R $dataPath $inputCsvFilenameWithoutExtension $familyName $missingDataIndicator

 
echo "Beginning processing $inputFileNameMainPart$inputFileNameDotSuffix with LepMAP2 modules." | tee $mainLog
echo " " >> $mainLog

### Filtering 
if [[ $doFiltering -eq 1 ]]
then
  printMainLogFileDivider  
  echo "Beginning LepMAP2 Filtering module. Output filename: $inputFileNameMainPart$filteredOutputSuffix." | tee -a $mainLog
  echo " " >> $mainLog

  $filterCommand > $dl$inputFileNameMainPart$filteredOutputSuffix 2>> $mainLog
else
  echo " " | tee -a $mainLog
  echo "doFiltering==0, therefore skipping LepMAP2 Filtering module." | tee -a $mainLog
fi

### SeparateChromosomes
if [[ $doSeparateChromosomes -eq 1 ]]
then
  printMainLogFileDivider  
  echo "Beginning LepMAP2 SeparateChromosomes module. Output filename: $inputFileNameMainPart$separateChromosomesOutputSuffix." | tee -a $mainLog
  echo " " >> $mainLog

  $separateChromosomesCommand > $dl$inputFileNameMainPart$separateChromosomesOutputSuffix 2>> $mainLog

  ### Count the number of markers in each Linkage Group in the SeparateChromosomes map file
  countMarkersInLinkageGroups $inputFileNameMainPart$separateChromosomesOutputSuffix
else
  echo " " | tee -a $mainLog
  echo "doSeparateChromosomes==0, therefore skipping LepMAP2 SeparateChromosomes module." | tee -a $mainLog
fi


### JoinSingles
if [[ $doJoinSingles -eq 1 ]]
then
  printMainLogFileDivider  
  echo "Beginning LepMAP2 JoinSingles module. Output filename: $inputFileNameMainPart$joinSinglesOutputSuffix." | tee -a $mainLog
  echo " " >> $mainLog

  $joinSinglesCommand > $dl$inputFileNameMainPart$joinSinglesOutputSuffix 2>> $mainLog

  ### Count the number of markers in each Linkage Group after JoinSingles
  countMarkersInLinkageGroups $inputFileNameMainPart$joinSinglesOutputSuffix
else
  echo " " | tee -a $mainLog
  echo "doJoinSingles==0, therefore skipping LepMAP2 JoinSingles module." | tee -a $mainLog
fi


### =========================================
### OrderMarkers SA for each linkage group, in parallel:
if [[ $doOrderMarkers -eq 1 ]]
then

  printMainLogFileDivider  
  echo "Beginning LepMAP2 OrderMarkers SA modules. See individual Order Markers log files for each Linkage Group." | tee -a $mainLog

  if [[ $doJoinSingles -ne 1 ]]
  then
    # Just for information purposes, count the number of markers in each Linkage Group
    # in the input map file, if we haven't already done it.
    countMarkersInLinkageGroups $inputFileNameMainPart$joinSinglesOutputSuffix
  fi

  getNumLinkageGroups $inputFileNameMainPart$joinSinglesOutputSuffix
  counter=0
  ncore="$(( ($(grep -c ^processor /proc/cpuinfo) -4) ))"

  echo "running parallel processing"

  for i in `seq 1 $numLinkageGroups`
  do
    if [ $counter -ge $ncore ]; then
      wait
      #echo "wait reset"
      counter=0
    fi 

    date > $dl$inputFileNameMainPart$orderMarkersLogFileNameSuffixPt1$i$orderMarkersLogFileNameSuffixPt2
    echo "Beginning LepMAP2 OrderMarkers SA module. Output filename: $dl$inputFileNameMainPart$orderMarkersOutputSuffixPt1$i$orderMarkersOutputSuffixPt2." | tee -a $dl$inputFileNameMainPart$orderMarkersLogFileNameSuffixPt1$i$orderMarkersLogFileNameSuffixPt2
    echo " " >> $dl$inputFileNameMainPart$orderMarkersLogFileNameSuffixPt1$i$orderMarkersLogFileNameSuffixPt2

    $orderMarkersCommandChooseLG$i > $dl$inputFileNameMainPart$orderMarkersOutputSuffixPt1$i$orderMarkersOutputSuffixPt2 2>> $dl$inputFileNameMainPart$orderMarkersLogFileNameSuffixPt1$i$orderMarkersLogFileNameSuffixPt2 &
    counter=$(( $counter + 1 ))
  done
  wait

else
  echo " " | tee -a $mainLog
  echo "doOrderMarkers==0, therefore skipping LepMAP2 OrderMarkers module." | tee -a $mainLog
fi


### =========================================



### Done
printMainLogFileDivider

echo "Finished processing." | tee -a $mainLog

echo  -e "\033[33;5;7mPipeline $name has finished\033[0m"


 

