#!/bin/bash

#####################################################
# IMPORTANT: This assumes but doesn't actually check that all the files with the given filename .extension in the source folder are csv files with one header line, and identical columns in the identical order. 

#####################################################
# Input Parameters: modify these parameters before running the script

# The root path. All subfolders are assumed to be relative to this path.
path="."

# The folder containing the source files, relative to path.
sourceFolder="kaksOut"

# All files in the source folder to be concatenated must have this filename extension:
fileExtension=".kaks"


outFileNameBase="RLK121-Pila408nt-KaKs" # don't include the filename extension
tabExtension=".tab"
csvExtension=".csv"

# true if the innput files are tab delimited, false if comma delimited
isInputTabDelimited=true

####################################################
# Execution:

echo "Running concatCsvFilesInFolder.sh"

if [ "$isInputTabDelimited" = true ] 
then
  outputFileName=$outFileNameBase$tabExtension
else 
  outputFileName=$outFileNameBase$csvExtension
fi



# read the header line of the first file in the source folder, and write it to the output 
# file, erasing any previous contents if it already exists.
ff=$(ls $path/$sourceFolder | grep -m 1 ".kaks")
head -1 $path/$sourceFolder/$ff > $path/$outputFileName

for fileName in $(ls $path/$sourceFolder | grep $fileExtension)
do
  echo "Processing "$fileName
  # append all lines after the first line to the output file
  tail -n +2 $path/$sourceFolder/$fileName >> $path/$outputFileName
	
done 


if [ "$isInputTabDelimited" = true ] 
then
  sed 's/\t/,/g' $path/$outputFileName > $path/$outFileNameBase$csvExtension
fi




echo " "
echo "============================================="
echo "Finished concatCsvFilesInFolder.sh"

