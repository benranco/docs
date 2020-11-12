#!/bin/bash

#####################################################
# Input Parameters: modify these parameters before running the script

# The root path. All subfolders are assumed to be relative to this path.
path="."

inFolder="../paraatOutput"
outFolder="kaksOut"
logFolder="kaksLog"

masterOutFileNameBase="RLK121-Pila408nt-KaKs" # don't include the filename extension
tabExtension=".tab"
csvExtension=".csv"

####################################################
# Execution:

echo "Running runKaksCalculatorOnFolder.sh"

for fileName in $(ls $path/$inFolder | grep ".axt")
do
	echo "Processing "$fileName
	KaKs_Calculator -i $path/$inFolder/$fileName -o $path/$outFolder/$fileName.kaks 2>&1 1> $path/$logFolder/$fileName.kaks.log
	
done 

echo " "
echo "============================================="
echo "Finished running KaKs_Calculator, now consolidating all the Kaks output files into a single tab delimited file, and then converting it to a .csv file."




####################################################
# Consolidating all the Kaks output files into a single tab delimited 
# file, and then converting it to a .csv file.

# The folder containing the source files (the KaKs output files), relative to path.
sourceFolder=$outFolder

# All files in the source folder to be concatenated must have this filename extension:
fileExtension=".kaks"


# true if the innput files are tab delimited, false if comma delimited
isInputTabDelimited=true



if [ "$isInputTabDelimited" = true ] 
then
  outputFileName=$masterOutFileNameBase$tabExtension
else 
  outputFileName=$masterOutFileNameBase$csvExtension
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
  sed 's/\t/,/g' $path/$outputFileName > $path/$masterOutFileNameBase$csvExtension
fi




echo " "
echo "============================================="
echo "Finished runKaksCalculatorOnFolder.sh"

