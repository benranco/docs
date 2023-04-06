#!/bin/bash

#####################################################
# Input Parameters: modify these parameters before running the script

# The root path. All subfolders are assumed to be relative to this path.
path="."

# The folder to write files, relative to path.
outFolder="output"

fileExtension=".csv"


####################################################
# Execution:

echo "Running extract_vcf_pt2.sh"

mkdir $path/$outFolder

for fileName in $(ls $path/ | grep "$fileExtension\$")
do
  echo "Processing "$fileName
  
  echo '"CHROM","POS","REF","ALT","GT","DP","RO","AO"' > $path/$outFolder/$fileName".new"
  grep -v -E ',"[0-9]/[2-9]",|,"[2-9]/[0-9]",' $path/$fileName >>  $path/$outFolder/$fileName".new"
	
done 


echo " "
echo "============================================="
echo "Finished."


