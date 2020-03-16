#!/bin/bash
echo "---------------------------------"

inputDir=$1
fileList=$2
outputDir=$3

#inputDir="./try1"
#outputDir="./redownloaded"
#fileList="./badFiles.txt"

echo ""

for f in $(cat $fileList)
do  
  mv $inputDir/$f $outputDir
done


