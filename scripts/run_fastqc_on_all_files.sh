#!/bin/bash
echo "---------------------------------"

pathToFastqcExecutable="/home/benrancourt/Desktop/FastQC/fastqc"
inputDir="./redownloaded"
fileExtension=".fastq.gz"
outputDir="./fastqc_reports_redownloaded"

echo ""

for f in $(ls $inputDir/*$fileExtension)
do  
  echo "---------------------------------"
  $pathToFastqcExecutable -o $outputDir $f
done


