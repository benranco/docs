#!/bin/bash

# ###########################################################


# these are the parameters you need to edit

inFolder="bams"
outFolder="bams-rg"
# Make sure this matches the input file naming conventions.
# If $inFilenameEnd begins with a dash, the dash must be escaped, eg. "\-hello.bam"
inFilenameEnd=".bam"
outFilenameEnd="_withRG.bam"

pathToSamtools="tools/samtools-1.3.1/samtools"

# These two parameters are used to extract the sample name from the filename, by splitting the filename into segments delimited by $fileNameSplitter and identifying the segment that uniquely identifies the sample.
fileNameSplitter="."
sampleNameSection=3  # array indexing starts at 0, not 1


# ############################
echo "Running add_samplenames_to_bams.sh on all "$inFilenameEnd" files in "$inFolder

for f in $(ls $inFolder | grep "$inFilenameEnd$")
do
  outFilename="$outFolder/${f/$inFilenameEnd/$outFilenameEnd}"
  echo "Processing "$f
  
  # see this for explanation of how this splits $f by "_": 
  # https://stackoverflow.com/questions/918886/how-do-i-split-a-string-on-a-delimiter-in-bash
  IFS=$fileNameSplitter read -ra ARRAY <<< "$f"
  sampleName=${ARRAY[$sampleNameSection]}
  readGroup="readGroup_"$sampleName
  
  
  $pathToSamtools addreplacerg -r "ID:$readGroup" -r "SM:$sampleName" --output-fmt BAM -o $outFilename $inFolder/$f
  
  $pathToSamtools index $outFilename
  
done

echo "Done!"

