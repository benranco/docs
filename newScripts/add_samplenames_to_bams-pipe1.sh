#!/bin/bash

# ###########################################################


# these are the parameters you need to edit

inFolder="/work3/pipeline1/SNPpipeline/dataTemp/single"
outFolder="bamsWithSNs_pipe1"
# Make sure this matches the input file naming conventions.
# If $inFilenameEnd begins with a dash, the dash must be escaped, eg. "\-hello.bam"
inFilenameEnd=".bam"
outFilenameEnd="_withSN_pipe1.bam"

pathToSamtools="tools/samtools-1.3.1/samtools"

# These two parameters are used to extract the sample name from the filename, by splitting the filename into segments delimited by $fileNameSplitter and identifying the segment that uniquely identifies the sample.
fileNameSplitter="_"
sampleNameSection=0  # array indexing starts at 0, not 1

curTime=`date +%Y%m%d%H%M%S`
addreplacerg_args=tmp_arg_list_addreplacerg_$curTime.txt
index_args=tmp_arg_list_index_$curTime.txt

freeCpus=2
ncpus="$(( ($(grep -c ^processor /proc/cpuinfo) - $freeCpus)))"


# ############################
echo "Running add_samplenames_to_bams.sh on all "$inFilenameEnd" files in "$inFolder
echo "Setting up the commands..."

# clear/create these two files:
> $addreplacerg_args
> $index_args
  
for f in $(ls $inFolder | grep "$inFilenameEnd$")
do
  outFilepath="$outFolder/${f/$inFilenameEnd/$outFilenameEnd}"
  echo "Setting up "$f
  
  # see this for explanation of how this splits $f by "_": 
  # https://stackoverflow.com/questions/918886/how-do-i-split-a-string-on-a-delimiter-in-bash
  IFS=$fileNameSplitter read -ra ARRAY <<< "$f"
  sampleName=${ARRAY[$sampleNameSection]}
  readGroup="readGroup_"$sampleName
  
  echo "ID:$readGroup" "SM:$sampleName" "$outFilepath" "$inFolder/$f" >> $addreplacerg_args
  echo "$outFilepath" >> $index_args
  
done

#$pathToSamtools addreplacerg --output-fmt BAM 
#$pathToSamtools index 

echo "Running the commands in parallel..."

parallel -k -j $ncpus --colsep ' ' "$pathToSamtools" addreplacerg --output-fmt BAM -r {1} -r {2} -o {3} {4} < $addreplacerg_args

parallel -k -j $ncpus "$pathToSamtools" index {} < $index_args

rm $addreplacerg_args
rm $index_args

echo "Done!"

