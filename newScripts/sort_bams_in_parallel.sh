#!/bin/bash

# ###########################################################


# these are the parameters you need to edit

inFolder="bam_files_from_all_RSEM_with_sampleNames"
outFolder="bam_files_from_all_RSEM_with_sampleNames_sorted"
# Make sure this matches the input file naming conventions.
# If $inFilenameEnd begins with a dash, the dash must be escaped, eg. "\-hello.bam"
inFilenameEnd=".bam"
outFilenameEnd="_sorted.bam"

pathToSamtools="/work2/SNPpipeline/tools/samtools-1.3.1/samtools"


curTime=`date +%Y%m%d%H%M%S`
sort_args=tmp_arg_list_sort_$curTime.txt

freeCpus=2
ncpus="$(( ($(grep -c ^processor /proc/cpuinfo) - $freeCpus)))"


# ############################
echo "Running sort_bams_in_parallel.sh on all "$inFilenameEnd" files in "$inFolder
echo "Setting up the commands..."

# clear/create this file:
> $sort_args
  
for f in $(ls $inFolder | grep "$inFilenameEnd$")
do
  outFilepath="$outFolder/${f/$inFilenameEnd/$outFilenameEnd}"
  echo "Setting up "$f
  
  echo "$outFilepath" "$inFolder/$f" >> $sort_args
  
done


echo "Running the commands in parallel..."

parallel -k -j $ncpus --colsep ' ' "$pathToSamtools" sort -o {1} {2} < $sort_args

rm $sort_args

echo "Done!"



