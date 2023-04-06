#!/bin/bash

# ###########################################################


# these are the parameters you need to edit

inFolder="bam_files_from_all_RSEM_with_sampleNames_sorted"

# Make sure this matches the input file naming conventions.
# If $inFilenameEnd begins with a dash, the dash must be escaped, eg. "\-hello.bam"
inFilenameEnd=".bam"


pathToSamtools="/work2/SNPpipeline/tools/samtools-1.3.1/samtools"


curTime=`date +%Y%m%d%H%M%S`
index_args=tmp_arg_list_index_$curTime.txt

freeCpus=1
ncpus="$(( ($(grep -c ^processor /proc/cpuinfo) - $freeCpus)))"


# ############################
echo "Running index_bams_in_parallel.sh on all "$inFilenameEnd" files in "$inFolder
echo "Setting up the commands..."

# clear/create this file:
> $index_args
  
for f in $(ls $inFolder | grep "$inFilenameEnd$")
do
  echo "Setting up "$f
  
  echo "$inFolder/$f" >> $index_args
  
done

echo "Running the commands in parallel..."

parallel -k -j $ncpus "$pathToSamtools" index {} < $index_args

rm $index_args

echo "Done!"


