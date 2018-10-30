#!/bin/bash
echo "running run_trimmomatic_on_all_pairs.sh"

# to install java:
#   sudo yum install java-1.8.0-openjdk-devel

# yum search string
## yum search is not a reliable way to find the packages. Use yum provides instead and use a filename that you know should be in the package that you want to find - e.g yum provides '*/File/Glob.pm'
# yum list available 'name*'

# NOTE: this is not the final version. I think I had to change from using a base filename to actual filenames for either the input or the output.

# this should be the only parameter you need to edit
folder="/home/benrancourt/Documents/temp/files"

# make sure the input file naming conventions match this and are consistent for all pairs
inFilenameEnd="R1.fastq.gz"
outFilenameBaseEnd="filtered_R.fastq.gz"
logFilenameEnd="filtered-log.txt"
testCommand="echo hello"

ncore="$(( ($(grep -c ^processor /proc/cpuinfo) -2) ))"

for f in $(ls $folder | grep $inFilenameEnd)
do
  outFilenameBase="${f/$inFilenameEnd/$outFilenameBaseEnd}"
  logFilename="${f/$inFilenameEnd/$logFilenameEnd}"
  command="java -jar /opt/bio/trimmomatic/trimmomatic.jar PE -threads "$ncore" -trimlog "$logFilename" -basein "$f" -baseout "$outFilenameBase
  echo "in : "$f
  echo "out: "$outFilenameBase
  echo "command: "$command
  echo "ncore: "$ncore
  $testCommand  # replace this with the actual command
done


