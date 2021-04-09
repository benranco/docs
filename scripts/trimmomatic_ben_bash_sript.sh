#!/bin/bash
echo "running run_trimmomatic_on_all_pairs.sh"

# to install java:
#   sudo yum install java-1.8.0-openjdk-devel



# these should be the only parameters you need to edit
inFolder="/work/WBP_MJ9-2013"
outFolder="/work/WBP_NGS_trimmomatic/WBP_MJ9-2013-trim"
trimmomaticFolder="/home/centos/software/trimmomatic/Trimmomatic-0.39"

# make sure the input file naming conventions match this and are consistent for all pairs
inFilenameEnd="R1.fq"
r2inFilenameEnd="R2.fq"
outFilenameBaseEnd="filtered_R.fq"
logFilenameEnd="filtered-log.txt"
testCommand="echo hello"

ncore="$(( ($(grep -c ^processor /proc/cpuinfo) -1) ))"

for f in $(ls $inFolder | grep $inFilenameEnd)
do
  r2InFilename="${f/$inFilenameEnd/$r2inFilenameEnd}"
  outFilenameBase="$outFolder/${f/$inFilenameEnd/$outFilenameBaseEnd}"
  logFilename="$outFolder/${f/$inFilenameEnd/$logFilenameEnd}"
  command="java -jar "$trimmomaticFolder"/trimmomatic-0.39.jar PE -threads "$ncore" "$inFolder"/"$f" "$inFolder"/"$r2InFilename" -baseout "$outFilenameBase" ILLUMINACLIP:"$trimmomaticFolder"/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
  echo "-----------"
  echo "r1 in: "$f
  echo "r2 in: "$r2InFilename
  echo "  out: "$outFilenameBase
  echo "command: "$command
  echo "ncore: "$ncore
#  $testCommand  # replace this with the actual command
  $command
done
