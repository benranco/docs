#!/bin/bash

# ###########################################################



baseDir=`pwd`

destinationDrr="/work2/tmp2"
#destinationDrr="/work2/bam_files_from_all_RSEM"

# For the two unusual portions of the shell commands below, see these two links:
#  https://stackoverflow.com/questions/5168071/list-sub-directories-with-ls
#  https://stackoverflow.com/questions/1371261/get-current-directory-or-folder-name-without-the-full-path


# ls -d */*/ lists only directories two levels deep
for d in $(ls -d */*/)
do
  echo "Processing in " $d
  cd $d
  mv bowtie2.bam $destinationDrr/${PWD##*/}.bowtie2.bam
  cd $baseDir
  
  
done


echo "Done!"

