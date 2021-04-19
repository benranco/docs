#!/bin/bash

# ###########################################################
# The conversion details are based on info in these four links: 
# http://seqanswers.com/forums/showthread.php?t=6419
#
# https://toolshed.g2.bx.psu.edu/repository/display_tool?repository_id=c6dafba9ea5df1a8&tool_config=database%2Fcommunity_files%2F000%2Frepo_39%2Fqseq2fastq%2Fqseq2fastq.xml&changeset_revision=6682236a1432
#
# Info on Line 1 of the fastq format available in these two links: 
# https://help.basespace.illumina.com/articles/descriptive/fastq-files/
#
# https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf  (page 17)
#
# Example of output of my own script:
# 
# The following data line in QSEQ file format:
# HS6     104     2       1101    1206    1887    0       2       ACCAGTAACATAAGCTTGAAGGAAGAAACCGAGCATTGAGTACATGGCCAACCTTCCATTTTTAATCTCCTTT.TT    caefeffdffeaadW``Nb]\cbcacecdd_feLfa_]]_BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB    1
# 
# Will be converted to FASTQ file format:
# @HS6_0104:2:1101:1206:1887:N:2
# ACCAGTAACATAAGCTTGAAGGAAGAAACCGAGCATTGAGTACATGGCCAACCTTCCATTTTTAATCTCCTTTNTT
# +
# caefeffdffeaadW``Nb]\cbcacecdd_feLfa_]]_BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
# 
######
# from file A05_N6ABXX_5_3_concat_qseq.txt:
# HS6     103     5       1101    1441    1822    0       3       CCCAGTAGTCCTAGCCGTAAACGATGGATACTAAGTGCTGTGCGTATCGACCCGCGCAGTGCTGTAGCTAACGCGT    gggdggcgege_eeeggegafefdNdfcddggegW`cY`bSdadcf]edeb`b_IYY\[]ddddJ__a`ad`fK_O    1
#
# Will be converted to FASTQ file format:
# @HS6:0103:N6ABXX:5:1101:1441:1822 2:N:0:05
# CCCAGTAGTCCTAGCCGTAAACGATGGATACTAAGTGCTGTGCGTATCGACCCGCGCAGTGCTGTAGCTAACGCGT
# +
# gggdggcgege_eeeggegafefdNdfcddggegW`cY`bSdadcf]edeb`b_IYY\[]ddddJ__a`ad`fK_O
#
# Of special note:
# - The second read in a pair in most of the qseq data is identified by a 3 instead of a 2, but I will be changing it to a 2 because that is required by the fasta specification.
# - I copied the N6ABXX from the file name, I assume it is the "flowcell id".
# - I copied the sample number 05 from the file name.
# - parameter 7 of the qseq data is not used.
#
# About formats:
# 
# QSEQ format: 
# QSEQ files are the output of the Illumina pipeline. These files contain the sequence, corresponding qualities, as well as lane, tile and X/Y position of clusters.
# 
# According to Illumina manual qseq files have the following format:
# 
# 1. Machine name: (hopefully) unique identifier of the sequencer.
# 2. Run number: (hopefully) unique number to identify the run on the sequencer.
# 3. Lane number: positive integer (currently 1-8).
# 4. Tile number: positive integer.
# 5. X: x coordinate of the spot. Integer (can be negative).
# 6. Y: y coordinate of the spot. Integer (can be negative).
# 7. Index: positive integer. No indexing should have a value of 1.
# 8. Read Number: 1 for single reads; 1 or 2 for paired ends.
# 9. Sequence
# 10. Quality: the calibrated quality string.
# 11. Filter: Did the read pass filtering? 0 - No, 1 - Yes.
# 
# FASTQ format: 
# A FASTQ file normally uses four lines per sequence. 
# Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line). 
# Line 2 is the raw sequence letters. 
# Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again. 
# Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence.
# 
# Info on Line 1 of the fastq format available here: 
# https://help.basespace.illumina.com/articles/descriptive/fastq-files/
# https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf  (page 17)
#
# NOTE: This script will remove all lines which didn't pass the Illumina Filter (see variable 11. in above documentation).
#
# ###########################################################

# Input Parameters:
# IMPORTANT: This script assumes the input filenames all begin with Axx_xxxx_restOfName, where xx=sample number and xxxx=flow cell id. E.g. here the sample number is 05 and the flow cell id is N6ABXX:
#     A05_N6ABXX_5_3_concat_qseq.txt

# these should be the only parameters you need to edit
#inFolder="WWP_PN6-Raw-data-2011"
#outFolder="WWP_PN6-fastq-data-2011"
inFolder="."
outFolder="."
# Make sure this matches the input file naming conventions.
# If $inFilenameEnd begins with a dash (e.g. "-qseq.txt") it must be escaped: "\-qseq.txt"
inFilenameEnd="_qseq.txt"
outFilenameEnd=".fastq"
logFilenameEnd="-convert.log"

# ############################
echo "running convert-qseq-to-fastq.sh on all "$inFilenameEnd" files in "$inFolder

for f in $(ls $inFolder | grep $inFilenameEnd)
do
  outFilename="$outFolder/${f/$inFilenameEnd/$outFilenameEnd}"
  logFilename="$outFolder/${f/$inFilenameEnd/$logFilenameEnd}"
  echo "Converting "$f
  
  # see this for explanation of how this splits $f by "_": 
  # https://stackoverflow.com/questions/918886/how-do-i-split-a-string-on-a-delimiter-in-bash
  IFS='_' read -ra ARRAY <<< "$f"
  flowcell=${ARRAY[1]}
  sampleNum=${ARRAY[0]}
  sampleNum=`echo $sampleNum | cut -c 2-`
  
  cat $inFolder/$f | awk -v flow=$flowcell -v sample=$sampleNum -F '\t' '{if ($11 > 0) { gsub(/3/,"2", $8); gsub(/\./,"N", $9); gsub(/0/,"Y", $11); gsub(/1/,"N", $11);  printf("@%s:%s:%s:%s:%s:%s:%s %s:%s:0:%s\n%s\n+\n%s\n",$1,$2,flow,$3,$4,$5,$6, $8,$11,sample,$9,$10) }}' 2>$logFilename 1>$outFilename 

done

echo "Done!"


