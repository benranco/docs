#!/bin/bash

# ###########################################################
# The conversion details are based on info in these two links: 
#    http://seqanswers.com/forums/showthread.php?t=6419
#
#    https://toolshed.g2.bx.psu.edu/repository/display_tool?repository_id=c6dafba9ea5df1a8&tool_config=database%2Fcommunity_files%2F000%2Frepo_39%2Fqseq2fastq%2Fqseq2fastq.xml&changeset_revision=6682236a1432
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
# NOTE: This script will remove all lines which didn't pass the Illumina Filter (see variable 11. in above documentation).
#
# ###########################################################

# Input Parameters:
# these should be the only parameters you need to edit
inFolder="WWP_PN6-Raw-data-2011"
outFolder="WWP_PN6-fastq-data-2011"
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
  
  # this will always remove those reads that didn't pass the Illumina filtering test, beause of: if ($11 > 0)
  cat $inFolder/$f | awk -F '\t' '{gsub(/0/,"N", $7); gsub(/\./,"N", $9); if ($11 > 0) printf("@%s_%04d:%s:%s:%s:%s:%s:%s\n%s\n+\n%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10)}' 2>$logFilename 1>$outFilename 
  
done

echo "Done!"


