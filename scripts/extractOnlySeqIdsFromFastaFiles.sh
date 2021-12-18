#!/bin/bash
echo "extractOnlySeqIdsFromFastaFiles.sh"


###########################################
# Editable parameters

inFolder="."
outFolder="."

# make sure the input file naming conventions match this and are consistent for all files
inFilenameEnd=".fasta"
outFilenamePrefix="seqIds_"

###########################################
# Execution code

for f in $(ls $inFolder | grep $inFilenameEnd"$")
do
  # extract all seq ids from the .fasta file. The grep -oE option only prints the segment of the matching lines that contains the matching string (\S matches any non-whitespace charcter). The sed command replaces all occurences of ">" with "".
  echo "Processing $f."
  grep -oE "^>\S*" $inFolder/$f | sed 's/>//g' > $outFolder/$outFilenamePrefix$f

done
