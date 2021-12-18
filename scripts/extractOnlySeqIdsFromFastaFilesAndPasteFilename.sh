#!/bin/bash
echo "extractOnlySeqIdsFromFastaFiles.sh"


###########################################
# Editable parameters

inFolder="."
outFolder="."

# make sure the input file naming conventions match this and are consistent for all files
inFilenameEnd=".fasta"
outFilenamePrefix="seqIds_"
outFilenameEnd=".csv"

requiredNumCommasPerFilename=7

###########################################
# Execution code

for f in $(ls $inFolder | grep $inFilenameEnd"$")
do
  newOutFilename=`echo $outFilenamePrefix$f | sed 's,'"$inFilenameEnd"','"$outFilenameEnd"',g' `

  numCommas=`echo $f | sed "s/[^,]//g" | tr -d "\n" | wc -m`
  padCommas=""
  i=$numCommas
  while [ $i -lt $requiredNumCommasPerFilename ]
  do
    i=$(($i+1))
    padCommas=$padCommas","
  done

  echo "$numCommas commas in filename, adding: $padCommas   Processing $f."
  
  # extract all seq ids from the .fasta file. The grep -oE option only prints the segment of the matching lines that contains the matching string (\S matches any non-whitespace charcter). The sed command replaces all occurences of ">" with "".
  grep -oE "^>\S*" $inFolder/$f | sed "s/>//g" | sed "s/$/,$f/g" | sed "s/$inFilenameEnd//g" | sed "s/$/$padCommas/g" > $outFolder/$newOutFilename

done
