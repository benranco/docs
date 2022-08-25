#!/bin/bash
echo "select_seqs_from_fasta.sh"

# NOTE: This is not very efficient compared to the selectSeqsFromFasta.R script that I created.

# these should be the only parameters you need to edit
inputSeqIdFile="./test/longestOrfSeqs-all.txt"
inputFasta="./trin-assemblies-renamedSeqIds/WWP18-nt675,367-Trinity.fasta"
outputFasta="./bigSelected.fasta"  # if this file already exists, it'll be overwritten.


################################################################
echo "creating fasta by extracting selected sequences by id from the input fasta file..." 

touch $outputFasta
echo -n "" >  $outputFasta  # make sure it is empty, since subsequent writes will just use >>

# For explanation of the below while-loop to iterate through the lines in a file, see these: 
# https://stackoverflow.com/questions/1521462/looping-through-the-content-of-a-file-in-bash
# https://unix.stackexchange.com/questions/402750/modify-global-variable-in-while-loop

while IFS="" read -r seq || [ -n "$seq" ]
do  
  startLine=`grep -n "^>$seq" $inputFasta | cut -f1 -d: `
  
  endLine="$" # sed interprets the $ to mean the last line. This is the default endLine.
  
  nextSeqLine=`tail -n +$(($startLine + 1)) $inputFasta | grep -n -m 1 "^>"  | cut -f1 -d: `
    
  if [ -z $nextSeqLine ]; then 
    endLine="$" # sed interprets the $ to mean the last line. This is the default endLine.
  else
    nextSeqLine=$(($startLine + $nextSeqLine))
    endLine=$(($nextSeqLine - 1))
  fi

  #echo "seq: $seq, startLine:$startLine, endLine:$endLine" 
  
  # To print out a range of lines from a file using sed:
  # sed -n 'startline,endline p; endline(+1?) q' in.txt > out.txt
  #  -n suppresses echoing the input as output
  #   p prints out the relevant lines
  #   q exit sed without processing rest of file
  
  sed -n "$startLine,$endLine p; $endLine q" $inputFasta >> $outputFasta
  
done < $inputSeqIdFile

echo "finished"


