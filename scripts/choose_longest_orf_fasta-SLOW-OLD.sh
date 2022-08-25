#!/bin/bash
echo "running choose_longest_orf_fasta.sh"

# What this script should do:
# - for each gene, check which transcript has the longest ORF (in .pep file)
# - if there is a tie for longest ORF, pick the transcript that has the longest DNA sequence (in .fasta file)
# - if no good ORF is generated from TransDecoder, pick the transcript that has the longest DNA sequence (in .fasta file)
# - After this processing, the total sequence number of 675,367 will decrease to 449,052 transcripts to represent 449,052 genes in WWP18 assembly.


##########################################################
# Input parameters. These should be the only parameters you need to edit.

# all input file paths must be relative to this outDir:
outDir="./test2"

# input file paths, all thes must e relative to outDir:
transdecoderFile="../../TransDecoderOutput/WWP18-nt675,367-Trinity.fasta.transdecoder.pep"
trinityFile="../../trin-assemblies-renamedSeqIds/WWP18-nt675,367-Trinity.fasta"
rScriptToSelectSeqsFromFasta="../selectSeqsFromFasta.R"

# this will be used to prefix all files created by this script except the finalOutputFasta file. 
outFilePrefix="test2-WWP18-"

# Final output fasta file. If this file already exists, it'll be overwritten.
finalOutputFasta="selected.fasta"  

##########################################################
# Output files created by this script:

allSeqIds="$outFilePrefix"allSeqIds.txt
allGeneIds="$outFilePrefix"allGeneIds.txt
geneIds100="$outFilePrefix"geneIds-100.txt  # only used for testing with smaller sample sets 
longestOrfSeqsAll="$outFilePrefix"longestOrfSeqs-all.txt
longestOrfSeqsFromPep="$outFilePrefix"longestOrfSeqs-fromPep.txt
longestDNASeqsFromFasta="$outFilePrefix"longestDNASeqs-fromFasta.txt

##########################################################
# Execution code. You shouldn't need to edit below this point.

echo "transdecoderFile: "$outDir/$transdecoderFile
echo "trinityFile: "$outDir/$trinityFile
echo "Location of the R script to extract the select longest sequences from the trinityFile: "$outDir/$rScriptToSelectSeqsFromFasta
echo "Final output fasta file name: "$outDir/$finalOutputFasta

# extract all seq ids from the .fasta file. The grep -oE option only prints the segment of the matching lines that contains the matching string (\S matches any non-whitespace charcter). The sed command replaces all occurences of ">" with "".
echo "Getting sequence ids."
grep -oE "^>\S*" $outDir/$trinityFile | sed 's/>//g' > $outDir/$allSeqIds

# remove the seqId portion of the seq Ids to leave only the gene ids, and keep only one occurene of each gene Id. The sort command is necessary before uniq because uniq only removes contiguous duplicate lines:
echo "Getting gene ids."
sed 's/_i[0-9]*$//g' $outDir/$allSeqIds | sort | uniq > $outDir/$allGeneIds

# For testing with smaller sample sets, uncomment the below line, and 
# replace $allGeneIds with $geneIds100 in the below for-loop.
##head -n 100 $outDir/$allGeneIds > $outDir/$geneIds100

touch $outDir/$longestOrfSeqsAll
touch $outDir/$longestOrfSeqsFromPep
touch $outDir/$longestDNASeqsFromFasta
echo -n "" >  $outDir/$longestOrfSeqsAll  # make sure it is empty, since subsequent writes will just use >>
echo -n "" >  $outDir/$longestOrfSeqsFromPep
echo -n "" >  $outDir/$longestDNASeqsFromFasta

echo "Beginning comparisons..."
for gene in $(cat $outDir/$allGeneIds)
do
  #echo "*************PROCESSING GENE $gene"
  grep "^>$gene" $outDir/$transdecoderFile | sed 's/>//g' > $outDir/pep.txt
  
  longestLen=0
  longestOrfSeq=""
  
  # see here for explanation of why I'm iterating through the file lines this way: 
  # https://stackoverflow.com/questions/1521462/looping-through-the-content-of-a-file-in-bash
  while IFS="" read -r line || [ -n "$line" ]
  do
    len=`echo "$line" | grep -oE "len:[0-9]*" | sed 's/len://' `
    #echo "*****len:$len,  PEP LINE: $line"
    
    if [ $len -gt $longestLen ]; then
      longestLen=$len
      # save the seq id. The ^ inside the [ ] means match any character except what's inside the [ ].
      longestOrfSeq=`echo "$line" | grep -oE "^[^\.]*" `
      #echo "+++++++ -gt longestLen:$longestLen,  longestOrfSeq: $longestOrfSeq"
    elif [ $len -eq $longestLen ]; then
      
      # get the length of the DNA sequence of the existing longestOrfSeq from the fasta file
      oldDNAlen=`grep -oE "$longestOrfSeq len=[0-9]*" $outDir/$trinityFile | sed "s/$longestOrfSeq len=//" `
      
      # get the length of the DNA sequence of the contending longestOrfSeq from the fasta file
      newSeq=`echo "$line" | grep -oE "^[^\.]*" `
      newDNAlen=`grep -oE "$newSeq len=[0-9]*" $outDir/$trinityFile | sed "s/$newSeq len=//" `
      
      if [ $newDNAlen -gt $oldDNAlen ]; then
        # save the seq id
        longestOrfSeq=$newSeq
        #echo "+++++++ -eq longestLen:$longestLen,  longestOrfSeq: $longestOrfSeq"
      fi
    fi 
    
  done < $outDir/pep.txt  # end inner while loop
  
  
  # ------------
  # if for some reason there's no result from the .pep file, just compare sequence lengths in the fasta file:
  if [ -z $longestOrfSeq ]; then 
    longestLen=0
    grep "^>$gene" $outDir/$trinityFile | sed 's/>//g' > $outDir/fa.txt
    
    while IFS="" read -r line || [ -n "$line" ]
    do
      len=`echo "$line" | grep -oE "len=[0-9]*" | sed 's/len=//' `
      #echo "*****len=$len,  FASTA LINE: $line"
      #echo "grep \">$gene\" ../TransDecoderOutput/WWP18-nt675,367-Trinity.fasta.transdecoder.pep"
      
      if [ $len -gt $longestLen ]; then
        longestLen=$len
        # save the seq id
        longestOrfSeq=`echo "$line" | grep -oE "^\S*" `
        #echo "+++++++ -gt longestLen=$longestLen,  longestFastaSeq: $longestOrfSeq"
      fi 
      
    done < $outDir/fa.txt  # end inner while loop    
  
    echo $longestOrfSeq >> $outDir/$longestDNASeqsFromFasta
  else 
    echo $longestOrfSeq >> $outDir/$longestOrfSeqsFromPep
  fi
  
  
  #echo "------- DONE GENE longestLen:$longestLen,  longestOrfSeq: $longestOrfSeq"
  echo $longestOrfSeq >> $outDir/$longestOrfSeqsAll
  
done  # end outer for-loop

echo "Finished selecting seq Ids based on longest ORF or longest DNA sequence (if no ORF result in .pep file)."
echo "Deleting temporary files."
rm $outDir/pep.txt $outDir/fa.txt


################################################################
################################################################

echo "Calling $rScriptToSelectSeqsFromFasta to extract the selected sequences from the original fasta file."
echo "==========================="

Rscript $outDir/$rScriptToSelectSeqsFromFasta $outDir $longestOrfSeqsAll $trinityFile $finalOutputFasta

echo "==========================="
echo "finished"


