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
outDir="./test3"

# input file paths, all thes must e relative to outDir:
#transdecoderFile="../../TransDecoderOutput/WWP18-nt675,367-Trinity.fasta.transdecoder.pep"
#trinityFile="../../trin-assemblies-renamedSeqIds/WWP18-nt675,367-Trinity.fasta"

transdecoderFile="../TransDecoderOutput/WWP18-nt675,367-Trinity.fasta.transdecoder.pep"
trinityFile="../trin-assemblies-renamedSeqIds/WWP18-nt675,367-Trinity.fasta"

# this will be used to prefix all files created by this script except the finalOutputFasta file. 
outFilePrefix="test3-wwp-"

# Final output fasta file. If this file already exists, it'll be overwritten.
finalOutputFasta="selected.fasta"  

# The R script to call that does all the comparisons and creates the final output.
rScriptToChooseLongestOrfs="../chooseLongestOrfs_fasta.R"


##########################################################
# Output files created by this script:

allGeneIds="$outFilePrefix"allGeneIds.txt
geneIds100="$outFilePrefix"geneIds-100.txt  # only used for testing with smaller sample sets 

fastaData="$outFilePrefix"fastaData.csv
pepData="$outFilePrefix"pepData.csv



longestOrfSeqsAll="$outFilePrefix"longestOrfSeqs-all.txt
longestOrfSeqsFromPep="$outFilePrefix"longestOrfSeqs-fromPep.txt
longestDNASeqsFromFasta="$outFilePrefix"longestDNASeqs-fromFasta.txt

##########################################################
# Execution code. You shouldn't need to edit below this point.

echo "transdecoderFile: "$outDir/$transdecoderFile
echo "trinityFile: "$outDir/$trinityFile
echo "Location of the R script to extract the sequences with longest orf or seq from the trinityFile: "$outDir/$rScriptToChooseLongestOrfs
echo "Final output fasta file name: "$outDir/$finalOutputFasta

# remove the seqId portion of the seq Ids to leave only the gene ids, and keep only one occurene of each gene Id. The sort command is necessary before uniq because uniq only removes contiguous duplicate lines:
echo "Getting all unique gene ids."
grep -oE "^>\S*" $outDir/$trinityFile | sed 's/>\([^ ]*\)_i[0-9]*$/\1/'  | sort | uniq > $outDir/$allGeneIds

# For testing with smaller sample sets, uncomment the below line, and 
# replace $allGeneIds with $geneIds100 in the below Rscript call.
#head -n 1000 $outDir/$allGeneIds > $outDir/$geneIds100

echo "geneId,seqId,seqLength" > $outDir/$fastaData
# extract only all header lines from a fasta file, but only containing: geneId geneId+seqId seqLenVal
grep "^>" $outDir/$trinityFile | sed "s/>\([^ ]*\)\(_i[0-9]*\) \(len=\)\([0-9]*\).*/\1,\1\2,\4/" >> $outDir/$fastaData

echo "geneId,seqId,seq_p,orfLength" > $outDir/$pepData
# extract only all header lines from a fasta file, but only containing: geneId geneId+seqId .p# orfLenVal
grep "^>" $outDir/$transdecoderFile | sed "s/>\([^ ]*\)\(_i[0-9]*\)\([^ ]*\) .*len:\([0-9]*\).*/\1,\1\2,\3,\4/" >> $outDir/$pepData


################################################################
################################################################

echo "Calling $rScriptToChooseLongestOrfs to extract the selected sequences from the original fasta file."
echo "==========================="

Rscript $outDir/$rScriptToChooseLongestOrfs $outDir $allGeneIds $fastaData $pepData $trinityFile $finalOutputFasta $outFilePrefix

echo "==========================="
echo "finished"


