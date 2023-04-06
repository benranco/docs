#!/bin/bash

if [ $# -lt 3 ];
then
    echo "usage: $0 [regions file] [ncpus] [freebayes arguments]"
    echo
    echo "Run freebayes in parallel over regions listed in regions file, using ncpus processors."
    echo "Will merge and sort output, producing a uniform VCF stream on stdout.  Flags to freebayes"
    echo "which would write to e.g. a particular file will obviously cause problms, so caution is"
    echo "encouraged when using this script."
    echo
    echo "examples:"
    echo
    echo "Run freebayes in parallel on 100000bp chunks of the ref (fasta_generate_regions.py is also"
    echo "located in the scripts/ directory in the freebayes distribution).  Use 36 threads."
    echo
    echo "    freebayes-parallel <(fasta_generate_regions.py ref.fa.fai 100000) 36 -f ref.fa aln.bam >out.vcf"
    echo
    echo "Generate regions that are equal in terms of data content, and thus have lower variance"
    echo "in runtime.  This will yield better resource utilization."
    echo
    echo "    bamtools coverage -in aln.bam | coverage_to_regions.py ref.fa 500 >ref.fa.500.regions"
    echo "    freebayes-parallel ref.fa.500.regions 36 -f ref.fa aln.bam >out.vcf"
    echo
    exit
fi

regionsfile=$1
shift
ncpus=$1
shift

curTime=`date +%Y%m%d%H%M%S`

tmpSortDir=tmpSortDir_$curTime
vcfSortPath="tools/vcftools/src/perl/vcf-sort"

tmpExtractedVcf=freebayesParallelBen_tmpExtractedVcf_$curTime.vcf
tmpSortedVcf=freebayesParallelBen_tmpSortedVcf_$curTime.vcf

maxOpenFiles=`ulimit -Hn`
# set this sessions "Soft" max open files to its currently allowed max open files 
ulimit -Sn $maxOpenFiles

command=("tools/freebayes/bin/freebayes" "$@")

(
#$command | head -100 | grep "^#" # generate header
# iterate over regions using gnu parallel to dispatch jobs
cat "$regionsfile" | parallel -k -j "$ncpus" "${command[@]}" --region {}
) | awk -F '\t' -v isHeader=1 '{ 
if ($1 ~ /^##/ && isHeader == 1) { print $0; }
else if ($1 ~ /^#[^#]/ && isHeader == 1) { print $0; isHeader=0; } 
else if ($1 ~ /^[^#]/ && $1 !~ /^parallel/ ) { print $0; }
}' 2>"$tmpExtractedVcf".log 1> $tmpExtractedVcf

# The above awk command does this (saving the first occurence of the header content but ignoring all other re-occurences of it): 
# isHeader=1
# if (line begins with ## && isHeader) print line
# if (line begins with # && isHeader) print line; isHeader=0
# if (line doesn't begin with # and doesn't begin with a parallel warning) print line


# sort the vcf data
mkdir $tmpSortDir
cat $tmpExtractedVcf | $vcfSortPath -c -t ./$tmpSortDir 2>"$tmpSortedVcf".log 1> $tmpSortedVcf
rm -r $tmpSortDir

#rm $tmpExtractedVcf
#rm "$tmpExtractedVcf".log

# Skip all lines that have the same CHROM, POS, REF and ALT as the line above them (assumes the input is a sorted VCF file).
# Note: if the ALT contains multiple comma-separated values, the full string, in the same order, must match in order to skip the line.
cat $tmpSortedVcf | awk -F '\t' -v chrom="" -v pos=0 -v ref="" -v alt="" '{ 
if (chrom != "" && $1 == chrom && $2 == pos && $4 == ref && $5 == alt) { next; }
else { chrom=$1; pos=$2; ref=$4; alt=$5; print $0; } 
}' 
#> $tmpOutputFileNameBase
# The output of the above command is what will be streamed to stdout.

#rm $tmpSortedVcf
#rm "$tmpSortedVcf".log

