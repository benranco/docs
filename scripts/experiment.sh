#!/bin/bash

input=$1 
single=$2
remove_indels=$3

export PERL5LIB=`pwd`"/tools/vcftools/src/perl"

echo "running post_process.sh"

filterInputFilename="$1.vcf"

if [[ $remove_indels == 1 ]]
then

  echo "removing indels"
  if [[ $single == 1 ]]
  then
    echo "   is single"
  fi

  if [[ $single == 2 ]]
  then
    echo "   is pooled"
  fi

  filterInputFilename="$1.recode.vcf"
fi

echo "filterInputFilename : "$filterInputFilename

echo "filtering using cutoffs"
if [[ $single == 1 ]]
then
  echo "   filtering for single"
fi

if [[ $single == 2 ]]
then
   echo "   filtering for pooled"
fi

echo "tab conversion"

if [[ $single == 1 ]]
then
  echo "   creating tab for single"
fi

if [[ $single == 2 ]]
then
  echo "   creating tab for pooled"
fi


