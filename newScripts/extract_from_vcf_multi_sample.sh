#!/bin/bash

naValue="-"

path="."
inputVcfFile="LP36_pipe1.vcf"

outputCsvFileFull="LP36_pipe1_full.csv"
outputCsvFileMin1="LP36_pipe1_minimal1.csv"
outputCsvFileMin2="LP36_pipe1_minimal2.csv"

echo "Running extract_from_vcf_multi_sample.sh."
echo "Processing "$path/$inputVcfFile


> $path/$outputCsvFileFull
> $path/$outputCsvFileMin1
> $path/$outputCsvFileMin2


cat $path/$inputVcfFile | awk -F '\t' -v naVal="$naValue" -v outFileFull="$path/$outputCsvFileFull" -v outFileMin1="$path/$outputCsvFileMin1" -v outFileMin2="$path/$outputCsvFileMin2" '{ 
if ($1 ~ /^##/ ) { 
  next; 
}
else if ($1 ~ /^#CHROM/ ) { 
  # write the header line
  headerFull=substr($1,2,5)","$2","$4","$5
  headerMin1=headerFull
  headerMin2=headerFull
  # In awk, the NF variable records the number of fields in the current record.
  # The sample records begin at the tenth tab-delimited field.
  for (i = 10; i <= NF; i++) {
    headerFull=headerFull","$i"_GT,"$i"_DP,"$i"_RO,"$i"_AO,"$i"_ratio"
    headerMin1=headerMin1","$i"_GT,"$i"_RO,"$i"_AO"
    headerMin2=headerMin2","$i"_GT,"$i"_ratio"
  }
  print headerFull >> outFileFull
  print headerMin1 >> outFileMin1
  print headerMin2 >> outFileMin2
} 
else if ($1 ~ /^[^#]/ ) { 
  # write the data line
  
  ref=$4
  altFull=$5  # there can be more than one alt value
  alt=altFull
  
  # if altFull is a comma-separated list of more than one value
  if (altFull ~ /.*,.*/ ) {
    numAlts=split(altFull,altArray,",")
    alt=altArray[1] # the alt will always be the first alt, because even if some samples use the second alt, the first might still be used in other samples
  }
  
  lineFull="\""$1"\","$2","ref","alt
  lineMin1=lineFull
  lineMin2=lineFull

  for (i = 10; i <= NF; i++) {
    numSampleFields=split($i,sampleData,":")
    
    if (numSampleFields == 8) {
      # this is a proper sample field containing data
    
      # The sample data is formatted like so: 
      # GT :DP:DPR:RO:QR:AO:QA:GL	
      # 0/0:2 :2,0:2 :77:0 :0 :0,-0.60206,-0.189968

      split(sampleData[1],gtBoth, /\/|\|/ ); # split the genotype gt field by either / or | 
      
      
      if (gtBoth[1] == 0) {
        gtFirst=ref
      }
      else if (gtBoth[1] == 1) {
        gtFirst=alt
      }
      else if (numAlts >= gtBoth[1]) {
        gtFirst=altArray[ gtBoth[1] ]
      }

      if (gtBoth[2] == 0) {
        gtSecond=ref
      }
      else if (gtBoth[2] == 1) {
        gtSecond=alt
      }
      else if (numAlts >= gtBoth[2]) {
        gtSecond=altArray[ gtBoth[2] ]
      }
      
      gtLetters=gtFirst"/"gtSecond
      dp=sampleData[2]
      ro=sampleData[4]
      aoFull=sampleData[6]
      ao=aoFull

      # if there are more than one ao values (because there are more than one alts)
      if (aoFull ~ /.*,.*/ ) {
        numAOs=split(aoFull,aoArray,",")
        ao=aoArray[1]
      } 
      
      # In the rare case when a second or third alt is used in the genotype and the first alt is not in the genotype, use its AO value instead of the first alt AO value 
      if (gtBoth[2] >= 2 && gtBoth[1] != 1 && numAOs >= gtBoth[2]) {
        ao=aoArray[ gtBoth[2] ]
      }
      # intentionally not using an else if here
      if (gtBoth[1] >= 2 && gtBoth[2] != 1 && numAOs >= gtBoth[1]) {
        ao=aoArray[ gtBoth[1] ]
      }
      
      divisor= ro + ao
      if (divisor > 0) {
        ratio = ro / (ro+ao)
      }
      else {
        ratio = 0
      }
      
    }
    else {
      # this sample field contains no data
      gtLetters=naVal
      dp=naVal
      ro=naVal
      ao=naVal
      ratio=naVal
    }
    
      # this looks weird but the quotes are around the commas
      lineFull=lineFull","gtLetters","dp","ro","ao","ratio
      lineMin1=lineMin1","gtLetters","ro","ao
      lineMin2=lineMin2","gtLetters","ratio

  } # end for-loop

  
  print lineFull >> outFileFull
  print lineMin1 >> outFileMin1
  print lineMin2 >> outFileMin2
}
}'


