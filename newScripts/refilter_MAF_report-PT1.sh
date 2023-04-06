#!/bin/bash

reportNameBase="MAF_cutoff_report.csv"
prefix="pipeAll_"
tempPostfix=".tmp.noIndels"
finalPostfix=".noIndels.max4NAs"

# remove lines with indels in the REF
grep -v -E '".*",[0-9]*,"[ACTG][ACTG][ACTG]*",' $prefix$reportNameBase >> $prefix$reportNameBase$tempPostfix

# remove lines with at least 5 NA occurences in the sample columns
grep -v -E '.*,.*,.*,NA.*,NA.*,NA.*,NA.*,NA' $prefix$reportNameBase$tempPostfix >> $prefix$reportNameBase$finalPostfix



