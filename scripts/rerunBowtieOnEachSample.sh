#!/bin/bash

#############################################################################################
# This script should be placed in the same folder as the start.sh script for the SNPpipeline.
# WARNING: Running this script will overwrite and delete any .sam files in the dataTemp folder,
# and will also overwrite any *_align.log files in the log folder.
# If this files are important, make backups of them first.
#############################################################################################

#(1) single, (2) pooled
single=1

# paired files (e.g. R1 & R2) vs unpaired files: (1) paired, (0) unpaired
paired=1

ncore="$(( ($(grep -c ^processor /proc/cpuinfo) -4)))"

# $1=sample, $2=ncore, $3=paired




if [[ single -eq 1 ]]
then
    echo "running single"
    if [[ paired -eq 0 ]]; then
	    for datapoint in $(ls "./data"); do
	        bash ./scripts/align.sh $datapoint $ncore $paired > ./logs/single/$1"_align.log" 2>&1 
          rm -f ./dataTemp/single/$datapoint.sam
	    done
    elif [[ paired -eq 1 ]]; then
	    for datapoint in $(ls "./data" | rev | cut -c 13- | rev | uniq)
	    do
	        bash ./scripts/align.sh $datapoint $ncore $paired > ./logs/single/$1"_align.log" 2>&1 
          rm -f ./dataTemp/single/$datapoint.sam
	        #sync
	        #echo 1 > /proc/sys/vm/drop_caches
	    done
    fi
elif [[ single -eq 2 ]]
then
    echo "running pooled"
    bash ./scripts/align.sh $datapoint $ncore $paired > ./logs/pooled/$1"_align.log" 2>&1 
    rm -f ./dataTemp/pooled/$datapoint.sam
fi


