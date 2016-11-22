#!/bin/bash

# ############################################################
# Set the following input parameters:

# make sure to end with a "/"
lepMap2Bin="/home/benrancourt/Desktop/LepMAP2/binary/bin/"

# 1==do the Filtering module, 0==skip the Filtering module
doFiltering=0

# make sure to end with a "/"
dataLocationWithFinalSlash="/home/benrancourt/Downloads/"

inputFileNameMainPart="R28-S91-M28K-Ben"
inputFileNameDotSuffix=".linkage"


# TODO: set lodlimit and sizeLimit, etc here rather than hardcoded into the command strings

# ############################################################
# ############################################################

logFileNameSuffix="-lepMAP2-log.txt"

# output file naming conventions
filteredOutputSuffix="-filtered.linkage"

mainInputFileSuffix=$inputFileNameDotSuffix
separateChromosomesOutputSuffix="-map.txt"
joinSinglesOutputSuffix="-map_js.txt"
orderMarkersOutputSuffix="-map_js-chr1.SA.txt"

orderMarkersOutputSuffixPt1="-map_js-chr"
orderMarkersOutputSuffixPt2=".SA.txt"
orderMarkersLogFileNameSuffixPt1="-orderMarkers-chr"
orderMarkersLogFileNameSuffixPt2="-lepMAP2-log.txt"

if [[ $doFiltering -eq 1 ]]
then
  mainInputFileSuffix=$filteredOutputSuffix
  separateChromosomesOutputSuffix="-filtered-map.txt"
  joinSinglesOutputSuffix="-filtered-map_js.txt"
  orderMarkersOutputSuffix="-filtered-map_js-chr1.SA.txt"
  orderMarkersMFOutputSuffix="-filtered-map_js-chr1.MF.txt"
fi

dl=$dataLocationWithFinalSlash

# ############################################################
#build the commands

#default dataTolerance is 0.01
filterCommand="java -cp "$lepMap2Bin" Filtering data="$dl$inputFileNameMainPart$inputFileNameDotSuffix

separateChromosomesCommand="java -cp "$lepMap2Bin" SeparateChromosomes data="$dl$inputFileNameMainPart$mainInputFileSuffix" lodLimit=6 sizeLimit=10"

joinSinglesCommand="java -cp "$lepMap2Bin" JoinSingles "$dl$inputFileNameMainPart$separateChromosomesOutputSuffix" data="$dl$inputFileNameMainPart$mainInputFileSuffix" lodLimit=3"

orderMarkersCommand="java -cp "$lepMap2Bin" OrderMarkers data="$dl$inputFileNameMainPart$mainInputFileSuffix" map="$dl$inputFileNameMainPart$joinSinglesOutputSuffix" alpha=0.1 chromosome=1 polishWindow=100 filterWindow=10  sexAveraged=1"

orderMarkersCommandChooseLG="java -cp "$lepMap2Bin" OrderMarkers data="$dl$inputFileNameMainPart$mainInputFileSuffix" map="$dl$inputFileNameMainPart$joinSinglesOutputSuffix" alpha=0.1 polishWindow=100 filterWindow=10  sexAveraged=1 chromosome="

#-----
#orderMarkersCommand1="java -cp "$lepMap2Bin" OrderMarkers data="$dl$inputFileNameMainPart$mainInputFileSuffix" map="$dl$inputFileNameMainPart$joinSinglesOutputSuffix" alpha=0.1 chromosome=1 sexAveraged=1 improveOrder=0 polishWindow=100 filterWindow=10"

#orderMarkersCommand2="java -cp "$lepMap2Bin" OrderMarkers data="$dl$inputFileNameMainPart$mainInputFileSuffix" map="$dl$inputFileNameMainPart$joinSinglesOutputSuffix" alpha=0.1 chromosome=1 sexAveraged=1 improveOrder=0"

#orderMarkersCommand3a="java -cp "$lepMap2Bin" OrderMarkers data="$dl$inputFileNameMainPart$mainInputFileSuffix" map="$dl$inputFileNameMainPart$joinSinglesOutputSuffix" alpha=0.1 chromosome=1 improveOrder=0 polishWindow=100 filterWindow=10"

#orderMarkersCommand3b="java -cp "$lepMap2Bin" OrderMarkers data="$dl$inputFileNameMainPart$mainInputFileSuffix" map="$dl$inputFileNameMainPart$joinSinglesOutputSuffix" alpha=0.1 sexAveraged=1 evaluateOrder="$dl$inputFileNameMainPart$orderMarkersMFOutputSuffix" improveOrder=0"

#orderMarkersCommand4a="java -cp "$lepMap2Bin" OrderMarkers data="$dl$inputFileNameMainPart$mainInputFileSuffix" map="$dl$inputFileNameMainPart$joinSinglesOutputSuffix" alpha=0.1 chromosome=1 polishWindow=100 filterWindow=10"

#orderMarkersCommand4b="java -cp "$lepMap2Bin" OrderMarkers data="$dl$inputFileNameMainPart$mainInputFileSuffix" map="$dl$inputFileNameMainPart$joinSinglesOutputSuffix" alpha=0.1 sexAveraged=1 evaluateOrder="$dl$inputFileNameMainPart$orderMarkersMFOutputSuffix" improveOrder=0"

#orderMarkersCommand5a="java -cp "$lepMap2Bin" OrderMarkers data="$dl$inputFileNameMainPart$mainInputFileSuffix" map="$dl$inputFileNameMainPart$joinSinglesOutputSuffix" alpha=0.1 chromosome=1 polishWindow=100 filterWindow=10 improveOrder=0"

#orderMarkersCommand5b="java -cp "$lepMap2Bin" OrderMarkers data="$dl$inputFileNameMainPart$mainInputFileSuffix" map="$dl$inputFileNameMainPart$joinSinglesOutputSuffix" alpha=0.1 sexAveraged=1 evaluateOrder="$dl$inputFileNameMainPart$orderMarkersMFOutputSuffix" improveOrder=0"

#orderMarkersCommand6="java -cp "$lepMap2Bin" OrderMarkers data="$dl$inputFileNameMainPart$mainInputFileSuffix" map="$dl$inputFileNameMainPart$joinSinglesOutputSuffix" alpha=0.1 chromosome=1 polishWindow=100 filterWindow=10  sexAveraged=1 improveOrder=0"

#orderMarkersCommand7="java -cp "$lepMap2Bin" OrderMarkers data="$dl$inputFileNameMainPart$mainInputFileSuffix" map="$dl$inputFileNameMainPart$joinSinglesOutputSuffix" alpha=0.1 chromosome=1 polishWindow=100 filterWindow=10  sexAveraged=1"



# ############################################################
function getNumLinkageGroups {
  # parameter $1 should be the name the map file we will be checking to determin number of linkage groups

  # number of linkage groups = the number of unique lines (minus the header) in the file
  numLinkageGroups=`sort -u $dl$1 | wc -l` 
  numLinkageGroups=$((numLinkageGroups-1)) # - 1 because of the header line
}


# ############################################################
function countMarkersInLinkageGroups {
  mapFileName=$1 # $1 retrieves the first input parameter to this function

  ### Count the number of markers in each Linkage Group
  echo " " | tee -a $dl$inputFileNameMainPart$logFileNameSuffix  
  echo "===========================================================================" >> $dl$inputFileNameMainPart$logFileNameSuffix
  echo "===========================================================================" >> $dl$inputFileNameMainPart$logFileNameSuffix
  date | tee -a $dl$inputFileNameMainPart$logFileNameSuffix

  getNumLinkageGroups $mapFileName

  echo "Number of Linkage Groups in "$mapFileName": "$numLinkageGroups" (might be only "$((numLinkageGroups-1))" if some markers are in 0=no linkage group)" | tee -a $dl$inputFileNameMainPart$logFileNameSuffix  

  for i in `seq 0 $numLinkageGroups`
  do 
    count=`grep -c "^$i$" $dataLocationWithFinalSlash$mapFileName`
    echo "Number of markers in LG "$i": "$count | tee -a $dl$inputFileNameMainPart$logFileNameSuffix  
  done

}




# ############################################################
#execute:
 
echo "Beginning processing $inputFileNameMainPart$inputFileNameDotSuffix with LepMAP2 modules." | tee $dl$inputFileNameMainPart$logFileNameSuffix
echo " " >> $dl$inputFileNameMainPart$logFileNameSuffix

### Filtering
if [[ $doFiltering -eq 1 ]]
then
  echo " " | tee -a $dl$inputFileNameMainPart$logFileNameSuffix  
  echo "===========================================================================" >> $dl$inputFileNameMainPart$logFileNameSuffix
  echo "===========================================================================" >> $dl$inputFileNameMainPart$logFileNameSuffix
  date | tee -a $dl$inputFileNameMainPart$logFileNameSuffix  
  echo "Beginning LepMAP2 Filtering module. Output filename: $inputFileNameMainPart$filteredOutputSuffix." | tee -a $dl$inputFileNameMainPart$logFileNameSuffix
  echo " " >> $dl$inputFileNameMainPart$logFileNameSuffix

  $filterCommand > $dl$inputFileNameMainPart$filteredOutputSuffix 2>> $dl$inputFileNameMainPart$logFileNameSuffix
else
  echo " " | tee -a $dl$inputFileNameMainPart$logFileNameSuffix
  echo "doFiltering==0, therefore skipping LepMAP2 Filtering module." | tee -a $dl$inputFileNameMainPart$logFileNameSuffix
fi

### SeparateChromosomes
echo " " | tee -a $dl$inputFileNameMainPart$logFileNameSuffix  
echo "===========================================================================" >> $dl$inputFileNameMainPart$logFileNameSuffix
echo "===========================================================================" >> $dl$inputFileNameMainPart$logFileNameSuffix
date | tee -a $dl$inputFileNameMainPart$logFileNameSuffix  
echo "Beginning LepMAP2 SeparateChromosomes module. Output filename: $inputFileNameMainPart$separateChromosomesOutputSuffix." | tee -a $dl$inputFileNameMainPart$logFileNameSuffix
echo " " >> $dl$inputFileNameMainPart$logFileNameSuffix

$separateChromosomesCommand > $dl$inputFileNameMainPart$separateChromosomesOutputSuffix 2>> $dl$inputFileNameMainPart$logFileNameSuffix

### Count the number of markers in each Linkage Group in the SeparateChromosomes map file
countMarkersInLinkageGroups $inputFileNameMainPart$separateChromosomesOutputSuffix

### JoinSingles
echo " " | tee -a $dl$inputFileNameMainPart$logFileNameSuffix  
echo "===========================================================================" >> $dl$inputFileNameMainPart$logFileNameSuffix
echo "===========================================================================" >> $dl$inputFileNameMainPart$logFileNameSuffix
date | tee -a $dl$inputFileNameMainPart$logFileNameSuffix  
echo "Beginning LepMAP2 JoinSingles module. Output filename: $inputFileNameMainPart$joinSinglesOutputSuffix." | tee -a $dl$inputFileNameMainPart$logFileNameSuffix
echo " " >> $dl$inputFileNameMainPart$logFileNameSuffix

$joinSinglesCommand > $dl$inputFileNameMainPart$joinSinglesOutputSuffix 2>> $dl$inputFileNameMainPart$logFileNameSuffix

### Count the number of markers in each Linkage Group after JoinSingles
countMarkersInLinkageGroups $inputFileNameMainPart$joinSinglesOutputSuffix

### OrderMarkers MF
#echo " " | tee -a $dl$inputFileNameMainPart$logFileNameSuffix  
#echo "===========================================================================" >> $dl$inputFileNameMainPart$logFileNameSuffix
#echo "===========================================================================" >> $dl$inputFileNameMainPart$logFileNameSuffix
#date | tee -a $dl$inputFileNameMainPart$logFileNameSuffix  
#echo "Beginning LepMAP2 OrderMarkers MF module. Output filename: $inputFileNameMainPart$orderMarkersMFOutputSuffix." | tee -a $dl$inputFileNameMainPart$logFileNameSuffix
#echo " " >> $dl$inputFileNameMainPart$logFileNameSuffix

#$orderMarkersCommand5a > $dl$inputFileNameMainPart$orderMarkersMFOutputSuffix 2>> $dl$inputFileNameMainPart$logFileNameSuffix


### OrderMarkers SA
echo " " | tee -a $dl$inputFileNameMainPart$logFileNameSuffix  
echo "===========================================================================" >> $dl$inputFileNameMainPart$logFileNameSuffix
echo "===========================================================================" >> $dl$inputFileNameMainPart$logFileNameSuffix
#date | tee -a $dl$inputFileNameMainPart$logFileNameSuffix  
#echo "Beginning LepMAP2 OrderMarkers SA module. Output filename: $inputFileNameMainPart$orderMarkersOutputSuffix." | tee -a $dl$inputFileNameMainPart$logFileNameSuffix
#echo " " >> $dl$inputFileNameMainPart$logFileNameSuffix

#$orderMarkersCommand > $dl$inputFileNameMainPart$orderMarkersOutputSuffix 2>> $dl$inputFileNameMainPart$logFileNameSuffix



### =========================================
date | tee -a $dl$inputFileNameMainPart$logFileNameSuffix  
echo "Beginning LepMAP2 OrderMarkers SA modules. See individual Order Markers log files for each Linkage Group." | tee -a $dl$inputFileNameMainPart$logFileNameSuffix

getNumLinkageGroups $inputFileNameMainPart$joinSinglesOutputSuffix
counter=0
ncore="$(( ($(grep -c ^processor /proc/cpuinfo) -4)/2 ))"

echo "running parallel processing"

for i in `seq 1 $numLinkageGroups`
do 
  if [ $counter -lt $ncore ]; then	
    date > -a $dl$inputFileNameMainPart$orderMarkersLogFileNameSuffixPt1$i$orderMarkersLogFileNameSuffixPt2
    echo "Beginning LepMAP2 OrderMarkers SA module. Output filename: $dl$inputFileNameMainPart$orderMarkersOutputSuffixPt1$i$orderMarkersOutputSuffixPt2." | tee -a $dl$inputFileNameMainPart$orderMarkersLogFileNameSuffixPt1$i$orderMarkersLogFileNameSuffixPt2
    echo " " >> $dl$inputFileNameMainPart$orderMarkersLogFileNameSuffixPt1$i$orderMarkersLogFileNameSuffixPt2

    $orderMarkersCommandChooseLG$i > $dl$inputFileNameMainPart$orderMarkersOutputSuffixPt1$i$orderMarkersOutputSuffixPt2 2>> $dl$inputFileNameMainPart$orderMarkersLogFileNameSuffixPt1$i$orderMarkersLogFileNameSuffixPt2 &
    counter=$(( $counter + 1 ))
  else
    wait
    #echo "wait reset"
    counter=0
  fi
done
wait

### =========================================



### Done
echo " " | tee -a $dl$inputFileNameMainPart$logFileNameSuffix  
echo "===========================================================================" >> $dl$inputFileNameMainPart$logFileNameSuffix
echo "===========================================================================" >> $dl$inputFileNameMainPart$logFileNameSuffix
date | tee -a $dl$inputFileNameMainPart$logFileNameSuffix  

echo "Finished processing." | tee -a $dl$inputFileNameMainPart$logFileNameSuffix

echo  -e "\033[33;5;7mPipeline $name has finished\033[0m"


 

