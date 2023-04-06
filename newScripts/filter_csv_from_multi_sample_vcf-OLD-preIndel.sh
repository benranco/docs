#!/bin/bash

# Filter pseudocode: 
#
# specify array of all genotype column names 
# specify array of RS genotype column names
# specify array of SS genotype column names
# numMissingRSData=2
# numMissingSSData=2
# mafCutoff=0.01
# for each row
#   if header row, get col nums of RS and SS samples
#   else
#     if REF length > 1, remove row
#     for each sample
#       if RS and missing, ++rsMissingTally
#       if SS and missing, ++ssMissingTally
#     if rsMissingTally + ssMissingTally > allowXmissingRSplusSSdata, remove row
#     calculate MAF
#     if maf < mafCutoff, remove row
# done
#     
# Calculate MAF pseudocode: 
# 
# specify array of Genotype column names
# for each row
#   tally A, C, T, G 
#   select max, second_max
#   MAF = second_max / second_max + max


######################################################
# Input parameters: 

path="."
inputCsvFile="LP36_pipe1_full-first1000.csv"

outputCsvFileNameBase="yLP36_pipe1-first1000" 

allowXmissingRSplusSSdata=4

mafCutoff=0.1

naValue="-"

refColName="REF"

# the format of these is just a string with names separated by commas, they'll be converted into arrays later:
genotypeColNames="C10_GT,C11_GT,C12_GT,C13_GT,C1_GT,C3_GT,C4_GT,C5_GT,C6_GT,C7_GT,C8_GT,C9_GT,LP-C2_GT,RN1_GT,RS10_GT,RS1_GT,RS2_GT,RS3_GT,RS4_GT,RS5_GT,RS6_GT,RS7_GT,RS8_GT,RS9_GT,SN1_GT,SN2_GT,SS10_GT,SS1_GT,SS2_GT,SS3_GT,SS4_GT,SS5_GT,SS6_GT,SS7_GT,SS8_GT,SS9_GT"
RSgenotypeColNames="RS10_GT,RS1_GT,RS2_GT,RS3_GT,RS4_GT,RS5_GT,RS6_GT,RS7_GT,RS8_GT,RS9_GT"
SSgenotypeColNames="SS10_GT,SS1_GT,SS2_GT,SS3_GT,SS4_GT,SS5_GT,SS6_GT,SS7_GT,SS8_GT,SS9_GT"


######################################################
# Execution code below this point. 

echo "Running filter_csv_from_multi_sample_vcf.sh."
echo "Processing "$path/$inputCsvFile



headerLine=`head -n 1 $path/$inputCsvFile`

cat $path/$inputCsvFile | awk -F ',' -v header="$headerLine" -v naVal="$naValue" -v refCol=$refColName -v allowMissing=$allowXmissingRSplusSSdata -v mafCut=$mafCutoff -v allGenotypes="$genotypeColNames" -v rsGenotypes="$RSgenotypeColNames" -v ssGenotypes="$SSgenotypeColNames" -v outFileBase="$path/$outputCsvFileNameBase" '

##################
function printDictionaryElements(dict, 
                                 #### local-scope variable below here
                                 key) {
  for (key in dict) { print key ": " dict[key] }
}

##################
function indexOf(element, dict, 
                 #### local-scope variable below here
                 key) {
  for (key in dict) {
    if (element == dict[key]) {
      return key
    }
  }
  return -1
}

##################
function calculateMAF(totalA, totalC, totalG, totalT, mafStatsArr,
                      ####  local-scope variables below here
                      max, secondMax) {

  # Calculate MAF pseudocode: 
  # 
  #   tally A, C, T, G (passed in as arguments)
  #   select max, secondMax
  #   MAF = secondMax / (secondMax + max)

  if (totalA >= totalC) {
    max=totalA
    secondMax=totalC
  }
  else {
    max=totalC
    secondMax=totalA
  }
  
  if (totalG >= max) {
    secondMax=max
    max=totalG
  }
  else if (totalG > secondMax) {
    secondMax=totalG
  }
  
  if (totalT >= max) {
    secondMax=max
    max=totalT
  }
  else if (totalT > secondMax) {
    secondMax=totalT
  }

  mafStatsArr["max"] = max
  mafStatsArr["second_max"] = secondMax
  mafStatsArr["maf"] = secondMax / (secondMax + max)
  
  return mafStatsArr["maf"]
}

##################
function generateGenoStats(allGenoIndexes, naValue, statsArr,
                           ####  local-scope variables below here
                           totalA, totalC, totalG, totalT, totalEmpty, i, maf) {

    totalA=0
    totalC=0
    totalG=0
    totalT=0
    totalEmpty=0
    
    for (i in allGenoIndexes) {
      # The $i accesses the ith element of the current input line
      if ($i == naValue) {
        totalEmpty++
      }
      else {
      	# use gsub to count the occurences of A, since it returns the number of substitutions
        totalA += gsub(/A/, "A", $i)
        totalC += gsub(/C/, "C", $i)
        totalG += gsub(/G/, "G", $i)
        totalT += gsub(/T/, "T", $i)
      }
    } # end for-loop
    
    maf = calculateMAF(totalA, totalC, totalG, totalT, mafStatsArr)

    statsArr["A"] = totalA
    statsArr["C"] = totalC
    statsArr["G"] = totalG
    statsArr["T"] = totalT
    statsArr["empty"] = totalEmpty
    statsArr["sum"] = totalA + totalC + totalG + totalT
    statsArr["max"] = mafStatsArr["max"]
    statsArr["second_max"] = mafStatsArr["second_max"]
    statsArr["MAF"] = maf
    
    return maf
}

##################
BEGIN { 

  outFileStats = outFileBase"-stats.csv"
  outFileFull = outFileBase"-full.csv"
  outFileMin1 = outFileBase"-minimal1.csv"
  outFileMin2 = outFileBase"-minimal2.csv"
  outFileGt = outFileBase"-genotype.csv"
  outFileRatio = outFileBase"-ratio.csv"
  
  # split saves a dictionary of values and numerical keys into the second input parameter
  numAllGenos=split(allGenotypes,allGenos,",")
  numRsGenos=split(rsGenotypes,rsGenos,",")
  numSsGenos=split(ssGenotypes,ssGenos,",")
  numAllCols=split(header,allCols,",")
  
  print "-------------------------"
  print "Log messages printed to stdout. Data output printed to the provided output file."
  print " "
  print "Number of genotype columns: "numAllGenos
  print "Number of RS genotype columns: "numRsGenos
  print "Number of SS genotype columns: "numSsGenos
  print "Total number of columns: "numAllCols
  print " "
  refIndex=indexOf(refCol,allCols)
  
  print "Index of the ref col: "refIndex
  print "MAF cutoff = "mafCut
  
  # get the column nums of the genotype columns from allCols
  
  for (key in ssGenos) {
    # this creates a new dictionary without values, ssGenoIndexes, whose keys are the indices of allCols whose values correspond to the values of the ssGenos dictionary
    ssGenoIndexes[indexOf(ssGenos[key],allCols)]   
  }
  
  for (key in rsGenos) {
    rsGenoIndexes[indexOf(rsGenos[key],allCols)]   
  }

  for (key in allGenos) {
    allGenoIndexes[indexOf(allGenos[key],allCols)]   
  }
  
#  print "--"
#  printDictionaryElements(ssGenoIndexes)
#  print "--"
#  printDictionaryElements(rsGenoIndexes)
#  print "--"
#  printDictionaryElements(allGenoIndexes)
#  
#  print "--"
#  if (175 in ssGenoIndexes) print "ssGenoIndexes TRUE"
#  if (85 in rsGenoIndexes) print "rsGenoIndexes TRUE"
#  if (30 in allGenoIndexes) print "allGenoIndexes TRUE"


  print " "
  #print "Lines removed listed below: "
  print "-------------------------"
  print "Processing..."
  
  lineNum=0
  numIndels=0
  skippedMissingData=0
  skippedMaf=0
  printedLines=0
}

################## Main
{ 
  lineNum++;
  
  # Pseudocode, for each line:
  #     if header line, print
  #     if REF length > 1, remove row
  #     for each sample
  #       if RS and missing, ++rsMissingTally
  #       if SS and missing, ++ssMissingTally
  #     if rsMissingTally + ssMissingTally > 4, remove row
  #     calculate MAF
  #     if maf < mafCutoff, remove row
  

  # header line
  if (lineNum == 1) {
    print $1","$2","$3","$4",A,C,G,T,empty,sum,max,second_max,MAF" > outFileStats
    print $0 > outFileFull

    # Get the indices of the columns to print for the additional csv output files.
    # This is an example of the input column naming convention for just one sample (beginning in column 5):
    # CHROM,POS,REF,ALT,sample1_GT,sample1_DP,sample1_RO,sample1_AO,sample1_ratio,
     
    outputMin1Indexes="1,2,3,4"
    outputMin2Indexes=outputMin1Indexes
    outputGtIndexes=outputMin1Indexes
    outputRatioIndexes=outputMin1Indexes

    min1Header=$1","$2","$3","$4
    min2Header=min1Header
    gtHeader=min1Header
    ratioHeader=min1Header
    
    for (i = 5; i <= NF; i++) {      
      if ($i ~ /_GT|_RO|_AO$/ ) {
        outputMin1Indexes=outputMin1Indexes","i
        min1Header=min1Header","$i
      }
      if ($i ~ /_GT|_ratio$/ ) {
        outputMin2Indexes=outputMin2Indexes","i
        min2Header=min2Header","$i
      }   
      if ($i ~ /_GT$/ ) {
        outputGtIndexes=outputGtIndexes","i
        gtHeader=gtHeader","$i
      } 
      if ($i ~ /_ratio$/ ) {
        outputRatioIndexes=outputRatioIndexes","i
        ratioHeader=ratioHeader","$i
      }    
    } # end for-loop
    
    print min1Header > outFileMin1
    print min2Header > outFileMin2
    print gtHeader > outFileGt
    print ratioHeader > outFileRatio
        
    # split saves a dictionary of values and numerical keys into the second input parameter
    numMin1Cols=split(outputMin1Indexes,outputMin1Dict,",")
    numMin2Cols=split(outputMin2Indexes,outputMin2Dict,",")
    numGtCols=split(outputGtIndexes,outputGtDict,",")
    numRatioCols=split(outputRatioIndexes,outputRatioDict,",")
    
    for (key in outputMin1Dict) {
      # this creates a new dictionary without values, outputMin1Cols, whose keys are the indices of the columns which should be included in the output
      outputMin1Cols[outputMin1Dict[key]]   
    }    
    for (key in outputMin2Dict) {
      outputMin2Cols[outputMin2Dict[key]]   
    }
    for (key in outputGtDict) {
      outputGtCols[outputGtDict[key]]   
    }
    for (key in outputRatioDict) {
      outputRatioCols[outputRatioDict[key]]   
    }
      
    next
  }
  # exclude all lines with indels
  else if (length($refIndex) > 1) {
    numIndels++
    #print lineNum" : indel, SKIP"
    next  # skip row if it is an indel
  }
  else {
    
    rsMissingTally=0 
    ssMissingTally=0 
    
    for (i in rsGenoIndexes) {
      if ($i == naVal) {
        rsMissingTally++        
      }
    }
    
    for (i in ssGenoIndexes) {
      if ($i == naVal) {
        ssMissingTally++        
      }
    }
    
    if (rsMissingTally + ssMissingTally > allowMissing) {
      skippedMissingData++
      #print lineNum" : too many missing genotypes: rs"rsMissingTally", ss"ssMissingTally", SKIP"
      next  # skip row if too much missing data
    }
    
    # Calculate MAF
        
    maf = generateGenoStats(allGenoIndexes, naVal, statsArr) 
    
    if (maf < mafCut) {
      skippedMaf++
      #print lineNum" : rs"rsMissingTally", ss"ssMissingTally", A="totalA", C="totalC", G="totalG", T="totalT", maf="maf", SKIP"
      next  # skip row if maf is below the maf cutoff
    }
    
    # the row passes all filtering, so print it to the output files in various formats
     
    # print the stats file. This looks weird but the quotes are around the commas
    print $1","$2","$3","$4","statsArr["A"]","statsArr["C"]","statsArr["G"]","statsArr["T"]","statsArr["empty"]","statsArr["sum"]","statsArr["max"]","statsArr["second_max"]","statsArr["MAF"] >> outFileStats
    
    
    print $0 >> outFileFull
    
    min1Line=""
    min2Line=""
    gtLine=""
    ratioLine=""
    
    divider=""
    
    for (i=1; i<=NF; i++) {
      if (i > 1) {
        divider=","
      }
      if (i in outputMin1Cols) {
        min1Line=min1Line divider $i
      }
      if (i in outputMin2Cols) {
        min2Line=min2Line divider $i
      }
      if (i in outputGtCols) {
        gtLine=gtLine divider $i
      }
      if (i in outputRatioCols) {
        ratioLine=ratioLine divider $i
      }
    } # end for-loop
    
    print min1Line >> outFileMin1
    print min2Line >> outFileMin2
    print gtLine >> outFileGt
    print ratioLine >> outFileRatio
    
    printedLines++
  }
    

}
##################
END {

  print "-------------------------"
  print "Finished filtering."
  print " "
  print lineNum-1" total input data lines."
  print printedLines" data lines printed to output."
  print " "
  print numIndels+skippedMissingData+skippedMaf" lines skipped in total."
  print numIndels" lines skipped because of indels."
  print skippedMissingData" lines skipped because of too much missing data."
  print skippedMaf" lines skipped because their MAF was below the MAF cutoff."

  
}'





