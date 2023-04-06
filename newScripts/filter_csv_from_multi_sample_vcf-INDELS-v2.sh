#!/bin/bash

# Filter pseudocode (might be out of date since adding indel capabilities): 
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

outputCsvFileNameBase="yLP36_pipe1-first1000-filtered" 

# 0=keep both indels and non-indels, 1=keep only non-indels, 2=keep only indels, 3=print indel and non-indel output to separate files
indelStatus=3

allowXmissingRSplusSSdata=4

mafCutoff=0.01

naValue="-"

refColName="REF"
altColName="ALT"

# the format of these is just a string with names separated by commas, they'll be converted into arrays later:
genotypeColNames="C10_GT,C11_GT,C12_GT,C13_GT,C1_GT,C3_GT,C4_GT,C5_GT,C6_GT,C7_GT,C8_GT,C9_GT,LP-C2_GT,RN1_GT,RS10_GT,RS1_GT,RS2_GT,RS3_GT,RS4_GT,RS5_GT,RS6_GT,RS7_GT,RS8_GT,RS9_GT,SN1_GT,SN2_GT,SS10_GT,SS1_GT,SS2_GT,SS3_GT,SS4_GT,SS5_GT,SS6_GT,SS7_GT,SS8_GT,SS9_GT"
RSgenotypeColNames="RS10_GT,RS1_GT,RS2_GT,RS3_GT,RS4_GT,RS5_GT,RS6_GT,RS7_GT,RS8_GT,RS9_GT"
SSgenotypeColNames="SS10_GT,SS1_GT,SS2_GT,SS3_GT,SS4_GT,SS5_GT,SS6_GT,SS7_GT,SS8_GT,SS9_GT"

# the number of fields/columns used per sample in the input .csv file
numFieldsPerSample=5

######################################################
# Execution code below this point. 

echo "Running filter_csv_from_multi_sample_vcf.sh."
echo "Processing "$path/$inputCsvFile



# this assumes the header line is the first line in the file 
headerLine=`head -n 1 $path/$inputCsvFile` 

cat $path/$inputCsvFile | awk -F ',' -v header="$headerLine" -v indelStatus="$indelStatus" -v naVal="$naValue" -v refCol=$refColName -v altCol=$altColName -v allowMissing=$allowXmissingRSplusSSdata -v mafCut=$mafCutoff -v allGenotypes="$genotypeColNames" -v rsGenotypes="$RSgenotypeColNames" -v ssGenotypes="$SSgenotypeColNames" -v nFieldsPerSample="$numFieldsPerSample" -v outFileBase="$path/$outputCsvFileNameBase" '

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
function numericDictStats(dict, statsResults, 
                                 #### local-scope variables below here
                                 key, max, secondMax, maxNum, secondMaxNum, total, numElems) {
  max=naVal
  secondMax=naVal
  
  maxNum=0
  secondMaxNum=0
  
  total=0
  numElems=0
  
  for (key in dict) { 
    numElems++
    total += dict[key]
    if (dict[key] > maxNum) {
      secondMax = max
      secondMaxNum = maxNum
      max = key
      maxNum = dict[key]
    }
    else if (dict[key] > secondMaxNum) {
      secondMax = key
      secondMaxNum = dict[key]
    }
  } # end for-loop
  
  statsResults["max"] = max
  statsResults["second_max"] = secondMax
  statsResults["max_num"] = maxNum
  statsResults["second_max_num"] = secondMaxNum
  statsResults["num_elems"] = numElems
  statsResults["total_all"] = total
  
}

##################
function calculateMAF(totalA, totalC, totalG, totalT, mafStatsArr, totalIndel1, totalIndel2,
                      ####  local-scope variables below here
                      max, secondMax) {

  # Calculate MAF pseudocode: 
  # 
  #   select max, secondMax from A, C, T, G, indel1, indel2 (passed in as arguments)
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


  if (totalIndel1) {
    if (totalIndel1 >= max) {
      secondMax=max
      max=totalIndel1
    }
    else if (totalIndel1 > secondMax) {
      secondMax=totalIndel1
    }
  }

  if (totalIndel2) {
    if (totalIndel2 >= max) {
      secondMax=max
      max=totalIndel2
    }
    else if (totalIndel2 > secondMax) {
      secondMax=totalIndel2
    }
  }


  mafStatsArr["max"] = max
  mafStatsArr["second_max"] = secondMax
  mafStatsArr["maf"] = secondMax / (secondMax + max)
  
  return mafStatsArr["maf"]
}

##################
function generateGenoStatsClassic(allGenoIndexes, naValue, statsArr,
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
function generateGenoStats(allGenoIndexes, naValue, isIndel, statsArr,
                           ####  local-scope variables below here
                           totalA, totalC, totalG, totalT, totalEmpty, i, maf, indelsDict, indelsStats, tmpArr, key) {
    
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
        # split saves a dictionary of values and numerical keys into the second input parameter
        split($i, tmpArr, "/")
        
        for (key in tmpArr) {
          if (length(tmpArr[key]) == 1) {
            if (tmpArr[key] == "A") { totalA++ }
            else if (tmpArr[key] == "C") { totalC++ }
            else if (tmpArr[key] == "G") { totalG++ }
            else if (tmpArr[key] == "T") { totalT++ }
          }
          else {
            # building a dictionary of indels used in this row, where the indels are the keys and the values are their total number of occurences in the row
            if (tmpArr[key] in indelsDict) {
              indelsDict[tmpArr[key]] += 1
            }
            else {
              indelsDict[tmpArr[key]] = 1
            }
          }
        } # end inner for-loop
        
        # for this, the indelsStats dictionary is built up with each successive call
        numericDictStats(indelsDict, indelsStats)
                
      }
    } # end for-loop

    if (isIndel) {
      maf = calculateMAF(totalA, totalC, totalG, totalT, mafStatsArr, indelsStats["max_num"], indelsStats["second_max_num"])
    }
    else {
      maf = calculateMAF(totalA, totalC, totalG, totalT, mafStatsArr)      
    }
    
    statsArr["A"] = totalA
    statsArr["C"] = totalC
    statsArr["G"] = totalG
    statsArr["T"] = totalT
    statsArr["empty"] = totalEmpty
    statsArr["sum"] = totalA + totalC + totalG + totalT
    statsArr["max"] = mafStatsArr["max"]
    statsArr["second_max"] = mafStatsArr["second_max"]
    statsArr["MAF"] = maf
    
    if (isIndel) {
      statsArr["indel1"] = indelsStats["max_num"]
      statsArr["indel2"] = indelsStats["second_max_num"]
      statsArr["other_indels"] = indelsStats["total_all"] - indelsStats["max_num"] - indelsStats["second_max_num"]
      statsArr["indel1_val"] = indelsStats["max"]
      statsArr["indel2_val"] = indelsStats["second_max"]      
      statsArr["sum"] += indelsStats["total_all"]
    }
    else {
      statsArr["indel1"] = 0
      statsArr["indel2"] = 0
      statsArr["other_indels"] = 0
      statsArr["indel1_val"] = naVal
      statsArr["indel2_val"] = naVal
    }
    
    return maf
}

##################
BEGIN { 

  OFS = FS # set the ouput field separator to equal the input field separator

  outFileStats = outFileBase"-stats.csv"
  outFileFull = outFileBase"-full.csv"
  outFileMin1 = outFileBase"-minimal1.csv"
  outFileMin2 = outFileBase"-minimal2.csv"
  outFileGt = outFileBase"-genotype.csv"
  outFileRatio = outFileBase"-ratio.csv"
  
  if (indelStatus == 2 || indelStatus == 3) {
    indelOutFileStats = outFileBase"-indels-stats.csv"
    indelOutFileFull = outFileBase"-indels-full.csv"
    indelOutFileMin1 = outFileBase"-indels-geno+ro+ao.csv"
    indelOutFileMin2 = outFileBase"-indels-geno+ratio.csv"
    indelOutFileGt = outFileBase"-indels-genotype.csv"
    indelOutFileRatio = outFileBase"-indels-ratio.csv"  
    
    indelsReplacedOutFileFull = outFileBase"-indels-replaced-full.csv"
    indelsReplacedOutFileMin1 = outFileBase"-indels-replaced-geno+ro+ao.csv"
    indelsReplacedOutFileMin2 = outFileBase"-indels-replaced-geno+ratio.csv"
    indelsReplacedOutFileGt = outFileBase"-indels-replaced-genotype.csv"
  }
  
  # split saves a dictionary of values and numerical keys into the second input parameter
  numAllGenos=split(allGenotypes,allGenos,",")
  numRsGenos=split(rsGenotypes,rsGenos,",")
  numSsGenos=split(ssGenotypes,ssGenos,",")
  numAllCols=split(header,allCols,",")

  if (indelStatus == 0) {
    indelMessage = "print both indels and non-indels to the same file."
  }
  else if (indelStatus == 1) {
    indelMessage = "keep only non-indels."
  }
  else if (indelStatus == 2) {
    indelMessage = "keep only indels."
  }
  else if (indelStatus == 3) {
    indelMessage = "print indel and non-indel output to separate files."
  }
  
  print "-------------------------"
  print "Log messages printed to stdout. Data output printed to the provided output file."
  print " "  
  print "Saving indels or non-indels?: "indelMessage
  print " "  
  print "Number of genotype columns: "numAllGenos
  print "Number of RS genotype columns: "numRsGenos
  print "Number of SS genotype columns: "numSsGenos
  print "Total number of columns: "numAllCols
  print " "
  refIndex=indexOf(refCol,allCols)
  altIndex=indexOf(altCol,allCols)
  
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
  mainPrintedLines=0
  indelPrintedLines=0
  nonStandardRows=0
  
  expectedNumFields=0
}

################## Main
{ 
  lineNum++;

  isIndel = (length($refIndex) > 1) 
    
  # Pseudocode, for each line (might be out of date since adding indel capabilities):
  #     if header line, print
  #     if REF length > 1, remove row
  #     for each sample
  #       if RS and missing, ++rsMissingTally
  #       if SS and missing, ++ssMissingTally
  #     if rsMissingTally + ssMissingTally > 4, remove row
  #     calculate MAF
  #     if maf < mafCutoff, remove row
  
  if (NF <= 1) {
    # skip row if it has one or zero fields
    nonStandardRows++
    next
  }
  # assume this is the header line because expectedNumFields has not been set yet
  else if (expectedNumFields == 0) {

    expectedNumFields = NF
    
    # Get the indices of the columns to print for the additional csv output files.
    # This is an example of the input column naming convention for just one sample (beginning in column 5):
    # CHROM,POS,REF,ALT,sample1_GT,sample1_DP,sample1_RO,sample1_AO,sample1_ratio,
     
    outputMin1Indexes="1,2,3,4" # this assumes four meta-columns
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
  
      
    if (indelStatus == 0 || indelStatus == 1 || indelStatus == 3) {
      if (indelStatus == 0) {
        # include indels in the main output stats file
        print $1","$2","$3","$4",A,C,G,T,indel1,indel2,other_indels,indel1_val,indel2_val,empty,sum,max,second_max,MAF" > outFileStats
      }
      else {
        print $1","$2","$3","$4",A,C,G,T,empty,sum,max,second_max,MAF" > outFileStats      
      }

      print $0 > outFileFull

      print min1Header > outFileMin1
      print min2Header > outFileMin2
      print gtHeader > outFileGt
      print ratioHeader > outFileRatio
    }
    
    if (indelStatus == 2 || indelStatus == 3) {
      print $1","$2","$3","$4",A,C,G,T,indel1,indel2,other_indels,indel1_val,indel2_val,empty,sum,max,second_max,MAF" > indelOutFileStats
      print $0 > indelOutFileFull

      print $0 > indelsReplacedOutFileFull

      print min1Header > indelOutFileMin1
      print min2Header > indelOutFileMin2
      print gtHeader > indelOutFileGt
      print ratioHeader > indelOutFileRatio

      print min1Header > indelsReplacedOutFileMin1
      print min2Header > indelsReplacedOutFileMin2
      print gtHeader > indelsReplacedOutFileGt

    }
        
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
  } # end of dealing with the header line
  # skip row if it does not have the expected number of fields
  else if (NF != expectedNumFields) {
    #print lineNum" : NF != expectedNumField, SKIP"
    nonStandardRows++
    next
  }
  # exclude all lines with indels
  else if (isIndel && indelStatus == 1) { 
    numIndels++   
    #print lineNum" : indel, SKIP"
    next  # skip row if it is an indel
  }
  # exclude all lines that are not indels
  else if (!isIndel && indelStatus == 2) {
    #print lineNum" : non-indel, SKIP"
    next  # skip row if it is a non-indel
  }
  else {
    
    if (isIndel) { 
      numIndels++
      
      # set all indels that do not match either the ref or the alt to missing data
      for (i in allGenoIndexes) {
        if ($i != naVal) {
          # split saves a dictionary of values and numerical keys into the second input parameter
          split($i, tmpArr, "/")
          
          isOtherIndel=0
          
          for (key in tmpArr) {
            if (tmpArr[key] != $refIndex && tmpArr[key] != $altIndex) { isOtherIndel=1; break }
          }
          
          if (isOtherIndel) {
            for (k=0; k < nFieldsPerSample; k++) {
              $(i+k) = naVal
            }
          }
          
        }          
      } # end outer for-loop

    }
    
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
    
    maf = generateGenoStats(allGenoIndexes, naVal, isIndel, statsArr) 

    if (maf < mafCut) {
      skippedMaf++
      #print lineNum" : maf="maf", below cutoff, SKIP"
      next  # skip row if maf is below the maf cutoff
    }
    
    # the row passes all filtering, so print it to the output files in various formats
     
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


    if (indelStatus == 0 || (!isIndel && (indelStatus == 1 || indelStatus == 3)) ) {
      if (indelStatus == 0) {
        # include indels in the main output stats file
        # print the stats file. This looks weird but the quotes are around the commas
        print $1","$2","$3","$4","statsArr["A"]","statsArr["C"]","statsArr["G"]","statsArr["T"]","statsArr["indel1"]","statsArr["indel2"]","statsArr["other_indels"]","statsArr["indel1_val"]","statsArr["indel2_val"]","statsArr["empty"]","statsArr["sum"]","statsArr["max"]","statsArr["second_max"]","statsArr["MAF"] >> outFileStats
      }
      else {
        # do not include indels in the main output stats file
        print $1","$2","$3","$4","statsArr["A"]","statsArr["C"]","statsArr["G"]","statsArr["T"]","statsArr["empty"]","statsArr["sum"]","statsArr["max"]","statsArr["second_max"]","statsArr["MAF"] >> outFileStats
      }

      print $0 >> outFileFull
      print min1Line >> outFileMin1
      print min2Line >> outFileMin2
      print gtLine >> outFileGt
      print ratioLine >> outFileRatio
      
      mainPrintedLines++
    }
       
    if (isIndel && (indelStatus == 2 || indelStatus == 3)) {
      # print the stats file. This looks weird but the quotes are around the commas
      print $1","$2","$3","$4","statsArr["A"]","statsArr["C"]","statsArr["G"]","statsArr["T"]","statsArr["indel1"]","statsArr["indel2"]","statsArr["other_indels"]","statsArr["indel1_val"]","statsArr["indel2_val"]","statsArr["empty"]","statsArr["sum"]","statsArr["max"]","statsArr["second_max"]","statsArr["MAF"] >> indelOutFileStats
      
      print $0 >> indelOutFileFull
      print min1Line >> indelOutFileMin1
      print min2Line >> indelOutFileMin2
      print gtLine >> indelOutFileGt
      print ratioLine >> indelOutFileRatio

      # now that all files which receive the non-replaced indels have been written to, replace the indels and write to those files: 
      
      for (i in allGenoIndexes) {
        if ($i != naVal) {
          nKeys = split($i, tmpArr, "/")
          
          for (key=1; key <= nKeys; key++) {
            if (tmpArr[key] == $refIndex) {
              tmpArr[key] = "R"
            }
            else if (tmpArr[key] == $altIndex) {
              tmpArr[key] = "A"
            }
            
            if (nKeys >= 2) {
              $i = tmpArr[1]"/"tmpArr[2]
            }
          } # end inner for-loop
          
        }
      } # end outer for-loop
      
      min1Line=""
      min2Line=""
      gtLine=""
      
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
      } # end for-loop

      print $0 >> indelsReplacedOutFileFull
      print min1Line >> indelsReplacedOutFileMin1
      print min2Line >> indelsReplacedOutFileMin2
      print gtLine >> indelsReplacedOutFileGt
      
      indelPrintedLines++
    }

  }
    

}
##################
END {

  print "-------------------------"
  print "Finished filtering."
  print "final lineNum="lineNum
  print " "
  print "("nonStandardRows" lines skipped because they were non-standard)."
  print lineNum-1-nonStandardRows" total input lines regarded as data lines."
  print numIndels" lines with indels were found."
  print mainPrintedLines" data lines printed to the main output files."
  print indelPrintedLines" data lines printed to the indel output files."
  print " "

  print nonStandardRows" lines skipped because they were non-standard."
  
  totalSkipped=nonStandardRows
  
  if (indelStatus == 1) {
    print numIndels" lines skipped because of indels."
    totalSkipped += numIndels
  }
  else if (indelStatus == 2) {
    numNonIndels=lineNum-1-nonStandardRows-numIndels
    print numNonIndels" lines skipped because they were not indels."
    totalSkipped += numNonIndels
  }
  
  print skippedMissingData" lines skipped because of too much missing data."
  print skippedMaf" lines skipped because their MAF was below the MAF cutoff."
  
  totalSkipped += skippedMissingData+skippedMaf
  
  print totalSkipped" lines skipped in total."

  
}'





