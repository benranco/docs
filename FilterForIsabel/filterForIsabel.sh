#!/bin/bash

basePath="/run/media/benrancourt/27a1f1ff-6590-4f38-aa2c-13f6f1e0e771/ben/SNPpipeline-isabel"

reportsSubDir="reports/reports-backup/done-exceptMissingData/filter_redo"

# part 1 input:
filledReportName="filled_report.csv"
#filledReportName="report_filtered_by_depth_recoded.csv"
# part 1 output:
part1FilteredReportName="report_filtered_pt1.csv"

# these files are output by getDepthStats-parallel.R (part 2), and read by part 3:
depthReportName=$part1FilteredReportName"_depth.csv"
depthDetailedReportName=$part1FilteredReportName"_depth_detailed.csv"

# part 3 output:
#filteredByDepthReportName="report_filtered_by_depth.csv"
#filteredByDepthRecodedReportName="report_filtered_by_depth_recoded.csv"
filteredByDepthReportName="report_filtered_more.csv"
filteredByDepthRecodedReportName="report_filtered_more_recoded.csv"
#filteredByDepthRecodedReportName=test2000lines.csv
# part 4 output:
finalFastaFileName="report_filtered_by_depth_recoded.fasta"
#finalFastaFileName=test2000lines.fasta

#Rscript ./filterForIsabel_part1.R $basePath $reportsSubDir $filledReportName $part1FilteredReportName

#Rscript ./getDepthStats-parallel.R $basePath $reportsSubDir $part1FilteredReportName

#Rscript ./filterForIsabel_part3.R $basePath $reportsSubDir $depthReportName $depthDetailedReportName $part1FilteredReportName $filteredByDepthReportName $filteredByDepthRecodedReportName 

Rscript ./filterForIsabel_part4.R $basePath $reportsSubDir $filteredByDepthRecodedReportName $finalFastaFileName

