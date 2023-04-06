
grep -v ",NA" MAF_cutoff_report.csv > MAF_cutoff_report-withoutNAs.csv

head -n 1 MAF_cutoff_report.csv > MAF_cutoff_report-withoutIndels.csv
grep -v -E '".*",[0-9]*,"[ACTG][ACTG][ACTG]*",' MAF_cutoff_report.csv >> MAF_cutoff_report-withoutIndels.csv

