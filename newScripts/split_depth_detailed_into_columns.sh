#!/bin/bash
echo "running split_depth_detailed_into_columns.sh"

awk 'NR == 1' MAF_cutoff_report_depth_stats_detailed.csv > tmp_header.csv
awk 'NR > 1' MAF_cutoff_report_depth_stats_detailed.csv > tmp_content.csv

# get rid of all quotes
sed -i 's/\"//g' tmp_header.csv 
sed -i 's/\"//g' tmp_content.csv 

# temporarily remove CHROM,POS,REF from the header line, we will add them back when finished
sed -i 's/CHROM,POS,REF//g' tmp_header.csv 

# It's easier to read the below regex in the sed command without all the backslashes before the ( and ) characters, so here it is without those backslashes:
#sed 's/((,[[:alnum:]_.-]*)(-DP);(DPB);(AO);(RO))/\2\3\2-\4\2-\5\2-\6/g' tmp_header.csv > MAF_cutoff_report_depth_stats_detailed_split.csv
# Here's an example of what it would do. Change this line: 
# ,C1_filtered-DP;DPB;AO;RO,C10_filtered-DP;DPB;AO;RO,C11_filtered-DP;DPB;AO;RO
# to this line: 
# ,C1-DP,C1-DPB,C1-AO,C1-RO,C10-DP,C10-DPB,C10-AO,C10-RO,C11-DP,C11-DPB,C11-AO,C11-RO
# It assumes that sample names are made of of alphanumeric or _ or . or - characters. 

sed 's/\(\(,[[:alnum:]_.-]*\)\(-DP\);\(DPB\);\(AO\);\(RO\)\)/\2\3\2-\4\2-\5\2-\6/g' tmp_header.csv > MAF_cutoff_report_depth_stats_detailed_split.csv

# add the CHROM,POS,REF back to the beginning of the header line
sed -i 's/.*$/CHROM,POS,REF&/' MAF_cutoff_report_depth_stats_detailed_split.csv


sed 's/;/,/g' tmp_content.csv >> MAF_cutoff_report_depth_stats_detailed_split.csv

rm tmp_header.csv
rm tmp_content.csv

