options(stringsAsFactors = FALSE, warn = 1)

path <- "/home/benrancourt/Desktop/holly/FilesForBenMarch14-2019"

table1 <- read.csv(paste(path.expand(path),"Up_2x+_2wD_8101_vs_susc_(2188).txt",sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

table2 <- read.csv(paste(path.expand(path),"54K_CA_DF_transcriptome_BLAST_annotated.csv",sep="/"),header=TRUE, check.names=FALSE)

table3 <- read.csv(paste(path.expand(path),"Fold-change&FDR_pval_tol&susc_all_54K.csv",sep="/"),header=TRUE, check.names=FALSE)

table4 <- read.csv(paste(path.expand(path),"RPKM_matrix_8101_8112_8206_8214.csv",sep="/"),header=TRUE, check.names=FALSE)

idColName <- colnames(table1)[1]
colnames(table2)[1] <- idColName
colnames(table3)[1] <- idColName
colnames(table4)[1] <- idColName

message("merging tables 1 and 2")
report <- merge(table1, table2, by = idColName, all.x = TRUE, all.y = FALSE)
message("merging with table 3")
report <- merge(report, table3, by = idColName, all.x = TRUE, all.y = FALSE)
message("merging with table 4")
report <- merge(report, table4, by = idColName, all.x = TRUE, all.y = FALSE)

write.csv(report, paste(path.expand(path), "merged_for_Holly.csv", sep = "/"), row.names=FALSE)

message("merge completed")
