
options(stringsAsFactors = FALSE, warn = 1)

write("Running removeRows.R.", stdout())


# ########################################################
# Input Parameters:

path <- "/home/benrancourt/Downloads/r38-remove86genes/5272-remove86genes-try2"

rowsToRemoveFileName <- "genesToRemove.csv"
rowsRemovedOutputFileName <- "rowsThatWereRemoved.csv"

#dataToModifyFileName <- "R-38-M32178-38329_p0.01-uniqueGenes-postLepMAP2-stage2-allChr-LG5andLG11stage3addedManually-afterLG13combinedToLG12.csv"
#mainOutputFileName <- "R-38-M32178-38329_p0.01-uniqueGenes-postLepMAP2-stage2-allChr-LG5andLG11stage3addedManually-afterLG13combinedToLG12-86genesRemoved.csv"

#dataToModifyFileName <- "R-38-M32178-38329_p0.01-uniqueGenes-postLepMAP2-stage2-allChr-LG5andLG11stage3addedManually-beforeLG13combinedToLG12.csv"
#mainOutputFileName <- "R-38-M32178-38329_p0.01-uniqueGenes-postLepMAP2-stage2-allChr-LG5andLG11stage3addedManually-beforeLG13combinedToLG12-86genesRemoved.csv"

#dataToModifyFileName <- "Data-5,092originalLGs+751-stage1output-postLepMAP2-allChr-stage2output.csv"
#mainOutputFileName <- "Data-5,092originalLGs+751-stage1output-postLepMAP2-allChr-stage2output-86genesRemoved.csv"

dataToModifyFileName <- "Data-5,092originalLGs+751-finalOutputOfLepMAP2-manuallyCombined.csv"
mainOutputFileName <- "Data-5,092originalLGs+751-finalOutputOfLepMAP2-manuallyCombined-86genesRemoved.csv"

# ########################################################


rowsToRemove <- read.csv(paste(path.expand(path),rowsToRemoveFileName, sep="/"),header=TRUE,col.names=c("gene","position","linkage_group"))
dataToModify <- read.csv(paste(path.expand(path),dataToModifyFileName, sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

rowsThatWereRemoved <- dataToModify[(dataToModify$gene %in% rowsToRemove$gene), ]
dataToModify <- dataToModify[!(dataToModify$gene %in% rowsToRemove$gene), ]

write.csv(rowsThatWereRemoved, paste(path.expand(path),rowsRemovedOutputFileName, sep="/"), row.names=FALSE)
write.csv(dataToModify, paste(path.expand(path),mainOutputFileName, sep="/"), row.names=FALSE)

write(paste0("FINISHED."), stdout())


