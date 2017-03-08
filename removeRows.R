
options(stringsAsFactors = FALSE, warn = 1)

write("Running removeRows.R.", stdout())


# ########################################################
# Input Parameters:

path <- "/home/benrancourt/Downloads/r38-remove86genes/5272-remove86genes-try2/"

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


rowsToRemove <- read.csv(paste0(path,rowsToRemoveFileName),header=TRUE,col.names=c("gene","position","linkage_group"))
dataToModify <- read.csv(paste0(path,dataToModifyFileName),header=TRUE)

rowsThatWereRemoved <- dataToModify[(dataToModify$gene %in% rowsToRemove$gene), ]
dataToModify <- dataToModify[!(dataToModify$gene %in% rowsToRemove$gene), ]

write.csv(rowsThatWereRemoved, paste0(path,rowsRemovedOutputFileName), row.names=FALSE)
write.csv(dataToModify, paste0(path,mainOutputFileName), row.names=FALSE)

write(paste0("FINISHED."), stdout())


