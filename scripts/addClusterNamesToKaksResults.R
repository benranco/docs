options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)

write("Running addClusternamesToKaksResults.R.", stdout())


#if(!require(stringr)) {
#  install.packages('stringr', repos='http://cran.us.r-project.org')
#}


#########################################
# This script assumes a few modifications to the a few of the columns in the input files 
# which are more convenient to do using search & replace in a text editor (to split the 
# values in a column by "-" or ".", for example).

# ########################################################
# Input Parameters:

path <- "/home/ben/pfc/junjun/ass4/step2/kaks"


inputFileKaks <- "NB-ARC158-lngth150-LPss-m-KaKs-seqIDsSplit.csv"
inputFileClusters <- "Envelop158-clusters-seqIds-toMerge.csv"

outputFile <- "NB-ARC158-lngth150-LPss-m-KaKs-withClusterIDs.csv"



# ########################################################
# Execution Code:

dataKaks <- read.csv(paste(path.expand(path),inputFileKaks, sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

dataCluster <- read.csv(paste(path.expand(path),inputFileClusters, sep="/"),header=TRUE, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.

dataCluster <- dataCluster[,1:2]


#seqPairs <- strsplit(dataKaks[,"Sequence"], split="-", fixed=TRUE)[[1]]
#dataKaks <- cbind(seqPairs[

withClusters <- merge(dataKaks, dataCluster, by.x="Sequence1", by.y="Name-CDS", all.x=TRUE, sort=FALSE)

# put the new last column (Cluster-ID) at the front:
ClusterID <- withClusters[,ncol(withClusters)]
withClusters <- cbind( ClusterID, withClusters[,1:ncol(withClusters)-1] )

write.csv(withClusters, file=paste(path.expand(path),outputFile,sep="/"), row.names=FALSE )

write(paste0("FINISHED."), stdout())


