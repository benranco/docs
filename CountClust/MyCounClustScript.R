

#########################################################
# installation. I entered the following commands to install CountClust, as guided by https://bioconductor.org/packages/release/bioc/vignettes/CountClust/inst/doc/count-clust.html:

#source("https://bioconductor.org/biocLite.R")
#biocLite("CountClust")
#biocLite("BiocUpgrade")
#biocLite("CountClust")
#install.packages("devtools")
#library(devtools)
#install_github("TaddyLab/maptpx")
#install_github('kkdey/CountClust')
#library(CountClust)

#########################################################
#options(stringsAsFactors = FALSE, warn = 1)
options(scipen = 999) # disable scientific notation for large numbers

#write("Running filterForHolly.R.", stdout())

#########################################################
# Initial Input Parameters:

#args <- commandArgs(trailingOnly = TRUE)

path <- "/home/benrancourt/Desktop/holly"

inputFileName <- "RPKM_matrix_8101_8112_8206_8214.csv"
#inputFileName <- "rpkm_small.csv"

outputFileNameBase <- "RPKM"

#rdaTempFileName <- "rpkm_small.rda"

input <- read.csv(paste(path.expand(path),inputFileName,sep="/"),header=TRUE) ##, check.names=FALSE) # using check.names=FALSE in case the column names have dashes (-) in them. This will prevent them from being converted to periods. However, a column name with a dash in it will not be able to be used as a variable name, so we'll have to refer to columns by their index if accessing them.)

geneNames <- input[,1]
input <- input[,-1]

library(CountClust)

library(devtools)

#gomOutput <- FitGoM(t(input),
#            K=4, tol=0.1,
#            path_rda=paste(path.expand(path),rdaTempFileName,sep="/"))

#gomOutput <- FitGoM(t(input),K=c(4,8,12), tol=0.1)

gomOutput <- FitGoM(t(input),K=c(15:35), tol=0.1)

for(i in 15:35) {
  clustName <- paste0("clust_",i)
  omg <- gomOutput[[clustName]]$omega
  tht <- gomOutput[[clustName]]$theta
  tht <- as.data.frame(tht)
  tht <- cbind(geneNames,tht)

  write.csv(omg, paste(path.expand(path), paste0(outputFileNameBase,"-omega",i,"-withRowNames.csv"),sep="/"), row.names=TRUE)
  write.csv(tht, paste(path.expand(path), paste0(outputFileNameBase,"-theta",i,"-withRowNames.csv"),sep="/"), row.names=TRUE)

}




