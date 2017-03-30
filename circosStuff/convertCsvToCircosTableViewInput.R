# ########################################################
# To install Circos and the Circos-Tools package (on CentOS), download the current version of each from the Circos download page and unzip them. Move the "tools" folder and all its contents from the unzipped circos-tools location into the "circos" folder of the unzipped circos location.
#
# To install perl modules on CentOS for Circos, from command-line:
# perl -v # to determine if your perl version is at least 5.8 or newer
# navigate to whichever directory you extracted circos to, and navigate to its bin directory  (eg. /home/benrancourt/Desktop/circos/circos-0.69-3/bin), then:
# ./circos modules
# copy the list of modules to a text file
# sudo yum install perl-CPAN
# sudo yum install gd-devel # this includes /usr/bin/gdlib-config, which is needed by some CPAN modules
# sudo cpan Module::Build
# sudo perl -MCPAN -e shell # opens the CPAN shell
# to install the missing modules from the list, for each module, type:
# install modulename
# (eg: install Math::Bezier)
# install module Statistics::Descriptive, which is needed by circos' tableviewer tool (distributed in the tools package) but which might not be included in the list of modules needed by circos:
# install Statistics::Descriptive 
# type q to quit
# make sure you now have all the necessary modules:
# ./circos modules
# if some modules are still missing, go back into the CPAN shell and try installing them again, then read the ouput of the install operation and look for errors, such as modules that are required by this module which for some reason weren't automatically installed, or required non-perl linux libraries that need to be installed via yum. Then try installing those required modules, or resolving those dependency issues and retry installing the module that didn't work.
# ########################################################

options(stringsAsFactors = FALSE, warn = 1)

write("Running convertCsvToCircosInput.R.", stdout())


# ########################################################
# Input Parameters:

path <- "/home/benrancourt/Downloads/Circos/"

inputFileName <- "LPmerged-9511-PitaEcht4492-for-CIRCOS-Feb8-2017-experimental.csv"

finalOutputFileName <- "LPmerged-9511-PitaEcht4492-for-CIRCOS-Feb8-2017-experimental.circos.txt"

input <- read.csv(paste0(path,inputFileName),header=TRUE)

LGsInGroup1 <- sort(unique(input[ , 1]))
LGsInGroup2 <- sort(unique(input[ , 2]))

write("LG's in group 1: ", stdout())
write(LGsInGroup1, stdout())
write("LG's in group 2: ", stdout())
write(LGsInGroup2, stdout())

if (length(LGsInGroup1) != length(LGsInGroup2))
{
  write("Warning: there is a different number of LG's in each species. Continuing anyway.", stdout())
}

tallies <- matrix(nrow = length(LGsInGroup1), ncol = length(LGsInGroup2))

for (i in 1:length(LGsInGroup1))
{
  # get a list of all elements in group2 that are matched with the current LG in group1
  matches <- input[ input[,1]==LGsInGroup1[i] , 2]
  
  for (k in 1:length(LGsInGroup2))
  {
    # count the number of times an LG from group2 is matched with the current group1 LG.
    numMatches <- sum(matches == LGsInGroup2[k])
    if(numMatches > 0)
    {
      tallies[LGsInGroup1[i], LGsInGroup2[k]] <- numMatches
    }
  }
}

tallies <- data.frame(tallies)

# arranged side by side:
group1_order <- c(1:length(LGsInGroup1))*2-1
group2_order <- c(1:length(LGsInGroup2))*2

# arranged like a basketball:
#group1_order <- c(1:length(LGsInGroup1))
#group2_order <- c(length(LGsInGroup2):1)+12

# create column names
#colnames(tallies) <- paste0(colnames(input)[2],"-LG",LGsInGroup2)
#colnames(tallies) <- paste0("lg-",LGsInGroup2)

#group2Labels <- paste0(colnames(input)[2],"-LG",LGsInGroup2)
group2Labels <- paste0("lg-",LGsInGroup2)

#group1Labels <- paste0(colnames(input)[1],"-LG",LGsInGroup1)
group1Labels <- paste0("LG",LGsInGroup1)

# add two new rows to the top, one to tell Circos the order to use for them, 
# and the other to tell Circos what to use for labels
tallies <- rbind(group2_order, group2Labels, tallies)

# add two new cols to the front, one to tell Circos the order to use for them, 
# and the other to tell Circos what to use for labels
tallies <- cbind(c("-","-",group1_order), c("-","-",group1Labels), tallies)

# replace all occurences of NA with "-"
tallies[is.na(tallies)] <- "-"

write.table(tallies, file= paste0(path,finalOutputFileName), append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# add col at beginning called data, with lg names
# create proper row and col names

# determine num unique lg's in col1
# determine num unique lg's in col2
# make sure they are the same = x
# create x by x dataframe
# for each lg in col1
#   make sublist of just those rows in that lg
#   for each lg in col2 of that sublist
#     count num
#     add num to dataframe at col1Num,col2num
# output dataframe as space or tab delimited file, of format:
#
# data lg1 lg2 lg3...
#  LG1 300 -   2
#  LG2 -   270 -
#  LG3 1   1   320
#  ...


write(paste0("================================================"), stdout())

write(paste0("FINISHED."), stdout())


