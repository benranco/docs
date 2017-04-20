This README contains:
- information on how to run the script used to generate Circos graphs
- information on how to install Circos
- helpful links for learning more about Circos

Note: much of this documentation is replicated from the documentation inside the run
script.

##########################################################################################
# Assuming a complete and working installation of Circos on your Linux computer, these are 
# the files and directory structure that should be present to run this script. In a folder, 
# called whatever you'd like it to be called (eg. CircosScript):
# Files:
# - README-CIRCOS-SCRIPT.TXT (not necessary, but good to keep here for reference)
# - run  (this is the main script that starts the whole process)
# - convertCsvToCircosInput.R (called by the run script, converts your .csv into Circos format)
# Sub-Folder:
# - config (this folder contains custom configuration files that will be used by Circos
#
# In the config sub-folder:
# Files:
# - circos-noTicks.conf
# - circos-withTicks.conf
# - ideogram-noTicks.conf
# - ideogram-withTicks.conf
# - ticks-noTicks.conf
# - ticks-withTicks.conf
#
# These .conf configuration files are set up to render the circos graphs in a certain way. 
# You shouldn't need to edit these files, but if you want to customize some element of the 
# graphs, read through these files and also browse some of the circos website links below.
#
# You will also supply a .csv file with the data you'd like to display in the graph. More 
# details on the format of that file below.
#
###########################################################################################
# The main script is called "run", and it should be found in the same folder as this 
# readme file.
#
# The script creates Circos graphs using Circos with some customized configuration
# files. It runs Circos twice, first to create a version without position ticks
# along the edges of the circle, and second to create a version with postion
# ticks. It then uses some command-line Inkscape and ImageMagick commands to 
# copy the Circos output images to 600dpi tiff images of different sizes.
#
# The custom Circos configuration files this script expects to see are:
#   circos-noTicks.conf, ideogram-noTicks.conf, ticks-noTicks.conf
#   circos-withTicks.conf, ideogram-withTicks.conf, ticks-withTicks.conf
#
# It also expects to see a .csv file (specified in the csvFilePath and csvFileName
# input parameters below), this is the format in which you will provide the input data
# that you want to graph. The csv must contains four columns, organized:
#   Group1 LG, Group1 position, Group2 LG, Group2 position.
# The two LG columns must be numeric (integer), indicating linkage group number (eg. 3, 
# not "LG3"). Do not use any non-numeric characters in the linkage group ids. The two 
# position columns must also be numeric.
# The four data columns will be expected to have headers, but their names are unimportant.
#
# Before running Circos, this script calls the R script "convertCsvToCircosInput.R", 
# which converts the input .csv file into two data files which will be read by Circos.
# These two data files are named:
#   mykaryotype.txt, mydata.txt
# They will be deleted after running Circos to generate the graphs. 
#
# To run this script:
# First, review and set the input parameters (open and edit the "run" script in a text 
# editor) appropriately.
# Then, from a command-line terminal, within the directory that the script is contained in, 
# type: ./run
# The custom Circos configuration files named above must be in a directory named "config"
# within the directory that contains the run script. 
#
# If you encounter errors while running, open the output text file run.out
# in a text editor to see if it logs any helpful details on the errors.
#
##########################################################################################

##########################################################################################
# NOTES ON INSTALLING CIRCOS (If it isn't already installed):
#
# This can be a bit of a pain.
#
# To install Circos and the Circos-Tools package (on CentOS), download the current version 
# of each from the Circos download page and unzip them. Move the "tools" folder and all 
# its contents from the unzipped circos-tools location into the "circos" folder of the 
# unzipped circos location.
#
# To install the perl modules on CentOS for Circos, from the command-line, follow the 
# following instructions:
#
# Type:
#    perl -v # to determine if your perl version is at least 5.8 or newer
# Navigate to whichever directory you extracted circos to, and navigate to its bin directory  
# (eg. /home/benrancourt/Desktop/circos/circos-0.69-3/bin), then type:
#    ./circos modules
# Copy the list of modules to a text file, so you can refer back to it to see which modules are missing.
# Execute the following command-line commands:
#    sudo yum install perl-CPAN
#    sudo yum install gd-devel # this includes /usr/bin/gdlib-config, which is needed by some CPAN modules
#    sudo cpan Module::Build
#    sudo perl -MCPAN -e shell # opens the CPAN shell
# To install the missing modules from the list, for each module, from within the CPAN shell, type:
#    install modulename
# (eg: install Math::Bezier)
# NOTE: You might want to install module Statistics::Descriptive, which is needed by circos' tableviewer tool 
# (distributed in the tools package) but which might not be included in the list of modules needed by circos
# (this script doesn't use the tableviewer tool, by the way):
#    install Statistics::Descriptive 
# Type q to quit.
# Make sure you now have all the necessary modules by re-entering this command in the command-line:
#    ./circos modules
# If some modules are still missing, go back into the CPAN shell and try installing them again, then read the 
# ouput of the install operation and look for errors, such as modules that are required by this module which 
# for some reason weren't automatically installed, or required non-perl linux libraries that need to be 
# installed via yum. Then try installing those required modules, or resolving those dependency issues and retry 
# installing the module that didn't work.
##########################################################################################

##########################################################################################
These are some links that I found helpful in learning about Circos, both for how to install 
it and how to use it.
 
Format of data input files for Circos:
http://circos.ca/documentation/tutorials/configuration/data_files/
Quickstart guide for Circos, very useful:
http://circos.ca/documentation/tutorials/quick_start/

Circos command-line parameters:
http://circos.ca/tutorials/lessons/configuration/runtime_parameters/
Circos configuration files:
http://circos.ca/tutorials/lessons/configuration/configuration_files/
Circos data file format:
http://circos.ca/documentation/tutorials/configuration/data_files/
Circos best practices:
http://circos.ca/documentation/tutorials/reference/best_practices/
Circos quick start:
http://circos.ca/documentation/tutorials/quick_start/

Circos image examples:
http://circos.ca/images/
Circos guide to use in genomics and sequence analysis:
http://circos.ca/guide/genomic/

Circos installing perl modules:
http://circos.ca/documentation/tutorials/configuration/perl_and_modules/

##########################################################################################

