This is Ben's little readme for BLAST.
--------------------------------------

INSTALLING BLAST

This version of BLAST needs to be installed on CentOS Linux.
If your virtual machine is based on an image that doesn't already have BLAST installed, there is an .tar.gz archive 
in the blast-2.11 folder. Simply unzip it (tar xvzf filename.tar.gz) in the folder you want to install it, and then 
add it's bin folder to your PATH environment variable. For example, if you unzipped the tar.gz file inside this 
folder in your home directory 
~/software/blast-2.11 
then just edit your ~/.bash_profile file and add the following to the end of your PATH variable:
:$HOME/software/blast-2.11/ncbi-blast-2.11.0+/bin

Then exit your virtual machine and restart it for thist to take effect.

--------------------------------------

USING THE NCBI NR DATABASE

The NCBI NR dabase has been downloaded to the ncbi_nr_db folder. 
To use it with your blast run, the most straightforward way to do it is to copy the .ncbirc file in this directory 
into the directory from which you'll be running blast, and edit it. If you don't see the .ncbirc file, use this 
command from the current directory to see it: 
ls -la

If you need to download other NCBI databases, use the perl script provided in the bin directory of your blast 
installation to do it. Run it from the folder in which you want to download the database. 
Here is the command that I used to download the NR database. The nohup command at the beginning is to ensure the 
process continues with "no hangups" even if the ssh session is interrupted: 
nohup ../blast-2.11/ncbi-blast-2.11.0+/bin/update_blastdb.pl --verbose --timeout 180 nr > update_blastdb.log &

To keep track of its progress, you can tail the update_blastdb.log:
tail -f update_blastdb.log

Or you can run the top command to see if the perl process is still running.

To get  a list of usage instructions, execute the update_blastdb.pl script without any input parameters.

--- 
To run blastp...
First, copy and edit the .ncbirc file as explained above.
An example basic blastp command, in multi-threaded mode using 16 threads, and running it through nohup to ensure
"no hangups" even if the ssh session is terminated:
nohup blastp -num_threads 16 -db nr -query inputdata.fasta -out results.out > console_out.log &

-----------------------------------------------

USER MANUEL, OTHER LINKS

Below are some links which might be useful.

---
Downloading and installing BLAST:
https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
Basic instructions, see CentOS (seems you just need to download the tar.gz, decompress, and add it's bin folder to your PATH):
http://www.metagenomics.wiki/tools/blast/install

---
BLAST+ user manual: 
https://www.ncbi.nlm.nih.gov/books/NBK279691/
https://www.ncbi.nlm.nih.gov/books/NBK279690/
Configuring: https://www.ncbi.nlm.nih.gov/books/NBK279695/
Overview of features: https://www.ncbi.nlm.nih.gov/books/NBK279668/
Quick start: https://www.ncbi.nlm.nih.gov/books/NBK279680/
The online blastp tool: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=MultiSensor

As for the RAM and processing power for BLAST, you may check the link as below. Thanks. Jun-Jun
https://doctorlib.info/medical/blast/13.html

Why is blast taking so long?:
https://secure.clcbio.com/helpspot/index.php?pg=kb.page&id=248

An example basic blastp command, in multi-threaded mode using 16 threads, and running it through nohup to ensure 
"no hangups" even if the ssh session is terminated: 
nohup blastp -num_threads 16 -db nr -query inputdata.fasta -out results.out > console_out.log &

