
---
Megagenome Analyzer 6 links:
from: https://bio.tools/megan
"Metagenome Analysis Software - MEGAN (MEtaGenome ANalyzer) is a new computer program that allows laptop analysis of large metagenomic datasets. In a preprocessing step, the set of DNA reads (or contigs) is compared against databases of known sequences using BLAST or another comparison tool. MEGAN can then be used to compute and interactively explore the taxonomical content of the dataset, employing the NCBI taxonomy to summarize and order the results."

https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/megan6/
https://software-ab.informatik.uni-tuebingen.de/download/megan6/welcome.html
https://software-ab.informatik.uni-tuebingen.de/download/megan6/manual.pdf
- requires jdk 11 or later (packed with it for windows and and mac)

Useful info for installing java-11: 
https://www.linode.com/docs/guides/how-to-install-openjdk-on-centos-8/

---

Ben's Notes: 
- MEGAN6 can be run from an openstack vm based on my CentOS7_GUI_VNC-Java14 image.
- I've already installed MEGAN6 in the following folders:
    public/genomics/junjun/Data-Dec2020/ben2021/a6/blastrun/forMegan6  (this was the original installation)
    public/genomics/junjun/MEGAN6   (copied from the original installation location, except for the megan-map files, which I'm not sure I need)
- to use, just mount either of these folders to /work3 or something: 
  sudo mount -t nfs 10.20.0.6:/Public/genomics/junjun/Data-Dec2020/ben2021/a6/blastrun/forMegan6 /work3
- to run MEGAN6, from within your VNC session, open a terminal and cd to /work3/megan6/megan, then run ./MEGAN
- MEGAN expects the input to be the output file of a BLASTp/BLASTx/BLASTn run, as well as the fasta reads file used for the BLAST run.
- MEGAN6 expects the BLAST input formats to be either in standard format (outfmt 0) or BLAST XML, and can also parse tabular format, with some caveats (see section 34.3 of the MEGAN6 manual).
- Below is some correspondence with Jun-Jun from my first time running MEGAN6 (originally in file ass6.instructions.txt): 


-----------------------------------
[This first usage of MEGAN6 was on data from BLASTp]
-----------------------------------
Hi Jun-Jun, 

That first screenshot was taken with all the default settings, including Min Support Percent = 0.05. 
I've attached a new screenshot after recalculating the LCA algorithm with the new Min Support Percent = 0.01, as you suggested.

I've left the other parameters as default: 
Min Score: 50.0
Max Expected: 0.01 
Min Percent Identity: 0.0
Top Percent: 10.0
Min Support Percent: (custom) 0.01
Min Support: 0
Use Min-Complexity Filter: off
LCA Algorithm: naive
Percent to cover: 100%
Read assignment mode: readCount
Use 16S Percent Identity Filter: off

These LCA parameters are described in section 28 Parameters Dialog of the manual (page 46):
https://software-ab.informatik.uni-tuebingen.de/download/megan6/manual.pdf

Let me know if you'd like to make any other customizations. 
---
Hello, Ben:

Great. The output is expected. Does it cover 100% of sequences. Can you trace each seq/read to a taxon, or as no hits, and not assigned?

Hello, Ben:

I suspect,  those seqs/reads that can not reach 'Min Support Percent' have gone to the category of 'not assigned'. If so, we still have a chance to lose those biodiversity with frequency below  'Min Support Percent'.
After you get seq IDs under each taxon category, we will see if it is true.

[
No Hits: 57239
Not Assigned: 123886
]
---

Now: 
- get list of sequence ids at each of these ten categories:
No hits (Assigned=57239)
Not assigned (Assigned=123886)
[Riboviria]
Bacteria (Assigned=852, Summarized=4637)
Amoebozoa (Assigned=12, Summarized=75)
Metazoa (Assigned=51, Summarized=4684)
Fungi (Assigned=2621, Summarized=52697)
Sar/Harosa (Assigned=68, Summarized=608)
chlorophyta (Assigned=35, Summarized=739)
embryophyta (Assigned=11738, Summarized=127769)

- it's easy to export fasta files of the seqs by clicking one of the buttons [I think it's the "AC" button, "Extract reads for selected nodes"], then do either of these to extract lists of seq ids:

grep ">" reads-Sar.fasta | sed 's/\>([^ ]*\) .*/\1/' | sed 's/>//' > seqIds-Sar.txt # using above, print out a list of seq Ids from a fasta file, one per line, with the starting ">" still there. The second sed command replaces all ">" with empty string.
grep ">" reads-Sar.fasta | sed 's/>\([^ ]*\) .*/\1/' > seqIds-Sar3.txt # does the same as above, but with just the one sed command, by moving the \( to exclude the initial ">" from the captured portion.

---
Hi Jun-Jun,

I've attached a zip file containing the reads (.fasta) and sequence Ids (.txt) in each of the categories you selected. The only category that wasn't represented in this version of the MEGAN6 tree based on the complete BLASTp output from the 399753 seqs is Riboviria. I assume that one was grouped into another category.

Here are the categories you selected along with statistics from MEGAN6. I'm not sure what the summarized means. The fasta files that I exported for each category contain only the assigned reads. [see above]

---
Hello, Ben:

Based on your data with Min Support Percent: (custom) 0.01, we need seq IDs under 11 categories as below. Each category with the assigned or summarized IDs with whatever number is larger. For examples, we need 4,637 seq IDs for Bacteria. Thanks.   Jun-Jun
---
Redo, but get Summarized: 
- get list of sequence ids at each of these eleven categories:
No hits (Assigned=57239, Summarized=57239)
Not assigned (Assigned=123886, Summarized=123886)
Pinus nigra virus1 (Assigned=129, Summarized=129)
Bacteria (Assigned=852, Summarized=4637)
Amoebozoa (Assigned=12, Summarized=75)
Discoba (Assigned=2, Summarized=57)
Metazoa (Assigned=51, Summarized=4684)
Fungi (Assigned=2621, Summarized=52697)
Sar/Harosa (Assigned=68, Summarized=608)
chlorophyta (Assigned=35, Summarized=739)
embryophyta (Assigned=11738, Summarized=127769)

grep ">" reads-Sar-includeSummarized.fasta | sed 's/>\([^ ]*\) .*/\1/' > seqIds-Sar-includeSummarized.txt

---
[from Jun-Jun]
The total number is 399,753. Below total ID number of listed categories is only 372,334. About 27,419 IDs are not included.
---

=57239+123886+129+4637+75+57+4684+52697+608+739+127769 = 372520
=399753-372520 = 27233

Hi Jun-Jun, 

I've solved why there is a discrepancy between the number of reads in our eleven categories (total 372520) and the total number of reads in our input (399753).
The discrepancy is: 399753-372520 = 27233

The missing 27233 reads are accounted for by these two groups:

Assigned reads in these parent categories, which are further up from our end nodes: 
    Opisthokonta (Assigned=579, Summarized=57960)
    Viridiplantae (Assigned=279, Summarized=129669)
    Eukaryota (Assigned=11251, Summarized=199620)
    cellular organisms (Assigned=13244, Summarized=217501)
    NCBI (Assigned=956, Summarized=399753)

Assigned reads in these categories which are in-line with no branching children, and therefore very hard to see in our tree:
Two unseen categories with assigned reads in-line between Viridiplantae and Embryophta:
    Streptophyta (Assigned=526, Summarized=128651)
    Streptophytina (Assigned=356, Summarized=128125)
Three unseen categories with assigned reads in-line between NCBI and Pinus nigra virus1:
    Viruses (Assigned=18, Summarized=171)
    Riboviria (Assigned=22, Summarized=153)
    Ortervirales (Assigned=2, Summarized=131)

I've attached a screenshot revealing all these categories.

When we add up all these missing assigned reads, they add to the exact number of unaccounted-for reads: 
=579+279+11251+13244+956+526+356+18+22+2 = 27233

How would you like me to proceed with extracting seq Ids from these categories, along with the eleven categories you've already identified? Should I take the summarized reads from your eleven categories, and the assigned reads from these intermediate categories?

-----------
NOW, REDOING THINGS A LITTLE, KEEPING TEN OF THE ELEVEN CATEGORIES, BUT RELACING: 
Pinus nigra virus1 (Assigned=129, Summarized=129)
WITH:
Viruses (Assigned=18, Summarized=171)

THEN, in its own folder outputting the reads and seq ids of those unaccounted-for reads to add up to a total of 399753.

The eleven main categories of interest:

No hits (Assigned=57239, Summarized=57239)
Not assigned (Assigned=123886, Summarized=123886)
Viruses (Assigned=18, Summarized=171)
Bacteria (Assigned=852, Summarized=4637)
Amoebozoa (Assigned=12, Summarized=75)
Discoba (Assigned=2, Summarized=57)
Metazoa (Assigned=51, Summarized=4684)
Fungi (Assigned=2621, Summarized=52697)
Sar/Harosa (Assigned=68, Summarized=608)
Chlorophyta (Assigned=35, Summarized=739)
Embryophyta (Assigned=11738, Summarized=127769)

These add up to: 
total = 57239+123886+171+4637+75+57+4684+52697+608+739+127769 = 372562
percent accounted for = 372562/399753 = 0.93198
unaccounted for = 399753-372562 = 27191 

The categories containing unaccounted-for assigned reads (only interested in assigned reads):

Assigned reads in these parent categories, which are further up from our end nodes: 
    Opisthokonta (Assigned=579, Summarized=57960)
    Viridiplantae (Assigned=279, Summarized=129669)
    Eukaryota (Assigned=11251, Summarized=199620)
    cellular organisms (Assigned=13244, Summarized=217501)
    NCBI (Assigned=956, Summarized=399753)
Two categories in-line between Viridiplantae and Embryophta with assigned reads, so very hard to see in our tree:
    Streptophyta (Assigned=526, Summarized=128651)
    Streptophytina (Assigned=356, Summarized=128125)

Create the seq id files using, for example:
grep ">" reads-Sar-includeSummarized.fasta | sed 's/>\([^ ]*\) .*/\1/' > seqIds-Sar-includeSummarized.txt


===========
Hello, Ben:

Your running outcome looks great, they are what we want. Next step, please consider move onto whitebark pine, starting from Trinity first, then transdecode, and megan6. Thanks.

Jun-Jun

Steps:
- raw data is in four subfolders of /public/genomics/Data-Dec2020/WBP-NGS-rawData
- trimmomatic for quality control on each raw data file
- combine the files into a master R1 and R2 file for Trinity?
- will I need 1 TB or 2 TB memory for Trinity? Hopefully can use 2 TB. Will I need to split it into two subsets and run on Trinity separately?
- Trinity
- Transdecode
- BLASTp (using a different strategy?)
- MEGAN6

...see ass6-instructions-pt2.txt to continue...

