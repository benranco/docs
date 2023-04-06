THESE ARE SOME NOTES AND EMAIL EXCHANGES WITH JUN-JUN FROM ass6-pt6.txt and ass8-tree503.txt on using RSEM with the align_and_estimate_abundance.pl (packaged with Trinity). I've pasted the ass8-tree503 notes first, because they're shorter and gives a more consise treatment of how to use RSEM. The notes from ass6-pt6 will be useful if RSEM needs to be re-installed, or to understand more fully how to use RSEM. 

SEE ALSO THE align_and_estimate_abundance.pl-USAGE.txt FILE FOR MORE NOTES ON USING THE align_and_estimate_abundance.pl SCRIPT WITH BOWTIE2 AND RSEM.

I ALSO WROTE THE TWO R SCRIPTS (filter_RSEM_output-genes.R, filter_RSEM_output-isoforms.R) to search for the output files (RSEM.genes.results or RSEM.isoforms.results) within the output subdirectories (one for each sample) of the align_and_estimate_abundance.pl/RSEM run and combine the results of the bowtie2/RSEM read mapping process from all samples into a single file each for genes and isoforms, and then filter them by TPM >= 1 in at least one of the samples. I've also included total_expected_count, total_TPM and total_FPKM columns.





=====================================================================
=====================================================================
THESE NOTES WERE COPIED FROM ass8-tree503.txt, showing concisely how to use RSEM:

I've pasted my Trinity and Trimmomatic notes which precede it because their output was used as input, and I used the same tree503_input_samples.txt samples file from Trinity for the RSEM run.


==========

sudo mount -t nfs 10.20.0.6:/Public/genomics/junjun/Data-Dec2020/nanuq-2022/Rust-tree503-mRNA-seq-march2022 /work


--
Getting ready for Trinity run: 

Based on:  https://github.com/trinityrnaseq/trinityrnaseq/issues/885
The contents of the --samples_file (tree503_input_samples.txt):
H1      H1      ../step1-trim/NS.1805.003.NEBNext_dual_i7_177---NEBNext_dual_i5_177.H1_filtered_R_1P.fastq.gz ../step1-trim/NS.1805.003.NEBNext_dual_i7_177---NEBNext_dual_i5_177.H1_filtered_R_2P.fastq.gz
H2      H2      ../step1-trim/NS.1805.003.NEBNext_dual_i7_178---NEBNext_dual_i5_178.H2_filtered_R_1P.fastq.gz ../step1-trim/NS.1805.003.NEBNext_dual_i7_178---NEBNext_dual_i5_178.H2_filtered_R_2P.fastq.gz
H3      H3      ../step1-trim/NS.1805.003.NEBNext_dual_i7_179---NEBNext_dual_i5_179.H3_filtered_R_1P.fastq.gz ../step1-trim/NS.1805.003.NEBNext_dual_i7_179---NEBNext_dual_i5_179.H3_filtered_R_2P.fastq.gz
T1      T1      ../step1-trim/NS.1805.003.NEBNext_dual_i7_180---NEBNext_dual_i5_180.T1_filtered_R_1P.fastq.gz ../step1-trim/NS.1805.003.NEBNext_dual_i7_180---NEBNext_dual_i5_180.T1_filtered_R_2P.fastq.gz
T3      T3      ../step1-trim/NS.1805.003.NEBNext_dual_i7_181---NEBNext_dual_i5_181.T3_filtered_R_1P.fastq.gz ../step1-trim/NS.1805.003.NEBNext_dual_i7_181---NEBNext_dual_i5_181.T3_filtered_R_2P.fastq.gz
B1      B1      ../step1-trim/NS.1805.003.NEBNext_dual_i7_183---NEBNext_dual_i5_183.B1_filtered_R_1P.fastq.gz ../step1-trim/NS.1805.003.NEBNext_dual_i7_183---NEBNext_dual_i5_183.B1_filtered_R_2P.fastq.gz
B2      B2      ../step1-trim/NS.1805.003.NEBNext_dual_i7_184---NEBNext_dual_i5_184.B2_filtered_R_1P.fastq.gz ../step1-trim/NS.1805.003.NEBNext_dual_i7_184---NEBNext_dual_i5_184.B2_filtered_R_2P.fastq.gz
B3      B3      ../step1-trim/NS.1805.003.NEBNext_dual_i7_185---NEBNext_dual_i5_185.B3_filtered_R_1P.fastq.gz ../step1-trim/NS.1805.003.NEBNext_dual_i7_185---NEBNext_dual_i5_185.B3_filtered_R_2P.fastq.gz
U1      U1      ../step1-trim/NS.1805.003.NEBNext_dual_i7_186---NEBNext_dual_i5_186.U1_filtered_R_1P.fastq.gz ../step1-trim/NS.1805.003.NEBNext_dual_i7_186---NEBNext_dual_i5_186.U1_filtered_R_2P.fastq.gz
U4      U4      ../step1-trim/NS.1805.003.NEBNext_dual_i7_187---NEBNext_dual_i5_187.U4_filtered_R_1P.fastq.gz ../step1-trim/NS.1805.003.NEBNext_dual_i7_187---NEBNext_dual_i5_187.U4_filtered_R_2P.fastq.gz
U5      U5      ../step1-trim/NS.1805.003.NEBNext_dual_i7_188---NEBNext_dual_i5_188.U5_filtered_R_1P.fastq.gz ../step1-trim/NS.1805.003.NEBNext_dual_i7_188---NEBNext_dual_i5_188.U5_filtered_R_2P.fastq.gz
T4      T4      ../step1-trim/NS.1812.001.NEBNext_dual_i7_A1---NEBNext_dual_i5_A1.T4_filtered_R_1P.fastq.gz   ../step1-trim/NS.1812.001.NEBNext_dual_i7_A1---NEBNext_dual_i5_A1.T4_filtered_R_2P.fastq.gz


Real Trinity command: 
nohup /home/centos/software/trinity/trinityrnaseq-Trinity-v2.8.5/Trinity --seqType fq --samples_file tree503_input_samples.txt --max_memory 944G --CPU 64 --output trinityOut 2>&1 1>trinity_run.log &


--
Hi Jun-Jun,

Trimmomatic has finished running on the tree503 mRNA data, and I have begun running Trinity on the Trimmomatic paired reads output for all twelve samples.

I've attached some statistics from the Trimmomatic log file.

Ben
tree503_mRNA_trimmomatic_stats.csv


--
Hi Jun-Jun,

Trinity has completed successfully on the Limber Pine Trimmomatic output. Here is the output fasta file: 
https://www.dropbox.com/s/5v0s6l16kr1miyy/tree503-mRNAseqs-Trimm-Trinity.zip?dl=1


Here are some stats: 

84499649 / 1416200431 = 5.97% reads selected during normalization.
0 / 1416200431 = 0.00% reads discarded as likely aberrant based on coverage profiles.
0 / 1416200431 = 0.00% reads discarded as below minimum coverage threshold=1

Statistics:
===========
Trinity Version:      Trinity-v2.8.5
Compiler:             GCC
Trinity Parameters:   --seqType fq --samples_file tree503_input_samples.txt --max_memory 944G --CPU 64 --output trinityOut
Paired mode
 Input data
  Left.fasta    11952 MByte
  Right.fasta   11894 MByte
  Number of unique KMERs: 1010397781
  Number of reads:        168999298 Output data
  Trinity.fasta 521 MByte

--
Hello, Ben:
Great!

The 2nd step, the assembled transcripts will be filtered by gene-level expression estimates with transcripts per million (TPM) ≥ 1 in at least one of 12 samples using RSEM v1.3.3 (Li and Dewey, 2011).

The 3rd step, the assembled transcripts will be run using TransDecoder in the Trinity software package and filtered at a minimum protein length of 50.

Thanks.
--

===================
RSEM stuff for Tree503:

-------------
Running RSEM:

Based on reading through my notes in ass6-pt6.txt and align_and_estimate_abundance.pl-USAGE.txt:

I used the same parameters as the WWP18, WBP47, Hemlock30, LP36 runs, and using the same samples input file that I used for the Trinity run (in hindsight I should have copied the samples file and the Trinity.fasta output file into the RSEM folder first, since it creates files in the same directory): 
Usage for bowtie2/RSEM (with samples file, default params):
nohup perl $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts ../step2-trinity/tree503-mRNAseqs-Trimm-Trinity.fasta --seqType fq --samples_file ../step2-trinity/tree503_input_samples.txt   --est_method RSEM --aln_method bowtie2 --bowtie2_RSEM "--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 " --trinity_mode --prep_reference --thread_count 62 --output_dir outRSEM 2>&1 1> tree503_align_and_estimate_abundance-rsem_bowtie2_end-to-end.log &

--
Hi Jun-Jun, 

I've begun running RSEM/bowtie2 via the align_and_estimate_abundance.pl script that comes packaged with Trinity, using the same input parameters (the default suggestion for the --bowtie2_RSEM parameter string) for RSEM/bowtie2 as the WWP18, WBP47, Hemlock30 and LP36 datasets). I used the Trinity output fasta file and the Trimmomatic-processed read files as input.

Here is the command I used: 
perl $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts ../step2-trinity/tree503-mRNAseqs-Trimm-Trinity.fasta --seqType fq --samples_file ../step2-trinity/tree503_input_samples.txt   --est_method RSEM --aln_method bowtie2 --bowtie2_RSEM "--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 " --trinity_mode --prep_reference --thread_count 62 --output_dir outRSEM

Ben
--

Hi Jun-Jun,

RSEM is progressing, but I'd forgotten it takes so long. It is about a quarter of the way through, so I expect it to finish around the end of next week. Is there something you'd like me to focus on in the meantime?

Holly has asked me to help her with installing BUSCO and Translig on a Boreal Cloud image, so I've done a bit of work on that. I can also try again to understand the Linkage Disequilibrium paper, and investigate software that might do something similar either in R or Python.
--
You can try to run ‘Transdecoder’ on the Tree503 race assembly, as well as to pick up the longest isoform per gene.
--
[While waiting for RSEM to finish, I did some TransDecoder stuff, which is documented in its own section below under "TransDecoder stuff for Tree503".
--
After RSEM completed, I ran my scripts: 
filter_RSEM_output-genes.R
filter_RSEM_output-isoforms.R
which use this sort of command to find all the RSEM output files that I'll be working with:
find . -mindepth 2 -maxdepth 2 -name RSEM.isoforms.results -type f
find . -mindepth 2 -maxdepth 2 -name RSEM.genes.results -type f

First, make sure the RSEM process is all done:
grep -c "Expression Results are written" wwp18-align_and_estimate_abundance-rsem_bowtie2_end-to-end.log

Then created a folder RSEM_filtered.
Then run my filtering scripts:
nohup ./run_filter_RSEM.sh 2>&1 1>run.log &

--
Then created a backup of the results in RSEM_filtered before replacing all occurences of TRINITY with T503 in the output files (I should have done this with the original TRINITY output before running RSEM).
Then: 
sed -i 's/TRINITY/T503/g' *.*

scp -J borealremote brancourt@borealjump.nfis.org:/public/genomics/junjun/Data-Dec2020/nanuq-2022/Rust-tree503-mRNA-seq-march2022/step3-ReadMapping/RSEM_filtered/T503_RSEM_results_filtered.zip .


--
Hi Jun-Jun,

RSEM has also completed running on the Tree503 data. 
I've run my script to combine the results of the bowtie2/RSEM read mapping process from all samples into a single file each for genes and isoforms, and then filter them by TPM >= 1 in at least one of the samples. I've also included total_expected_count, total_TPM and total_FPKM columns. Here's a link to the results:
https://www.dropbox.com/s/j1ausy2wpqx787e/T503_RSEM_results_filtered.zip?dl=1

The .selected.tab files just show the selected genes or isoforms after combining and then filtering the results.

The .filtered.tab files are the combined results in tab-delimited format. I couldn't use comma-delimited because some of the fields contained commas. If you open it in Excel as a text or csv file, but choose tab as the delimiting field instead of comma, it should display properly. 

There are more selected isoform ids (125491) than gene ids (73279) after filtering, but more isoform ids were filtered out, and I remember from the previous datasets there were some selected gene ids that weren't found in the selected isoform ids. So that might also be the case this time. 
 
I ran my TPM filtering script on the RSEM.genes.results files and the RSEM.isoforms.results files separately, so if the TPM values were different in the genes vs isoforms files, or if the RSEM output recorded different gene selections in the genes vs isoforms files then there might have been a different set of ids selected after the filtering. Also, maybe the TMP values were generally lower per isoform than per gene, which might explain why more isoform ids were filtered out than gene ids. 

Ben
--






=====================================================================
=====================================================================
THESE NOTES WERE COPIED FROM ass6-pt6.txt, showing how I installed and learned how to use RSEM. Search for "Installing RSEM" for where the RSEM content begins.



===================================

---
Hello, Ben:

I like to filter the Trinity assembly further to decrease the total seq number of the assembly before running Megan6 based on the reads mapped back to assembly using CLC or Botwie. For examples, those seq with mapped reads less than one reads per one or 10 million of total mapped reads can be filtered out due to their low abundance inside the total sets, or using RPKM/FPKM at one.

For, using RPKM=1 as filter, we can cut WWP18 assembly from 675K sequence into 90K.

Thanks for your consideration.

---
Hello, Ben:

Please the process as below. I’m not sure if your original Trinity assembly process include this step or not. To ‘quantify read counts for each gene/isoform’ will allow us to filter out those genes/isoforms (the assembled sequences in output fasta file) with low numbers of reads in a samples. For example, we just keep those genes/isoforms with at least 10 mapped reads in at least one of total samples for MEGAN6 run.

Thanks.
Jun-Jun

https://southgreenplatform.github.io/trainings/trinityTrinotate/TP-trinity/

3.2 Reads mapping back rate and abundance estimation using the trinity script align_and_estimate_abundance.pl
Read congruency is an important measure in determining assembly accuracy. Clusters of read pairs that align incorrectly are strong indicators of mis-assembly. A typical “good” assembly has ~80% reads mapping to the assembly and ~80% are properly paired.

Several tools can be used to calculate reads mapping back rate over Trinity.fasta assembly : bwa, bowtie2 (mapping), kallisto, salmon (pseudo mapping). Quantify read counts for each gene/isoform can be calculate. Mapping and quantification can be obtained by using the –est_method argument into the align_and_estimate_abundance.pl script.

We will performing this analyses step successively with the align_and_estimate_abundance.pl script :

Pseudomapping methods (kallisto or salmon) are faster than mapping based. So firstly we will use salmon to pseudoalign reads from sample to the reference and quantify abondance.
Then, we will use bowtie2 and RSEM to align and quantify read counts for each gene/isoform.

---
Hi Jun-Jun,

I've read through the instructions in that link you sent me, and I think I understand what I need to do. It will require the read files (Trimmomatic output) for each dataset and the align_and_estimate_abundance.pl utility script that comes packaged with Trinity, along with the software it depends on. For me to run it on your Z800 in room 349 I would need to copy the Trimmomatic output data from the Boreal Cloud to your computer and install the Trinity package on your computer. It is quite complicated to install the Trinity package and its dependencies. I would prefer to use the Boreal Cloud where I have it already installed if possible, but I can install it on your computer if necessary. 

In the example on the website, they used Salmon first, which is much faster than Bowtie + RSEM. Is there a possibility Salmon would be good enough, or do we need to run bowtie2? 

Shall we have a chat about how to proceed? Feel free to call me when convenient.
---

Hi Jun-Jun,

Which fasta file would you like me to use for the Read Mapping? Should I use the -longest.fasta file after already having filtered it by longest ORF/sequence, or the original output of Trinity before filtering by longest ORF/sequence?
---
Please use “the original output of Trinity before filtering by longest ORF/sequence’
The longest ones may not have the highest numbers of mapped reads.
---
Also please check the parameters and setting for the read mapping, what are the default settings in bowtie2
---
See the link for settings
http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#setting-function-options

The bowtie2 aligner
End-to-end alignment versus local alignment

I prefer to use ‘local alignment’ with the default minimum score threshold is 20 + 8.0 * ln(L),
Preset options in --local mode
--very-fast-local       Same as: -D 5 -R 1 -N 0 -L 25 -i S,1,2.00
--fast-local            Same as: -D 10 -R 2 -N 0 -L 22 -i S,1,1.75
--sensitive-local       Same as: -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default in --local mode)
--very-sensitive-local  Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50

I prefer to use ‘--very-fast-local’ and ‘--very-sensitive-local’
Any others?

---
Thanks Jun-Jun, 

I am about to start an initial test run of the align_and_estimate_abundance.pl script using Salmon because it is much faster. After that, if I don't run into problems I'll try with RSEM + bowtie2. 

I think I might need to install RSEM first. I'll look at the RSEM parameters and bowtie2 parameters while I'm figuring this out (while the test with Salmon is running), and I'll let you know if I have questions.
---


------------

Trimmomatic data folders: 

/work/LP-NGS-RawData/LP-RNAseq-2015-trimmomatic
Cankered
Needle
Res
Sus

/work/Hemlock-Heterobasidion-2017-trimmomatic
Control
Cured
Pot-15
Wild

/work/WBP-NGS-trimmomatic
MJ9
Organ9
rust
SA9

/work/WWP-NGS-trimmomatic
HB
SB
H-Canker
H-SSfree
PN6

------------

First test with WWP18: 


Usage for salmon (with samples file):
nohup perl $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts trin-assemblies-renamedSeqIds/WWP18-nt675367-Trinity.fasta --seqType fq --samples_file wwp18-samples.tab   --est_method salmon --trinity_mode --prep_reference --thread_count 62 --output_dir outSalmonWWP18 2>&1 1> salmon_align_and_estimate_abundance-WWP18.log &


Usage for salmon:
nohup perl $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts trin-assemblies-renamedSeqIds/WWP18-nt675367-Trinity.fasta --seqType fq --left <string>   --right <string>   --est_method salmon --trinity_mode --prep_reference --thread_count 62 --output_dir outSalmonWWP18 2>&1 1> salmon_align_and_estimate_abundance-WWP18.log &


Usage for salmon:
nohup perl $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts trin-assemblies-renamedSeqIds/WWP18-10000.fasta --seqType fq --left WWP18/H-Canker/S0008EF_GCCAAT_L003_filtered_R_1P.fq,WWP18/H-Canker/S0008F2_CAGATC_L003_filtered_R_1P.fq,WWP18/H-Canker/S0008F6_CTTGTA_L003_filtered_R_1P.fq,WWP18/H-SSfree/S0008F1_ACTTGA_L003_filtered_R_1P.fq,WWP18/H-SSfree/S0008F7_TTAGGC_L003_filtered_R_1P.fq,WWP18/H-SSfree/S0008F9_ATCACG_L003_filtered_R_1P.fq,WWP18/SB/S0008F3_GGCTAC_L002_filtered_R_1P.fq,WWP18/SB/S0008F5_GATCAG_L002_filtered_R_1P.fq,WWP18/SB/S0008F8_TAGCTT_L002_filtered_R_1P.fq,WWP18/HB/S0008F0_ACAGTG_L002_filtered_R_1P.fq,WWP18/HB/S0008F4_CGATGT_L002_filtered_R_1P.fq,WWP18/HB/S0008FA_TGACCA_L002_filtered_R_1P.fq,WWP18/PN6/A03_N6ABXX_4_filtered_1P.fastq,WWP18/PN6/A04_N6ABXX_4_filtered_1P.fastq,WWP18/PN6/A05_N6ABXX_5_filtered_1P.fastq,WWP18/PN6/A06_N6ABXX_6_filtered_1P.fastq,WWP18/PN6/A07_U8ABXX_2_filtered_1P.fastq,WWP18/PN6/A08_GTABXX_2_filtered_1P.fastq   --right WWP18/H-Canker/S0008EF_GCCAAT_L003_filtered_R_2P.fq,WWP18/H-Canker/S0008F2_CAGATC_L003_filtered_R_2P.fq,WWP18/H-Canker/S0008F6_CTTGTA_L003_filtered_R_2P.fq,WWP18/H-SSfree/S0008F1_ACTTGA_L003_filtered_R_2P.fq,WWP18/H-SSfree/S0008F7_TTAGGC_L003_filtered_R_2P.fq,WWP18/H-SSfree/S0008F9_ATCACG_L003_filtered_R_2P.fq,WWP18/SB/S0008F3_GGCTAC_L002_filtered_R_2P.fq,WWP18/SB/S0008F5_GATCAG_L002_filtered_R_2P.fq,WWP18/SB/S0008F8_TAGCTT_L002_filtered_R_2P.fq,WWP18/HB/S0008F0_ACAGTG_L002_filtered_R_2P.fq,WWP18/HB/S0008F4_CGATGT_L002_filtered_R_2P.fq,WWP18/HB/S0008FA_TGACCA_L002_filtered_R_2P.fq,WWP18/PN6/A03_N6ABXX_4_filtered_2P.fastq,WWP18/PN6/A04_N6ABXX_4_filtered_2P.fastq,WWP18/PN6/A05_N6ABXX_5_filtered_2P.fastq,WWP18/PN6/A06_N6ABXX_6_filtered_2P.fastq,WWP18/PN6/A07_U8ABXX_2_filtered_2P.fastq,WWP18/PN6/A08_GTABXX_2_filtered_2P.fastq   --est_method salmon --trinity_mode --prep_reference --thread_count 62 --output_dir outSalmonWWP18 2>&1 1> salmon_align_and_estimate_abundance-WWP18.log &



--left parameter content WWP18: 

WWP18/H-Canker/S0008EF_GCCAAT_L003_filtered_R_1P.fq,WWP18/H-Canker/S0008F2_CAGATC_L003_filtered_R_1P.fq,WWP18/H-Canker/S0008F6_CTTGTA_L003_filtered_R_1P.fq,WWP18/H-SSfree/S0008F1_ACTTGA_L003_filtered_R_1P.fq,WWP18/H-SSfree/S0008F7_TTAGGC_L003_filtered_R_1P.fq,WWP18/H-SSfree/S0008F9_ATCACG_L003_filtered_R_1P.fq,WWP18/SB/S0008F3_GGCTAC_L002_filtered_R_1P.fq,WWP18/SB/S0008F5_GATCAG_L002_filtered_R_1P.fq,WWP18/SB/S0008F8_TAGCTT_L002_filtered_R_1P.fq,WWP18/HB/S0008F0_ACAGTG_L002_filtered_R_1P.fq,WWP18/HB/S0008F4_CGATGT_L002_filtered_R_1P.fq,WWP18/HB/S0008FA_TGACCA_L002_filtered_R_1P.fq,WWP18/PN6/A03_N6ABXX_4_filtered_1P.fastq,WWP18/PN6/A04_N6ABXX_4_filtered_1P.fastq,WWP18/PN6/A05_N6ABXX_5_filtered_1P.fastq,WWP18/PN6/A06_N6ABXX_6_filtered_1P.fastq,WWP18/PN6/A07_U8ABXX_2_filtered_1P.fastq,WWP18/PN6/A08_GTABXX_2_filtered_1P.fastq

--right parameter content WWP18: 

WWP18/H-Canker/S0008EF_GCCAAT_L003_filtered_R_2P.fq,WWP18/H-Canker/S0008F2_CAGATC_L003_filtered_R_2P.fq,WWP18/H-Canker/S0008F6_CTTGTA_L003_filtered_R_2P.fq,WWP18/H-SSfree/S0008F1_ACTTGA_L003_filtered_R_2P.fq,WWP18/H-SSfree/S0008F7_TTAGGC_L003_filtered_R_2P.fq,WWP18/H-SSfree/S0008F9_ATCACG_L003_filtered_R_2P.fq,WWP18/SB/S0008F3_GGCTAC_L002_filtered_R_2P.fq,WWP18/SB/S0008F5_GATCAG_L002_filtered_R_2P.fq,WWP18/SB/S0008F8_TAGCTT_L002_filtered_R_2P.fq,WWP18/HB/S0008F0_ACAGTG_L002_filtered_R_2P.fq,WWP18/HB/S0008F4_CGATGT_L002_filtered_R_2P.fq,WWP18/HB/S0008FA_TGACCA_L002_filtered_R_2P.fq,WWP18/PN6/A03_N6ABXX_4_filtered_2P.fastq,WWP18/PN6/A04_N6ABXX_4_filtered_2P.fastq,WWP18/PN6/A05_N6ABXX_5_filtered_2P.fastq,WWP18/PN6/A06_N6ABXX_6_filtered_2P.fastq,WWP18/PN6/A07_U8ABXX_2_filtered_2P.fastq,WWP18/PN6/A08_GTABXX_2_filtered_2P.fastq


------------
Installing RSEM: 

https://github.com/bli25broad/RSEM_tutorial
https://github.com/deweylab/RSEM

wget https://github.com/deweylab/RSEM/archive/refs/tags/v1.3.3.tar.gz
mv v1.3.3.tar.gz RSEM-v1.3.3.tar.gz
tar xvzf RSEM-v1.3.3.tar.gz 
cd RSEM-1.3.3/
make
make ebseq
pwd
nano ~/.bash_profile   ## edit the PATH environment variable to include RSEM directory, then log out & log back in

When I first tested RSEM with the align_and_estimate_abundance.pl script, I got this error:
Can't locate Env.pm in @INC (@INC contains: /home/centos/software/RSEM-1.3.3 /usr/local/lib64/perl5 /usr/local/share/perl5 /usr/lib64/perl5/vendor_perl /usr/share/perl5/vendor_perl /usr/lib64/perl5 /usr/share/perl5 .) at /home/centos/software/RSEM-1.3.3/rsem-prepare-reference line 10.

I found this solution: 
https://superuser.com/questions/1181310/perl-script-cant-locate-env-pm-in-inc
sudo yum install 'perl(Env)'

------------
Tests of RSEM/bowtie2 (actual runs in the next ==== section below: 

First test of bowtie2/RSEM with WWP18: 

ACTUAL RUN: Usage for bowtie2/RSEM (with samples file, default params):
nohup perl $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts trin-assemblies-renamedSeqIds/WWP18-nt675367-Trinity.fasta --seqType fq --samples_file wwp18-samples.txt   --est_method RSEM --aln_method bowtie2 --bowtie2_RSEM "--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 " --trinity_mode --prep_reference --thread_count 62 --output_dir outRSEM 2>&1 1> align_and_estimate_abundance-rsem_bowtie2_end-to-end.log &


Usage for bowtie2/RSEM (with samples file, default params):
nohup perl $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts trin-assemblies-renamedSeqIds/WWP18-1000.fasta --seqType fq --samples_file wwp18-samples.txt   --est_method RSEM --aln_method bowtie2 --bowtie2_RSEM "--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 " --trinity_mode --prep_reference --thread_count 62 --output_dir outRSEM 2>&1 1> bowtie-end-to-end-1000-rsem_align_and_estimate_abundance.log &


Usage for bowtie2/RSEM (with new samples file, custom params very-fast-local):
nohup perl $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts trin-assemblies-renamedSeqIds/WWP18-10000.fasta --seqType fq --samples_file wwp18-samples.txt   --est_method RSEM --aln_method bowtie2 --bowtie2_RSEM "--no-mixed --no-discordant --gbar 1000 --very-fast-local -k 200 " --trinity_mode --prep_reference --thread_count 62 --output_dir outRSEM 2>&1 1> bowtie2-customParams-rsem_align_and_estimate_abundance.log &

Usage for bowtie2/RSEM (with new samples file, custom params very-sensitive-local):
nohup perl $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts trin-assemblies-renamedSeqIds/WWP18-1000.fasta --seqType fq --samples_file wwp18-samples.txt   --est_method RSEM --aln_method bowtie2 --bowtie2_RSEM "--no-mixed --no-discordant --gbar 1000 --very-sensitive-local -k 200 " --trinity_mode --prep_reference --thread_count 62 --output_dir outRSEM 2>&1 1> bowtie2-customParams-rsem_align_and_estimate_abundance.log &


This seems to be the bowtie2 parameters it's using by default?:
CMD: set -o pipefail && bowtie2 --no-mixed --no-discordant --gbar 1000 --end-to-end -k 200  -q -X 800 -x /work/Bens-Trinity-assembly/ReadMappingCounts/trin-assemblies-renamedSeqIds/WWP18-10000.fasta.bowtie2 -1 /work/Bens-Trinity-assembly/ReadMappingCounts/WWP18/H-Canker/S0008EF_GCCAAT_L003_filtered_R_1P.fq -2 /work/Bens-Trinity-assembly/ReadMappingCounts/WWP18/H-Canker/S0008EF_GCCAAT_L003_filtered_R_2P.fq -p 62 | samtools view -@ 62 -F 4 -S -b | samtools sort -@ 62 -n -o bowtie2.bam

CMD: touch RSEM.isoforms.results.ok
CMD: set -o pipefail && bowtie2 --no-mixed --no-discordant --gbar 1000 --end-to-end -k 200  -q -X 800 -x /work/Bens-Trinity-assembly/ReadMappingCounts/trin-assemblies-renamedSeqIds/WWP18-10000.fasta.bowtie2 -1 /work/Bens-Trinity-assembly/ReadMappingCounts/WWP18/H-Canker/S0008F2_CAGATC_L003_filtered_R_1P.fq -2 /work/Bens-Trinity-assembly/ReadMappingCounts/WWP18/H-Canker/S0008F2_CAGATC_L003_filtered_R_2P.fq -p 62 | samtools view -@ 62 -F 4 -S -b | samtools sort -@ 62 -n -o bowtie2.bam 


============
Actual runs of bowtie2/RSEM (with samples file, default params):

WWP18 ACTUAL RUN: Usage for bowtie2/RSEM (with samples file, default params):
nohup perl $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts trin-assemblies-renamedSeqIds/WWP18-nt675367-Trinity.fasta --seqType fq --samples_file wwp18-samples.txt   --est_method RSEM --aln_method bowtie2 --bowtie2_RSEM "--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 " --trinity_mode --prep_reference --thread_count 62 --output_dir outRSEM 2>&1 1> align_and_estimate_abundance-rsem_bowtie2_end-to-end.log &

WBP47 ACTUAL RUN: Usage for bowtie2/RSEM (with samples file, default params):
nohup perl $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts trin-assemblies-renamedSeqIds/WBP47-nt813669-Trinity.fasta --seqType fq --samples_file wbp47-samples.txt   --est_method RSEM --aln_method bowtie2 --bowtie2_RSEM "--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 " --trinity_mode --prep_reference --thread_count 30 --output_dir wbp47-outRSEM 2>&1 1> wbp47-align_and_estimate_abundance-rsem_bowtie2_end-to-end.log &

Hemlock30 ACTUAL RUN: Usage for bowtie2/RSEM (with samples file, default params):
nohup perl $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts trin-assemblies-renamedSeqIds/Hemlock30-nt1193147-Trinity.fasta --seqType fq --samples_file hemlock30-samples.txt   --est_method RSEM --aln_method bowtie2 --bowtie2_RSEM "--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 " --trinity_mode --prep_reference --thread_count 30 --output_dir hemlock30-outRSEM 2>&1 1> hemlock30-align_and_estimate_abundance-rsem_bowtie2_end-to-end.log &

LP36 ACTUAL RUN: Usage for bowtie2/RSEM (with samples file, default params):
nohup perl $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts trin-assemblies-renamedSeqIds/LP36-nt884751-Trinity.fasta --seqType fq --samples_file lp36-samples.txt   --est_method RSEM --aln_method bowtie2 --bowtie2_RSEM "--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 " --trinity_mode --prep_reference --thread_count 30 --output_dir lp36-outRSEM 2>&1 1> lp36-align_and_estimate_abundance-rsem_bowtie2_end-to-end.log &



------------
Hi Jun-Jun, 

I've installed RSEM and I'm testing the align_and_estimate_abundance.pl script with RSEM/bowtie2, just using the default parameters for now. Bowtie2 is used by RSEM. So far it seems to be running smoothly.

I'll be looking at the parameter options for both RSEM and bowtie2 along with the ones you sent me for bowtie2.

If you'd like to look at RSEM as well here is some documentation:
Official RSEM readme:
https://github.com/deweylab/RSEM#table-of-contents
Tutorial of RSEM: 
https://github.com/bli25broad/RSEM_tutorial

The align_and_estimate_abundance.pl script gives an option to pass a string of parameters into RSEM. I've pasted the defaults below along with a couple other paramters for RSEM. I think this --bowtie2_RSEM string is where we'd include all bowtie2 parameters.

########################################
#  RSEM opts:
#
#  --bowtie_RSEM <string> :
#         if using 'bowtie', default: "--all --best --strata -m 300 --chunkmbs 512"
#
#  --bowtie2_RSEM <string> :
#         if using 'bowtie2', default: "--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 "
#
#  ** if you change the defaults, specify the full set of parameters to use! **
#
#  --include_rsem_bam              provide the RSEM enhanced bam file including posterior probabilities of read assignments.
#  --rsem_add_opts <string>        additional parameters to pass on to rsem-calculate-expression
#
########################################

I have a question for you. The align_and_estimate_abundance.pl script wants me to pass the input sample filenames to it in a tab-delimited samples.txt file with four parameters per line. The third and fourth parameters are the R1 and R2 filenames, but I'm not sure what to do for the first two, so I've just been putting a piece of the filename for both of those. I've attached my wwp18-samples.tab file for you to look at. Can you let me know if I should be specifying the first two parameters differently? Here's what the align_and_estimate_abundance.pl documentation says about it: 

tab-delimited text file indicating biological replicate relationships. ex:
cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq

Thanks,

Ben
------
Hello, Ben:

RSEM is a further run step on the basis of Bowtie-2. Bowtie2 only calculates read numbers while RSEM wants to compare samples based on the read numbers to identify genes (sequences) with significant differences among samples.

In the sample filenames with four parameters  , the 1st is for sample types, the 2nd is for biological replicates.

For WWP18, we have 18 samples: PN-6, HB-3, SB-3, HC-3, and HF-3, representing 5 types of samples (conditions), with replicates at 6, 3, 3, 3, and 3 tomes, respectively.

I'll get this table for you for three other species. Let us have a discussion if you have additional questions.

Jun-Jun
-----
Hi Jun-Jun, 

I've updated the wwp18-samples file (see attached), based on what you said. Can you tell me if this is correct?

It looks like there is no option with the align_and_estimate_abundance.pl script to run bowtie2 without RSEM. It calls RSEM, and RSEM then uses bowtie2. Is that okay with you? Would RSEM take a long time? Unless RSEM takes a really long time it might be faster to just use this script rather than writing my own script to call bowtie2 directly.

Here are the default bowtie2 parameters used by the script: 
"--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 "

Based on your preferred parameters, I will change this to: 
"--no-mixed --no-discordant --gbar 1000 --very-fast-local -k 200 "
or 
"--no-mixed --no-discordant --gbar 1000 --very-sensitive-local -k 200 "

Should I try it first with the --very-fast-local setting or the --very-sensitive-local setting?

Do these parameters look okay to you? I've pasted below the descriptions of the other parameters: 

--no-mixed
By default, when bowtie2 cannot find a concordant or discordant alignment for a pair, it then tries to find alignments for the individual mates. This option disables that behavior.

--no-discordant
By default, bowtie2 looks for discordant alignments if it cannot find any concordant alignments. A discordant alignment is an alignment where both mates align uniquely, but that does not satisfy the paired-end constraints (--fr/--rf/--ff, -I, -X). This option disables that behavior.

--gbar <int>
Disallow gaps within <int> positions of the beginning or end of the read. Default: 4.

-k mode: search for one or more alignments, report each
In -k mode, Bowtie 2 searches for up to N distinct, valid alignments for each read, where N equals the integer specified with the -k parameter. That is, if -k 2 is specified, Bowtie 2 will search for at most 2 distinct alignments. It reports all alignments found, in descending order by alignment score. The alignment score for a paired-end alignment equals the sum of the alignment scores of the individual mates. Each reported read or pair alignment beyond the first has the SAM 'secondary' bit (which equals 256) set in its FLAGS field. Supplementary alignments will also be assigned a MAPQ of 255. See the SAM specification for details.
Bowtie 2 does not "find" alignments in any specific order, so for reads that have more than N distinct, valid alignments, Bowtie 2 does not guarantee that the N alignments reported are the best possible in terms of alignment score. Still, this mode can be effective and fast in situations where the user cares more about whether a read aligns (or aligns a certain number of times) than where exactly it originated.


------------

"--no-mixed --no-discordant --gbar 1000 --very-sensitive-local -k 200 "

bowtie2 -p 62 -x ref.fasta -1 R1.fq -2 R2.fq --very-sensitive-local --met-file bowtie2-metrics.txt -S output.sam

bowtie2-build /work/Bens-Trinity-assembly/ReadMappingCounts/trin-assemblies-renamedSeqIds/WWP18-10000.fasta /work/Bens-Trinity-assembly/ReadMappingCounts/trin-assemblies-renamedSeqIds/WWP18-10000.fasta.bowtie2

nohup bowtie2 --no-mixed --no-discordant --gbar 1000 --very-fast-local -k 200  -q -X 800 -x /work/Bens-Trinity-assembly/ReadMappingCounts/trin-assemblies-renamedSeqIds/WWP18-10000.fasta.bowtie2 -1 /work/Bens-Trinity-assembly/ReadMappingCounts/WWP18/H-Canker/S0008EF_GCCAAT_L003_filtered_R_1P.fq -2 /work/Bens-Trinity-assembly/ReadMappingCounts/WWP18/H-Canker/S0008EF_GCCAAT_L003_filtered_R_2P.fq -p 62 -S bowtie2-output-S0008EF_GCCAAT.sam 2>&1 1> bowtie2-S0008EF_GCCAAT.log &

------------
============
AFTER successfully finishing the align....pl script for WWP18:
---
Hi Jun-Jun,

The WWP18 bowtie2/RSEM read mapping process has finished, but the others are still running.

I've attached a summary of the alignment rates for WWP18 in wwp18-overall-alignment-rate.txt.

There is an output folder for each sample/rep (18 in total). I've attached a zip file containing all the output files except the .bam files (which are too large) from one of these folders (HB-3-rep1-S0008F0_ACAGTG). Here's a list of the complete contents of the output folder for sample HB-3-rep1-S0008F0_ACAGTG of WWP18: 

bowtie2.bam  (4 GB)
bowtie2.bam.for_rsem.bam  (4 GB)
bowtie2.bam.ok
RSEM.genes.results
RSEM.isoforms.results
RSEM.isoforms.results.ok
RSEM.stat/
RSEM.stat/RSEM.cnt
RSEM.stat/RSEM.model
RSEM.stat/RSEM.theta

In your earlier instructions you said: "To ‘quantify read counts for each gene/isoform’ will allow us to filter out those genes/isoforms (the assembled sequences in output fasta file) with low numbers of reads in a samples. For example, we just keep those genes/isoforms with at least 10 mapped reads in at least one of total samples for MEGAN6 run."

I'm not sure which file to extract this information from, or how to find it. Do I need to search the BAM file, or does one of the RSEM output files (attached in the .zip file) contain the information?

I understand that my next step will be to go through all the sequences in the Trinity .fasta output file, and, for each sequence, check the number of mapped reads in each of the 18 samples. If the sequence has at least 10 mapped reads in at least one of these samples, keep it.
Thanks,

Ben

------------
The primary output of RSEM consists of two files, one for isoform-level estimates, and the other for gene-level estimates. Abundance estimates are given in terms of two measures. The first is an estimate of the number of fragments that are derived from a given isoform or gene.

outputs should be in these two files, but we need to put them in a table to review.

RSEM.genes.results
RSEM.isoforms.results
---
Yes. The data we want are shown in those two files at gene or isoform level. for each sample, there are three values/column
'expected_count',	'TPM',	'FPKM'. Please combine all 18 samples together in one table based on gene IDs, or isoform IDs, and then add three more column as 'total expected-count', 
'total-TPM', and 'total FPKM' per gene or per isoform by summing the values of all 18 samples. Thirdly, we can filter gene/isoform numbers based on one of these three values.

Great progress!

Have a nice weekend!
---
Hi Jun-Jun, 
If I put all 18 samples together in one table for RSEM.gene.results and one table for RSEM.isoforms.results, it will create two very large text/csv files, probably more than 700 MB each. I think it would be quite hard to open those for viewing in most programs. I'll begin writing a script to do that, but I can also just make it save the total-expected-count, total-TPM and total_FPKM if you prefer. 

I can also write a script that searches these files without combining them if you want to use their information to filter genes/isoforms. Just let me know what criteria you'd like to use to filter the genes/isoforms.
---
Let us filter the data at TPM =, or > 1 in at least one of 18 samples.

We are more interested in gene model than isoform model. Gene model will be linked to the longest isoform tor represent the genes with multiple isoforms per genes.
---

inPath

My filtering logic:
for each gene .csv file
  save to curr_gene_id the gene_id col of only those rows whose TPM >= 1
  write out the filtered csv file
  unique(c(master_gen_id,curr_gene_id))

--
Mods:
- filter individually only to get the master list of gene ids, and individual filtered output file
- merge the full tables
- then filter by the master list of gene ids
Then:
- add three more column as 'total expected-count', 'total-TPM', and 'total FPKM'


------------
Make sure the RSEM process is all done:
grep -c "Expression Results are written" wwp18-align_and_estimate_abundance-rsem_bowtie2_end-to-end.log

Then run my filtering scripts:
nohup ./run_filter_RSEM.sh 2>&1 1>run.log &

---
Hi Jun-Jun,

I've written and run a script to combine the results of the bowtie2/RSEM read mapping process from all samples into a single file each for genes and isoforms, and then filter them by TPM >= 1 in at least one of the samples. I've also included total_expected_count, total_TPM and total_FPKM columns. Here's a link to the results for WWP18, WBP47, LP36, Hemlock30:
https://www.dropbox.com/s/ycg8yw8w6fywjbr/RSEM_results_filtered.zip?dl=1

The .txt files just show the selected genes or isoforms after combining and then filtering the results.

The .tab files are the combined results in tab-delimited format. I couldn't use comma-delimited because some of the fields contained commas. If you open it in Excel as a text or csv file, but choose tab as the delimiting field instead of comma, it should display properly. At least, that's how it works in OpenOffice/LibreOffice.

These .tab files are all quite large, so they'll take a long time to load, and you'll need to make sure you have a decent amount of free RAM. You might want to start with the WWP18 genes file, which is the smallest at around 208 MB.

Let me know if you'd like me to make any changes to how I've done this.

Ben

---
---
Hello, Ben:

Your account issue is solved.

BTW, I reviewed the Hemlock data about the filtering transcript by TPM=1 and then picking up longest transcript (IDs with i) per gene (IDs with g). You get 90,730 seqs as the longest with TPM at lest 1 in at least one of 30 samples.

When I go back to see  how many transcripts (isoforms with i IDs) for 90,730 g IDs, you file “Hemlock30ALLsamples.RSEM.isoform.results.filtered.tab” contained only 123,545 transcripts (isoforms with I IDs), which corresponding to only 82,929 longest g IDs, different the 90,730 g IDs from another file.  Please have a confirmation.

Thanks!
---
Hi Jun-Jun, 

I ran my filtering script on the RSEM.genes.results files and the RSEM.isoforms.results files separately, so if the TPM values were different in the genes vs isoforms files, or if the RSEM output recorded different gene selections in the genes vs isoforms files then there might have been a different set of ids selected after the filtering.

Maybe the TMP values were generally lower per isoform than per gene? If this is so then that would explain why the filtered results were 90,730 gene ids but only 123,545 isoform ids with only 82,929 of the selected gene ids accounted for. Is this what you would expect?

I've confirmed that the Hemlock30_ALLsamples.RSEM.genes.results.filtered.tab file has 90730 gene ids and the Hemlock30_ALLsamples.RSEM.isoforms.results.filtered.tab file has 123546 isoform ids.

Ben

---

