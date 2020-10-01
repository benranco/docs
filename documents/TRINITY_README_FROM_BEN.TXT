Trinity and all supporting software have already been installed on this image.

Trinity documentation is at:
http://trinityrnaseq.github.io
https://github.com/trinityrnaseq/trinityrnaseq/wiki

I followed these instructions to install Trinity:
https://github.com/trinityrnaseq/trinityrnaseq/wiki/Installing-Trinity

The Trinity executable used by default is in the 'trinityrnaseq-Trinity-v2.8.5'.
I have set up a TRINITY_HOME environment variable to point to this location. 
To run Trinity from a bash script, you can type $TRINITY_HOME/Trinity or 
${TRINITY_HOME}/Trinity, and it'll execute. An example bash script that does
this is the runMe.sh test script in:
/home/centos/trinityrnaseq-Trinity-v2.8.5/sample_data/test_Trinity_Assembly/

Do not delete the 'tools' folder or its contents. Some of these tools have been 
copied in /usr/local/bin, but some have symbolic links from /usr/local/bin 
back to this location.



For possible future reference, below are the steps I went through to install 
Trinity on this image:

===================================================
# First, I started with my CentOS7_SNPpipeline image that I had already 
# created. The steps I had used to make this image were:

sudo yum update
sudo yum install nano
sudo yum install git
sudo yum install wget
sudo yum install epel-release
sudo yum install R libcurl-devel.x86_64 libxml2-devel.x86_64

# these were necessary to install so that some required R packages could 
# be successfully installed:
sudo yum install openssl-devel
sudo yum install mariadb-devel

# from wherever you want the SNPpipeline to be downloaded to, enter:
git clone https://github.com/benranco/SNPpipeline.git

# this is where I'll mount an external volume:
mkdir /work
cd ~
ln -s /work work

# after setting the VM up the way you want it, exit, shut it down, and do 
# this to save it as an image:
openstack server image create --name image_name server_name

===================================================
# Second, the steps to make the CentOS7_Trinity image:

# To install Trinity for Kangakola (after above Centos setup), based on 
# these install instructions:
# https://github.com/trinityrnaseq/trinityrnaseq/wiki/Installing-Trinity

# starting with the setup from the CentOS7_SNPpipeline iamge, I then did:
sudo yum update
sudo yum install python-pip
sudo yum install python-devel
sudo yum groupinstall 'development tools'
sudo yum install numpy scipy
  # for the above I used these instructions:
  # https://linuxize.com/post/how-to-install-pip-on-centos-7/
  # https://stackoverflow.com/questions/47316178/python3-how-to-install-numpy-on-centos-7

# then I did a basic Trinity install:
wget -c https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.8.5.tar.gz
tar -xvzf Trinity-v2.8.5.tar.gz
cd trinityrnaseq-Trinity-v2.8.5/
sudo yum install cmake3
ls /usr/bin/cmake*
sudo ln -s /usr/bin/cmake3 /usr/bin/cmake
make
make plugins
sudo make install 
# the 'sudo make install' command copies everything to /usr/local/bin, but I actually 
# decided to set $TRINITY_HOME to the version install in home (see below).

nano ~/.bashrc
  # added this line to ./bashrc, which executes every time you boot your computer:
  export TRINITY_HOME=/home/centos/trinityrnaseq-Trinity-v2.8.5
# logged out, stopped and restarted the virtual machine, then logged back in.

# after this I tested Trinity at each of the following steps, installing what was required
# after each test run (which also was indicted in the installation instructions webpage above).
# To test Trinity, I went into the 
# /home/centos/trinityrnaseq-Trinity-v2.8.5/sample_data/test_Trinity_Assembly
# directory, and entered the command:
./runMe.sh
# Finally, after I had installed everything necessary to complete this command successfully, 
# I also ran:
./runMe.sh 1
# which indicated that I need to also install blat, which I installed (see below).

# install samtools, required by Trinity:
cd ~/tools/
wget -c https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -xjvf samtools-1.9.tar.bz2 
cd samtools-1.9
sudo yum install ncurses-devel # ./configure said this was required for samtools tview
./configure
make
sudo make install # copies things to /usr/local by default (executables to usr/local/bin, etc)

# install Jellyfish version 2:
cd ~/tools/
wget -c https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz
tar -xvzf jellyfish-2.2.10.tar.gz
cd jellyfish-2.2.10
./configure
make
sudo make install

# dowloaded bowtie linux binary from 
# https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.5.1/
cd ~/tools/
unzip bowtie2-2.3.5.1-linux-x86_64.zip
cd /usr/local/bin
sudo ln -s /home/centos/tools/bowtie2-2.3.5.1-linux-x86_64/bowtie2 bowtie2
sudo ln -s /home/centos/tools/bowtie2-2.3.5.1-linux-x86_64/bowtie2-build bowtie2-build
cd ~

# download salmon linux binary:
cd ~/tools/
wget -c https://github.com/COMBINE-lab/salmon/releases/download/v0.14.1/salmon-0.14.1_linux_x86_64.tar.gz
tar -xvzf salmon-0.14.1_linux_x86_64.tar.gz
cd /usr/local/bin
sudo ln -s /home/centos/tools/salmon-latest_linux_x86_64/bin/salmon salmon
cd ~

# download blat executable:
cd ~/tools
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat
chmod 555 blat
sudo mv blat /usr/local/bin/

# here's some info on blat:
http://genome.ucsc.edu/goldenPath/help/blatSpec.html
http://genome.ucsc.edu/FAQ/FAQblat.html
http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads
http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/
# if I ever need to compile blat from source rather than just downloading the executable, 
# here's a helpful blog:
http://andy09360639.blogspot.com/2017/06/how-to-install-blat-on-centos7.html

# after setting the VM up the way you want it, exit, shut it down, and do 
# this to save it as an image:
openstack server image create --name image_name server_name

[centos@ben-a ~]$ 

