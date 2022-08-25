This  README file was created by Ben to help someone get started with BUSCO.
I haven't actually used BUSCO myself, I've only installed it.

==========
BUSCO website: 
https://busco.ezlab.org/

BUSCO userguide: 
https://busco.ezlab.org/busco_userguide.html


I installed BUSCO by first installing miniconda to /home/centos/software/miniconda, 
then creating a conda environment called busco_env, then installing the busco 
package (which also installs its dependencies) to that environment.

One thing to note is that the conda package dependencies required by BUSCO take 
up a lot of space, so there's only about 2 GB of free space left on this VM image.
So make sure you mount an external folder or drive to hold your data.
 
This system is set up to automatically activate the conda busco_env at start-up, so 
it might take a little longer to launch your VM than normal. But to run BUSCO, you 
shouldn't need to do anything more than follow the instructions in the BUSCO 
userguide.

==========
IMPORTANT: 

If you want to run a BUSCO process that will take a long time and you're afraid your 
ssh session might time out, then you can launch a "screen" session first, and then 
run BUSCO from within that. "screen" is a program that lets you create a terminal 
session that will stay active even if your ssh session terminates (you can use 
as an alternative to running a command via "nohup", for those who are familiar with
nohup).  

So you can create a screen session, start a process running within it, and it will 
keep running even if your ssh session terminates.

Here's how to use screen: 
from: https://linuxize.com/post/how-to-use-linux-screen/

To start a screen session called "mysession", use this command: 
    screen -S mysession
This will launch the new session. I've set things up so that the conda busco_env 
will be automatically activated in your screen session, so it might take a little
while to launch.
Once your screen session has launched, you can use it like a normal terminal 
session.
IMPORTANT: Don't use the "exit" command to leave your screen session unless you 
also want to terminate it and anything running it it. Instead, just detach from 
it.
To detach from a screen session you are currently in, use this command:
    Ctrl-a d
To see a list of screen commands from within a session, use this command: 
    Ctrl-a
To terminate a screen session from within the session, use: 
    exit
To re-attach to an active screen session, use: 
    screen -r sessionname
For example, using the session we created above: screen -r mysession
If you can't remember the name of your screen session, to see a list of active 
screen sessions, use this command: 
    screen -ls  

Those are basically all the screen commands you need to know.



==========================================
Below this line are some notes I made when installing Miniconda and BUSCO.
There are some helpful links for learning how to use conda, but other than that 
it shouldn't be of interest unless you need to re-install conda and BUSCO.



-----
Miniconda links: 

Anaconda vs Miniconda: 
https://www.educative.io/edpresso/anaconda-vs-miniconda
https://stackoverflow.com/questions/45421163/anaconda-vs-miniconda
- I will go with Miniconda, it seems if I need the Anacondan packages I can later do: conda install anaconda

How to install Miniconda on CentOS 7: 
https://deeplearning.lipingyang.org/2018/12/24/install-miniconda-on-centos-7-redhat-7/
https://docs.conda.io/en/latest/miniconda.html
https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html  
(also how to remove miniconda)
- it is best to install Miniconda for the local user, which does not require 
administrator permissions and is the most robust type of installation.

Anaconda FAQ: 
https://docs.anaconda.com/anaconda/user-guide/faq/

Commonly used conda commands: 
https://deeplearning.lipingyang.org/2018/12/25/conda-commands-create-virtual-environments-for-python-with-conda/

Conda cheatsheet: 
https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf

Conda searchable documentation: 
https://docs.conda.io/en/latest/

Editing the .condarc file (config): 
https://docs.conda.io/projects/conda/en/latest/user-guide/configuration/use-condarc.html

Difference between source activate myenv and conda activate myenv, and other questions:
https://quick-adviser.com/what-is-the-difference-between-conda-activate-and-source-activate/

Potential issues with using a conda environment within screen (shouldn't be a problem for us): 
https://stackoverflow.com/questions/50591901/screen-inside-the-conda-environment-doesnt-work
  
-----
Process to install miniconda: 

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
- then follow the prompts
- I chose to specify this as its install folder: /home/centos/software/miniconda3

After the installation script finished, it said: 
"==> For changes to take effect, close and re-open your current shell. <==
If you'd prefer that conda's base environment not be activated on startup, 
   set the auto_activate_base parameter to false: 
conda config --set auto_activate_base false
Thank you for installing Miniconda3!"

- I exited my shell and logged back in
- testing conda:
conda -V
conda update conda   # there were some updates
conda list

-----
Installing BUSCO: 

- in the below command the "-c bioconda", etc, specifies an additional channel to 
search for packages

BUSCO installation instructions: 
https://busco.ezlab.org/busco_userguide.html#conda-package

This command to install BUSCO to the base conda environment didn't work for some reason: 
conda install -c conda-forge -c bioconda busco=5.3.2
So I installed it to its own environment called busco_env, see below: 

conda create -n busco_env
conda install -n busco_env -c conda-forge -c bioconda busco=5.3.2
conda activate busco_env
conda deactivate
conda config --set auto_activate_base false
conda config --set auto_update_conda false
- there is no way in .condarc to auto-activate a different environment than base, 
so I put this line in the .bashrc file:
conda activate busco_env
- I then exited the VM, rebooted it, and logged back in to test. Also tested with screen.


