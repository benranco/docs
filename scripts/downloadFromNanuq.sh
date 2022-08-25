#!/bin/bash

# This is a comment because it begins with a #

# If you want to run a bash script, you either run it like this:
#    bash exampleBashScript.sh
# Or, if you have set the permissions of the file to be executable (you can do 
# that with this command: chmod 777 exampleBashScript.sh), you can run it like this: 
#    ./exampleBashScript.sh
#
# If this script will take a long time to finish and you want to run it using
# nohup so that it doesn't terminaate if the ssh session times out, you don't 
# need to use nohup on any of the commands inside the script, but when you run the
# script, run it like this:
#    nohup bash exampleBashScript.sh 2>&1 1>log.txt &
# Or like this: 
#    nohup ./exampleBashScript.sh 2>&1 1>log.txt &
#
# The 2>&1 means to pipe the stderr output to stdout.
# The 1>log.txt means to pipe the stdout to the log.txt file.
# The final & at the end means to run the script in the background so you can use the 
# terminal for other commands.
#
# You can then check on the contents of log.txt as it is being written to by using: 
#    tail -f log.txt
#
# The echo "hello" commands represent commands you can replace with your own.
# The "====" commands are separators so you can more easily find where one command ends
# and another begins in the log.txt output.


echo "hello 1"
echo "======================END OF COMMAND==========================="
echo ""

echo "hello 2"
echo "======================END OF COMMAND==========================="
echo ""

echo "hello 3"
echo "======================END OF COMMAND==========================="
echo ""

echo "hello 4"
echo "======================END OF COMMAND==========================="
echo ""

echo "hello 5"
echo "======================END OF COMMAND==========================="
echo ""


echo "Running the command to download stuff from Nanuq."

login="jjliu" # make sure this is your Nanuq login id
password="p9u6q6qd" # put your password between the quotes

##read -p "Login: " login && read -p "Password: " -s password && ## THIS PART NOT NEEDED IF USING THIS FILE
echo -n "j_username=$login&j_password=$password" > .auth.txt && chmod 600 .auth.txt && wget -O - "https://ces.genomequebec.com/nanuqMPS/readsetList?projectId=21146&tech=NovaSeq" --no-cookies --no-check-certificate --post-file .auth.txt | wget --no-cookies --no-check-certificate --post-file .auth.txt -ci -; rm -f .auth.txt

echo "======================END OF COMMAND==========================="
echo ""


## To confirm the download with the md5sum: md5sum -c readSets.md5



