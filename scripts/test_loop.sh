#!/bin/bash

counter=0
ncore="$(( ($(grep -c ^processor /proc/cpuinfo) -4)/2 ))"


for i in `seq 1 13`
do
  if [ $counter -ge $ncore ]; then
    #wait
    echo "wait reset"
    counter=0
  fi 
  echo $i

  counter=$(( $counter + 1 ))
done
#wait

counter=0

echo "running parallel processing"

for f in $(ls /home/benrancourt/Desktop/grace/SNPpipeline-ben-2016-08-25/reporttemp | grep -v "filled.Rds")
do    
	if [ $counter -ge $ncore ]; then	
		#wait
		echo "wait reset"
		counter=0
	fi
                    echo $f
	counter=$(( $counter + 1 ))
done
#wait





