#!/bin/bash

echo $0
echo $1
echo $2

dd="there"
echo "Processing" $dd "---" `date`
date

freeCpus=4

ncore="$(( ($(grep -c ^processor /proc/cpuinfo) ) ))"
echo $ncore

ncore="$(( ($(grep -c ^processor /proc/cpuinfo) )/2 ))"
echo $ncore

ncore="$(( ($(grep -c ^processor /proc/cpuinfo) - $freeCpus)/2 ))"
echo $ncore

