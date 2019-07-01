#!/bin/bash

inDir=$1

for i in $(seq 0 14);
do
    ./bin/Z_EE_Channel.exe Z_ee_fileList.txt $i 15 > unmergedOutputs/job_$i.log 2>&1 & 
done
