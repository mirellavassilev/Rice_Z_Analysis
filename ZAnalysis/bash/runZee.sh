#!/bin/bash

inDir=$1

for i in $(seq 0 11);
do
    ./bin/Z_EE_Channel.exe Z_ee_fileList_v3_June14.txt $i 12 > unmergedOutputs/job_$i.log 2>&1 & 
done
