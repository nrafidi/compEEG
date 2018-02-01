#!/bin/sh

subjectList=(AA BB CC DD EE F FF GG H HH I J JJ K L M N O P R S T U V W X Y Z)

winList=(1 2 3 4 5 6 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40)


for i in ${subjectList[@]}
do
for j in ${winList[@]}
do

jobName=${i}_${j}

Scripts/singlenode.pl ${jobName} "Scripts/matlab_batcher.sh \"addpath ~/compEEG/Classification; runCompClass_PermTest_cortex('${i}', ${j});\" \"${jobName}\""

done
done