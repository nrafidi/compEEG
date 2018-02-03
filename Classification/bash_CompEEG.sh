#!/bin/sh

subjectList=(AA AAA BB BBB CC CCC DD DDD EE EEE F FF FFF GG GGG H HH HHH I III J JJ JJJ K KKK L M MM N O OO P PP QQ R RR S SS T TT U V W WW X Y YY Z)

winList=(1 2 3 4 5 6 41 42 43 44 45 46 47 48 49 50 51 52 53 54)


for i in ${subjectList[@]}
do
for j in ${winList[@]}
do

jobName=${i}_${j}

Scripts/singlenode.pl ${jobName} "Scripts/matlab_batcher.sh \"addpath ~/compEEG/Classification; runCompClass_PermTest_cortex('${i}', ${j});\" \"${jobName}\""

done
done