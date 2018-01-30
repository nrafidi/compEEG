#!/bin/sh

subjectList=(AAA BBB CCC DDD EEE FFF GGG HHH III JJJ KKK MM OO PP QQ RR SS TT WW YY)

winList=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40)


for i in ${subjectList[@]}
do
for j in ${winList[@]}
do

jobName=${i}_${j}

Scripts/singlenode.pl ${jobName} "Scripts/matlab_batcher.sh \"addpath ~/compEEG/Classification; runCompClass_PermTest_cortex('${i}', ${j});\" \"${jobName}\""

done
done