#!/bin/sh

subjectList=(AA AAA BB BBB CC CCC DD DDD EE EEE F FF FFF GG GGG H HH HHH I III J JJ JJJ K KKK L M MM N O OO P PP QQ R RR S SS T TT U V W WW X Y YY Z)

winList=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54)


for i in ${subjectList[@]}
do

jobName=${i}

Scripts/singlenode.pl ${jobName} "Scripts/matlab_batcher.sh \"addpath ~/compEEG/Classification; runCompClassification_SlidingFeat_Models_cortex('${i}');\" \"${jobName}\""

done