#!/bin/bash

for ((seed=1; seed<=500; seed=seed+1))
do
	filebase="svpc-e1-d92-"
	fileend=".txt"
	filename=$filebase$seed$fileend
	logend=".log"
	logfile=$filebase$logend

	export OMP_NUM_THREADS=32
	echo $filename
	cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "tinput/$filename" --delta 0.999 --beta 52 --sheight 10 --amax --amax 2442 --trials 5  >> $logfile
done
