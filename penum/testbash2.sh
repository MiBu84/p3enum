#!/bin/bash

for ((seed=776; seed<=4000; seed=seed+1))
do
	filebase="svpc-e1-d75-"
	fileend=".txt"
	filename=$filebase$seed$fileend
	logend=".log"
	logfile=$filebase$logend

	#export OMP_NUM_THREADS=32
	echo $filename
        cd /home/mb31simi/svpgenerator && ./generate_random --dim 75 --seed $seed --bit 750 --bitexponent 1.0 &> /home/mb31simi/NewPenum/penum/tinput/$filename
	cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "tinput/$filename" --delta 0.9999 --beta 40 --sheight 10 --amax 2055 >> $logfile
	#rm tinput/$filename
done
--
