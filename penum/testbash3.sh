#!/bin/bash

for ((seed=12; seed<=12; seed=seed+1))
do
	filebase="svpc-e1-d100-"
	fileend=".txt"
	filename=$filebase$seed$fileend
	logend=".log"
	logfile=$filebase$logend

	#export OMP_NUM_THREADS=32
	echo $filename
        cd /home/mb31simi/svpgenerator && ./generate_random --dim 100 --seed $seed --bit 1000 --bitexponent 1.0 &> /home/mb31simi/NewPenum/penum/tinput/$filename
	cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "tinput/$filename" --delta 0.9999 --beta 64 --sheight 10 --amax 2425
	#rm tinput/$filename
done
--
