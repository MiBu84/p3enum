#!/bin/bash


#Set the blocksize for BKZ of this testsuite
bsize=44
dim=88

for ((seed=5001; seed<=8000; seed=seed+1))
do
	filebase="enum"
	filebase=$filebase$dim
	strpart="-"$bsize"-"
	filebase=$filebase$strpart$seed	
	strpart=".sh"
	filebase=$filebase$strpart
	echo $filebase

	sbatch $filebase
done
