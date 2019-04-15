#!/bin/bash

bsize=44
rhour=0
rmin=15
dim=88

for ((seed=5001; seed<=8000; seed=seed+1))
do
	sed -e "s;SEED;$seed;g" -e "s;BSIZE;$bsize;g" -e "s;DIM;$dim;g" -e "s;MAXHOUR;$rhour;g" -e "s;MAXMINUTE;$rmin;g" template2.sh > batches/"enum$dim-$bsize-$seed".sh
done
