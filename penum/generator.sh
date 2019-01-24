#!/bin/bash

bsize=56
rhour=0
rmin=45
dim=95

for ((seed=1122; seed<=2122; seed=seed+1))
do
	sed -e "s;SEED;$seed;g" -e "s;BSIZE;$bsize;g" -e "s;DIM;$dim;g" -e "s;MAXHOUR;$rhour;g" -e "s;MAXMINUTE;$rmin;g" template.sh > batches/"enum$dim-$bsize-$seed".sh
done
