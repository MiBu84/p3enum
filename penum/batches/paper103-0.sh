#!/bin/sh
#SBATCH -J paper-103-0
#SBATCH -e /home/mb31simi/NewPenum/penum/errs/paper-dim103-seed0.err.%j
#SBATCH -o /home/mb31simi/NewPenum/penum/logs/paper-dim103-seed0.out.%j
#SBATCH -n 1
#SBATCH -c 60
#SBATCH -C avx
#SBATCH --mem-per-cpu=5000
#SBATCH --exclusive
#SBATCH -t 02:09:59

SEED=0
DIM=103
BIT=1030

BETA=52
PREBETA=36
PRUNEPARAM=0

AMAX=2606
FILENAME=svpc-e1-d$DIM-$SEED


cd /home/mb31simi/NewPenum/penum/

mkdir -p /home/mb31simi/NewPenum/penum/200

cd /home/mb31simi/NewPenum/penum/200 

env | tee -a $FILENAME.log

hostname | tee -a paper-$FILENAME.log

source /etc/profile.d/modules.sh
module purge
module load gcc/8.2
module load fplll/2017.12_avx

cd /home/mb31simi/NewPenum/penum/

date
#/home/mb31simi/svpgenerator/generate_random --dim $DIM --seed $SEED --bit $BIT --bitexponent 1.0 &> /home/mb31simi/NewPenum/penum/inppap/$FILENAME.txt
/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED.log
date
