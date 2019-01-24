#!/bin/sh
#SBATCH -J paper-90-6880
#SBATCH -e /home/mb31simi/NewPenum/penum/errs/paper-dim90-seed6880.err.%j
#SBATCH -o /home/mb31simi/NewPenum/penum/logs/paper-dim90-seed6880.out.%j
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -C avx2
#SBATCH --mem-per-cpu=2000
#SBATCH --exclusive
#SBATCH -t 6:32:59

SEED=6880
DIM=90
BIT=900

BETA=48
PREBETA=38
PRUNEPARAM=06

AMAX=2436
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
export OMP_NUM_THREADS=24
#/home/mb31simi/svpgenerator/generate_random --dim $DIM --seed $SEED --bit $BIT --bitexponent 1.0 &> /home/mb31simi/NewPenum/penum/inppap/$FILENAME.txt

for ((dim=1; dim<=15; dim=dim+1))
do
	/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-b3848.log
done