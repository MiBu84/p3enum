#!/bin/sh
#SBATCH -J scale-90
#SBATCH -e /home/mb31simi/NewPenum/penum/errs/scale-dim90-seed0.err.%j
#SBATCH -o /home/mb31simi/NewPenum/penum/logs/scale-dim90-seed0.out.%j
#SBATCH -n 1
#SBATCH -c 60
#SBATCH -C avx
#SBATCH --mem-per-cpu=5000
#SBATCH --exclusive
#SBATCH -t 00:29:59

SEED=0
DIM=90
BIT=900

BETA=32
PREBETA=22
PRUNEPARAM=02

AMAX=-1
FILENAME=svpc-e1-d$DIM-$SEED


cd /home/mb31simi/NewPenum/penum/

mkdir -p /home/mb31simi/NewPenum/penum/300

cd /home/mb31simi/NewPenum/penum/300 

env | tee -a $FILENAME.log

hostname | tee -a paper-$FILENAME.log

source /etc/profile.d/modules.sh
module purge
module load gcc/8.2
module load fplll/2017.12_avx

cd /home/mb31simi/NewPenum/penum/

date
export OMP_NUM_THREADS=60
/home/mb31simi/NewPenum/penum/penum --basisfile "input/svpc-e1-d90.txt" --delta 0.999 --prebeta 22 --pruneparam 0.02 --beta 32 --sheight 10 --amax -1 >> paper/scale-e1-$DIM-$SEED-r4.log
export OMP_NUM_THREADS=50
/home/mb31simi/NewPenum/penum/penum --basisfile "input/svpc-e1-d90.txt" --delta 0.999 --prebeta 22 --pruneparam 0.02 --beta 32 --sheight 10 --amax -1 >> paper/scale-e1-$DIM-$SEED-r4.log
export OMP_NUM_THREADS=40
/home/mb31simi/NewPenum/penum/penum --basisfile "input/svpc-e1-d90.txt" --delta 0.999 --prebeta 22 --pruneparam 0.02 --beta 32 --sheight 10 --amax -1 >> paper/scale-e1-$DIM-$SEED-r4.log
export OMP_NUM_THREADS=30
/home/mb31simi/NewPenum/penum/penum --basisfile "input/svpc-e1-d90.txt" --delta 0.999 --prebeta 22 --pruneparam 0.02 --beta 32 --sheight 10 --amax -1 >> paper/scale-e1-$DIM-$SEED-r4.log
export OMP_NUM_THREADS=20
/home/mb31simi/NewPenum/penum/penum --basisfile "input/svpc-e1-d90.txt" --delta 0.999 --prebeta 22 --pruneparam 0.02 --beta 32 --sheight 10 --amax -1 >> paper/scale-e1-$DIM-$SEED-r4.log
export OMP_NUM_THREADS=10
/home/mb31simi/NewPenum/penum/penum --basisfile "input/svpc-e1-d90.txt" --delta 0.999 --prebeta 22 --pruneparam 0.02 --beta 32 --sheight 10 --amax -1 >> paper/scale-e1-$DIM-$SEED-r4.log
export OMP_NUM_THREADS=5
/home/mb31simi/NewPenum/penum/penum --basisfile "input/svpc-e1-d90.txt" --delta 0.999 --prebeta 22 --pruneparam 0.02 --beta 32 --sheight 10 --amax -1 >> paper/scale-e1-$DIM-$SEED-r4.log
export OMP_NUM_THREADS=1
/home/mb31simi/NewPenum/penum/penum --basisfile "input/svpc-e1-d90.txt" --delta 0.999 --prebeta 22 --pruneparam 0.02 --beta 32 --sheight 10 --amax -1 >> paper/scale-e1-$DIM-$SEED-r4.log



#/home/mb31simi/svpgenerator/generate_random --dim $DIM --seed $SEED --bit $BIT --bitexponent 1.0 &> /home/mb31simi/NewPenum/penum/inppap/$FILENAME.txt

#for ((runi=0; runi<15; runi=runi+1))
#do
#	/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-b3848.log
#done
date
