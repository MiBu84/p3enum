#!/bin/sh
#SBATCH -J paper-70-98937
#SBATCH -e /home/mb31simi/NewPenum/penum/errs/paper-dim70-seed98937.err.%j
#SBATCH -o /home/mb31simi/NewPenum/penum/logs/paper-dim70-seed98937.out.%j
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -C avx
#SBATCH --mem-per-cpu=1800
#SBATCH --exclusive
#SBATCH -t 00:22:59

SEED=98937
DIM=70
BIT=700

BETA=32
PREBETA=24
PRUNEPARAM=2

AMAX=2145
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
export OMP_NUM_THREADS=16
/home/mb31simi/svpgenerator/generate_random --dim $DIM --seed $SEED --bit $BIT --bitexponent 1.0 &> /home/mb31simi/NewPenum/penum/inppap/$FILENAME.txt

/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-v2.log
/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-v2.log
/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-v2.log
/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-v2.log
/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-v2.log
/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-v2.log
/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-v2.log
/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-v2.log
/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-v2.log
/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-v2.log
/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-v2.log
/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-v2.log
/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-v2.log
/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-v2.log
/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA --sheight 10 --amax $AMAX >> paper/svpc-e1-$DIM-$SEED-v2.log

%/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA--sheight 10 --amax -1 --itenum --aend $AMAX >> paper/svpc-e1-$DIM-$SEED.log
%/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA--sheight 10 --amax -1 --itenum --aend $AMAX >> paper/svpc-e1-$DIM-$SEED.log
%/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA--sheight 10 --amax -1 --itenum --aend $AMAX >> paper/svpc-e1-$DIM-$SEED.log
%/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA--sheight 10 --amax -1 --itenum --aend $AMAX >> paper/svpc-e1-$DIM-$SEED.log
%/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA--sheight 10 --amax -1 --itenum --aend $AMAX >> paper/svpc-e1-$DIM-$SEED.log
%/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA--sheight 10 --amax -1 --itenum --aend $AMAX >> paper/svpc-e1-$DIM-$SEED.log
%/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA--sheight 10 --amax -1 --itenum --aend $AMAX >> paper/svpc-e1-$DIM-$SEED.log
%/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA--sheight 10 --amax -1 --itenum --aend $AMAX >> paper/svpc-e1-$DIM-$SEED.log
%/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA--sheight 10 --amax -1 --itenum --aend $AMAX >> paper/svpc-e1-$DIM-$SEED.log
%/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA--sheight 10 --amax -1 --itenum --aend $AMAX >> paper/svpc-e1-$DIM-$SEED.log
%/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA--sheight 10 --amax -1 --itenum --aend $AMAX >> paper/svpc-e1-$DIM-$SEED.log
%/home/mb31simi/NewPenum/penum/penum --basisfile "./inppap/svpc-e1-d$DIM-$SEED.txt" --delta 0.9999 --prebeta $PREBETA --pruneparam 0.$PRUNEPARAM --beta $BETA--sheight 10 --amax -1 --itenum --aend $AMAX >> paper/svpc-e1-$DIM-$SEED.log
date