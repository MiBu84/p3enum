#!/bin/sh
#SBATCH -J enum-SEED-DIM
#SBATCH -e /home/mb31simi/NewPenum/penum/errs/enum-bsBSIZE-dimDIM-seedSEED.err.%j
#SBATCH -o /home/mb31simi/NewPenum/penum/logs/enum-bsBSIZE-dimDIM-seedSEED.out.%j
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -C avx
#SBATCH --mem-per-cpu=1200
#SBATCH --exclusive
#SBATCH -t MAXHOUR:MAXMINUTE:59

cd /home/mb31simi/NewPenum/penum/

mkdir -p /home/mb31simi/NewPenum/penum/500

cd /home/mb31simi/NewPenum/penum/500 

env | tee -a enum-bsBSIZE-dimDIM-seedSEED.log

hostname | tee -a enum-bsBSIZE-dimDIM-seedSEED.log

source /etc/profile.d/modules.sh
module purge
module load gcc/8.2
module load boost/1.66.0
module load fplll/2017.12_avx

cd /home/mb31simi/NewPenum/penum/

date
export OMP_NUM_THREADS=16
/home/mb31simi/svpgenerator/generate_random --dim DIM --seed SEED --bit 880 --bitexponent 1.0 &> /home/mb31simi/NewPenum/penum/tinput/svpc-e1-dDIM-SEED.txt
/home/mb31simi/NewPenum/penum/penum --basisfile "./tinput/svpc-e1-dDIM-SEED.txt" --delta 0.9999 --prebeta 34 --pruneparam 0.1 --beta BSIZE --sheight 10 --amax 2194 &> log88/enum-bsBSIZE-dimDIM-seedSEED.log
date
