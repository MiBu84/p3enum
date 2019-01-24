#!/bin/sh
#SBATCH -J enum-500-94
#SBATCH -e /home/mb31simi/NewPenum/penum/errs/enum-bs56-dim94-seed500.err.%j
#SBATCH -o /home/mb31simi/NewPenum/penum/logs/enum-bs56-dim94-seed500.out.%j
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -C avx2
#SBATCH --mem-per-cpu=1800
#SBATCH --exclusive
#SBATCH -t 1:59:59

cd /home/mb31simi/NewPenum/penum/

mkdir -p /home/mb31simi/NewPenum/penum/100

cd /home/mb31simi/NewPenum/penum/100 

env | tee -a enum-bs56-dim94-seed500.log

hostname | tee -a enum-bs56-dim94-seed500.log

source /etc/profile.d/modules.sh
module purge
module load gcc/8.2
module load mpfr/3.1.6
module load gmp/6.1.2
module load boost
module load fplll

cd /home/mb31simi/NewPenum/penum/

date
export OMP_NUM_THREADS=24
/home/mb31simi/svpgenerator/generate_random --dim 94 --seed 500 --bit 940 --bitexponent 1.0 &> /home/mb31simi/NewPenum/penum/tinput/svpc-e1-d94-500.txt
/home/mb31simi/NewPenum/penum/penum --basisfile "./tinput/svpc-e1-d94-500.txt" --delta 0.9999 --beta 56 --sheight 10 --amax 2450 &> log94/enum-bs56-dim94-seed500.log
date
