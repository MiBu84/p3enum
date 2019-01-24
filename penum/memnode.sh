#!/bin/sh
#SBATCH -J enum-500-94
#SBATCH -e /home/mb31simi/NewPenum/penum/errs/bigmem32.err.%j
#SBATCH -o /home/mb31simi/NewPenum/penum/logs/bigmem32.out.%j
#SBATCH -n 1
#SBATCH -c 60
#SBATCH -C avx
#SBATCH --mem-per-cpu=10000
#SBATCH --exclusive
#SBATCH -t 00:02:59

cd /home/mb31simi/NewPenum/penum/

mkdir -p /home/mb31simi/NewPenum/penum/100

cd /home/mb31simi/NewPenum/penum/100 

env | tee -a enum-bs56-dim94-seed500.log

hostname | tee -a enum-bs56-dim94-seed500.log

source /etc/profile.d/modules.sh
module purge
module load gcc/8.2 boost fplll/2017.12_avx
cd /home/mb31simi/NewPenum/penum/

date
export OMP_NUM_THREADS=32
/home/mb31simi/NewPenum/penum/penum --basisfile "tinput/svpc-e1-d100-0beta56.txt" --delta 0.999 --beta -1 --sheight 10 --amax 2672 >> stinkipinki.txt
date
