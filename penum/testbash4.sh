#!/bin/bash

module load gcc/8.2 boost fplll/2017.12_avx

export OMP_NUM_THREADS=5
cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "tinput/svpc-e1-d100-0beta56.txt" --delta 0.999 --beta -1 --sheight 10 --amax 2672 >> ScaleTestTo60T.txt
export OMP_NUM_THREADS=1
cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "tinput/svpc-e1-d100-0beta56.txt" --delta 0.999 --beta -1 --sheight 10 --amax 2672 >> ScaleTestTo60T.txt

for ((run=1; run<=3; run=run+1))
do
export OMP_NUM_THREADS=60
cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "tinput/svpc-e1-d100-0beta56.txt" --delta 0.999 --beta -1 --sheight 10 --amax 2672 >> ScaleTestTo60T.txt
export OMP_NUM_THREADS=50
cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "tinput/svpc-e1-d100-0beta56.txt" --delta 0.999 --beta -1 --sheight 10 --amax 2672 >> ScaleTestTo60T.txt
export OMP_NUM_THREADS=40
cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "tinput/svpc-e1-d100-0beta56.txt" --delta 0.999 --beta -1 --sheight 10 --amax 2672 >> ScaleTestTo60T.txt
export OMP_NUM_THREADS=30
cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "tinput/svpc-e1-d100-0beta56.txt" --delta 0.999 --beta -1 --sheight 10 --amax 2672 >> ScaleTestTo60T.txt
export OMP_NUM_THREADS=20
cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "tinput/svpc-e1-d100-0beta56.txt" --delta 0.999 --beta -1 --sheight 10 --amax 2672 >> ScaleTestTo60T.txt
export OMP_NUM_THREADS=10
cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "tinput/svpc-e1-d100-0beta56.txt" --delta 0.999 --beta -1 --sheight 10 --amax 2672 >> ScaleTestTo60T.txt
done

export OMP_NUM_THREADS=5
cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "tinput/svpc-e1-d100-0beta56.txt" --delta 0.999 --beta -1 --sheight 10 --amax 2672 >> ScaleTestTo60T.txt
export OMP_NUM_THREADS=1
cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "tinput/svpc-e1-d100-0beta56.txt" --delta 0.999 --beta -1 --sheight 10 --amax 2672 >> ScaleTestTo60T.txt
