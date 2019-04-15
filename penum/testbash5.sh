#!/bin/bash

export OMP_NUM_THREADS=32


    cd /home/mb31simi/NewPenum/penum
    cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "input/svpc-e1-d96.txt" --delta 0.9999 --beta 60 --sheight 10 --amax -1 --itenum --aend 2525 >> TestPerf.txt
    cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "input/svpc-e1-d96.txt" --delta 0.9999 --beta 58 --sheight 10 --amax -1 --itenum --aend 2525 >> TestPerf.txt
    cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "input/svpc-e1-d96.txt" --delta 0.9999 --beta 56 --sheight 10 --amax -1 --itenum --aend 2525 >> TestPerf.txt
    cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "input/svpc-e1-d96.txt" --delta 0.9999 --beta 54 --sheight 10 --amax -1 --itenum --aend 2525 >> TestPerf.txt
    cd /home/mb31simi/NewPenum/penum && ./penum --basisfile "input/svpc-e1-d96.txt" --delta 0.9999 --beta 62 --sheight 10 --amax -1 --itenum --aend 2525 >> TestPerf.txt
    
