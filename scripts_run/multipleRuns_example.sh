#!/bin/bash

#module purge
#module load gcc
#module load python/3.6.1

#                             runId  mStar  bInit   mdot0  FCF  tCycle nCycles

# fiducial outbursts
python3 ../python/netFieldEvoSSExors.py 10   0.8   -1.e-3  1.e-8  1.0  20   20  # no

#python3 netFieldEvoSSExors.py 10001   0.8   -1.e-3  1.e-9  1.0  20   20

#python3 netFieldEvoSSExors.py 10002   0.8   -1.e-3  1.e-10  1.0  20   20
