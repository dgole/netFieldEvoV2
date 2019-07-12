#!/bin/bash

#module purge
#module load gcc
#module load python/3.6.1

# runId, mStar, bStar, mdot0, tCycle, nCycles

#python3 netFieldEvoSSExors.py ${1}00 $2 -1.e3 1.e-8 2.e1 100
#python3 workspaceExors.py     ${1}01

#python3 netFieldEvoSSExors.py ${1}10 $2 -1.e3 1.e-8 2.e1 100
#python3 workspaceExors.py     ${1}11

#python3 netFieldEvoSSExors.py ${1}20 $2 -1.e3 1.e-9 2.e1 100
#python3 workspaceExors.py     ${1}21

#python3 netFieldEvoSSExors.py ${1}30 $2 -1.e3 1.e-10 2.e1 100
#python3 workspaceExors.py     ${1}31

python3 netFieldEvoSSExors.py ${1}00 $2 -1.e3 1.e-7 1.0 1

python3 netFieldEvoSSExors.py ${1}10 $2 -1.e3 1.e-8 1.0 1

python3 netFieldEvoSSExors.py ${1}20 $2 -1.e3 1.e-9 1.0 1

python3 netFieldEvoSSExors.py ${1}30 $2 -1.e3 1.e-10 1.0 1


#
