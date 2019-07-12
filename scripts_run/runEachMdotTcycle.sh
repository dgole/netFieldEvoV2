#!/bin/bash

#module purge
#module load gcc
#module load python/3.6.1

# runId, mStar, bStar, mdot0, tCycle, nCycles

################################################################################
#python3 netFieldEvoSSExors.py ${1}00 $2 -1.e3 1.e-7 6 40
#python3 workspaceExors.py     ${1}00

#python3 netFieldEvoSSExors.py ${1}01 $2 -1.e3 1.e-8 6 40
#python3 workspaceExors.py     ${1}01

#python3 netFieldEvoSSExors.py ${1}02 $2 -1.e3 1.e-9 6 40
#python3 workspaceExors.py     ${1}02

#python3 netFieldEvoSSExors.py ${1}03 $2 -1.e3 1.e-10 6 40
#python3 workspaceExors.py     ${1}03
################################################################################
#python3 netFieldEvoSSExors.py ${1}10 $2 -1.e3 1.e-7 20 20
#python3 workspaceExors.py     ${1}10

#python3 netFieldEvoSSExors.py ${1}11 $2 -1.e3 1.e-8 20 20
#python3 workspaceExors.py     ${1}11

python3 netFieldEvoSSExors.py ${1}12 $2 -1.e3 1.e-9 20 20
python3 workspaceExors.py     ${1}12

python3 netFieldEvoSSExors.py ${1}13 $2 -1.e3 1.e-10 20 20
python3 workspaceExors.py     ${1}13
################################################################################
#python3 netFieldEvoSSExors.py ${1}20 $2 -1.e3 1.e-7 60 10
#python3 workspaceExors.py     ${1}20

#python3 netFieldEvoSSExors.py ${1}21 $2 -1.e3 1.e-8 60 10
#python3 workspaceExors.py     ${1}21

#python3 netFieldEvoSSExors.py ${1}22 $2 -1.e3 1.e-9 60 10
#python3 workspaceExors.py     ${1}22

#python3 netFieldEvoSSExors.py ${1}23 $2 -1.e3 1.e-10 60 10
#python3 workspaceExors.py     ${1}23
################################################################################














#
