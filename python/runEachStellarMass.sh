#!/bin/bash

module purge
module load gcc
module load python/3.6.1

# runIdBase, bStar, mdot0, tCycle, nCycles

# runId, mStar, bStar, mdot0, tCycle, nCycles

python3 netFieldEvoSSExors.py ${1}0 1.0 $2 $3 $4 $5
python3 workspaceExors.py     ${1}0

#python3 netFieldEvoSSExors.py ${1}1 0.9 $2 $3 $4 $5
#python3 workspaceExors.py     ${1}1

#python3 netFieldEvoSSExors.py ${1}2 0.8 $2 $3 $4 $5
#python3 workspaceExors.py     ${1}2

python3 netFieldEvoSSExors.py ${1}3 0.7 $2 $3 $4 $5
python3 workspaceExors.py     ${1}3

#python3 netFieldEvoSSExors.py ${1}4 0.6 $2 $3 $4 $5
#python3 workspaceExors.py     ${1}4

python3 netFieldEvoSSExors.py ${1}5 0.5 $2 $3 $4 $5
python3 workspaceExors.py     ${1}5

#python3 netFieldEvoSSExors.py ${1}6 0.4 $2 $3 $4 $5
#python3 workspaceExors.py     ${1}6

python3 netFieldEvoSSExors.py ${1}7 0.3 $2 $3 $4 $5
python3 workspaceExors.py     ${1}7

#python3 netFieldEvoSSExors.py ${1}8 0.2 $2 $3 $4 $5
#python3 workspaceExors.py     ${1}8

python3 netFieldEvoSSExors.py ${1}9 0.1 $2 $3 $4 $5
python3 workspaceExors.py     ${1}9

























