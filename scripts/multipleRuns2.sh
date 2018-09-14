#!/bin/bash

module purge
module load gcc
module load python/3.6.1

python3 ../python/netFieldEvoSSExors.py 20500
python3 ../python/netFieldEvoSSExors.py 20501
python3 ../python/netFieldEvoSSExors.py 20502

























