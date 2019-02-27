#!/bin/bash

module purge
module load gcc
module load python/3.6.1

# means still needs to be run
## means done


# rDz=0.1, mdot=1.e-6, binit=-1.e-2
python3 netFieldEvoSSExors.py 2000
python3 netFieldEvoSSExors.py 2010

# rDz=0.2, mdot=1.e-6, binit=-1.e-2
python3 netFieldEvoSSExors.py 2001
python3 netFieldEvoSSExors.py 2011

# rDz=0.3, mdot=1.e-6, binit=-1.e-2
python3 netFieldEvoSSExors.py 2002
python3 netFieldEvoSSExors.py 2012

# rDz=0.4, mdot=1.e-6, binit=-1.e-2
python3 netFieldEvoSSExors.py 2003
python3 netFieldEvoSSExors.py 2013

# rDz=0.5, mdot=1.e-6, binit=-1.e-2
python3 netFieldEvoSSExors.py 2004
python3 netFieldEvoSSExors.py 2014

# rDz=0.6, mdot=1.e-6, binit=-1.e-2
python3 netFieldEvoSSExors.py 2005
python3 netFieldEvoSSExors.py 2015

# rDz=0.7, mdot=1.e-6, binit=-1.e-2
python3 netFieldEvoSSExors.py 2006
python3 netFieldEvoSSExors.py 2016

# rDz=0.8, mdot=1.e-6, binit=-1.e-2
python3 netFieldEvoSSExors.py 2007
python3 netFieldEvoSSExors.py 2017

# rDz=0.9, mdot=1.e-6, binit=-1.e-2
python3 netFieldEvoSSExors.py 2008
python3 netFieldEvoSSExors.py 2018

# rDz=1.0, mdot=1.e-6, binit=-1.e-2
python3 netFieldEvoSSExors.py 2009
python3 netFieldEvoSSExors.py 2019


# 201? is same but with tCycle=30




















