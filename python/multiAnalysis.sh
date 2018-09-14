#!/bin/bash

module purge
module load gcc
module load python/3.6.1

#python3 workspaceExors.py 1.0    10200 10210 10220 10230 10240
#python3 workspaceExors.py 3.2    10201 10211 10221 10231 10241
#python3 workspaceExors.py 10.0   10202 10212 10222 10232 10242
#python3 workspaceExors.py 32.0   10203 10213 10223 10233 10243
#python3 workspaceExors.py 100.0  10204 10214 10224 10234 10244

python3 workspaceExors.py 1.0    10240 
python3 workspaceExors.py 3.2    10241 
python3 workspaceExors.py 10.0   10242 
python3 workspaceExors.py 32.0   10243 
python3 workspaceExors.py 100.0  10244 

python3 workspaceExors.py 1.0    20000 
python3 workspaceExors.py 3.2    20001 
python3 workspaceExors.py 10.0   20002 
python3 workspaceExors.py 32.0   20003 
python3 workspaceExors.py 100.0  20004 

python3 workspaceExors.py 1.0    20010 
python3 workspaceExors.py 3.2    20011 
python3 workspaceExors.py 10.0   20012 
python3 workspaceExors.py 32.0   20013 
python3 workspaceExors.py 100.0  20014

python3 workspaceExors.py 10.0   20022 
python3 workspaceExors.py 10.0   20032 
