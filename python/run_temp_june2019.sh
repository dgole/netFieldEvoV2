#!/bin/bash

#module purge
#module load gcc
#module load python/3.6.1

#                             runId  mStar  bInit   mdot0  FCF  tCycle nCycles

# fiducial outbursts
#python3 netFieldEvoSSExors.py 10000   0.8   -1.e-3  1.e-8  1.0  20   20  # no
#python3 workspaceExors.py     10000

#python3 netFieldEvoSSExors.py 10001   0.8   -1.e-3  1.e-9  1.0  20   20 # big
#python3 workspaceExors.py     10001

#python3 netFieldEvoSSExors.py 10002   0.8   -1.e-3  1.e-10  1.0  20   20 # big/redundant
#python3 workspaceExors.py     10002

#python3 netFieldEvoSSExors.py 10003   0.8   -1.e-3  1.e-8  1.0  20   10 # no
#python3 workspaceExors.py     10003

#python3 netFieldEvoSSExors.py 10004   0.8   -1.e-3  1.e-9  1.0  40   10 # big
#python3 workspaceExors.py     10004

#python3 netFieldEvoSSExors.py 10005   0.8   -1.e-3  1.e-10  1.0  40   10 # big/redundant
#python3 workspaceExors.py     10005

#python3 netFieldEvoSSExors.py 10006   0.8   -1.e-3  1.e-8  1.0  6   20 # no
#python3 workspaceExors.py     10006

#python3 netFieldEvoSSExors.py 10007   0.8   -1.e-3  1.e-9  1.0  6   20 # ok, kinda weird
#python3 workspaceExors.py     10007



#python3 netFieldEvoSSExors.py 20000   0.5   -1.e-3  1.e-8  1.0  20   20  # no
#python3 workspaceExors.py     20000

#python3 netFieldEvoSSExors.py 20001   0.5   -1.e-3  1.e-9  1.0  20   20 # good marginal
#python3 workspaceExors.py     20001

#python3 netFieldEvoSSExors.py 20002   0.5   -1.e-3  1.e-10  1.0  20   20 # good
#python3 workspaceExors.py     20002

#python3 netFieldEvoSSExors.py 20003   0.5   -1.e-3  1.e-8  1.0  20   10 # no
#python3 workspaceExors.py     20003

#python3 netFieldEvoSSExors.py 20004   0.5   -1.e-3  1.e-9  1.0  40   10 # good, different
#python3 workspaceExors.py     20004

#python3 netFieldEvoSSExors.py 20005   0.5   -1.e-3  1.e-10  1.0  40   10 # big, little different
#python3 workspaceExors.py     20005

#python3 netFieldEvoSSExors.py 20006   0.5   -1.e-3  1.e-9  1.0  6   20
#python3 workspaceExors.py     20006

#python3 netFieldEvoSSExors.py 20007   0.5   -1.e-3  1.e-10  1.0  6   20
#python3 workspaceExors.py     20007



#python3 netFieldEvoSSExors.py 30000   0.3   -1.e-3  1.e-8  1.0  20   20 # no
#python3 workspaceExors.py     30000

#python3 netFieldEvoSSExors.py 30001   0.3   -1.e-3  1.e-9  1.0  20   20 # marginal and nice
#python3 workspaceExors.py     30001

#python3 netFieldEvoSSExors.py 30002   0.3   -1.e-3  1.e-10  1.0  20   20 # big and nice
#python3 workspaceExors.py     30002

#python3 netFieldEvoSSExors.py 30003   0.3   -1.e-3  1.e-8  1.0  20   10 # no
#python3 workspaceExors.py     30003

#python3 netFieldEvoSSExors.py 30004   0.3   -1.e-3  1.e-9  1.0  40   10 # good, different
#python3 workspaceExors.py     30004

#python3 netFieldEvoSSExors.py 30005   0.3   -1.e-3  1.e-10  1.0  40   10 # good, different
#python3 workspaceExors.py     30005

#python3 netFieldEvoSSExors.py 30006   0.3   -1.e-3  1.e-9  1.0  6   20
#python3 workspaceExors.py     30006

#python3 netFieldEvoSSExors.py 30007   0.3   -1.e-3  1.e-10  1.0  6   20
#python3 workspaceExors.py     30007




#python3 netFieldEvoSSExors.py 40000   0.1   -1.e-3  1.e-8  1.0  20   20 # no
#python3 workspaceExors.py     40000

#python3 netFieldEvoSSExors.py 40001   0.1   -1.e-3  1.e-9  1.0  20   20 # good, different
#python3 workspaceExors.py     40001

#python3 netFieldEvoSSExors.py 40002   0.1   -1.e-3  1.e-10  1.0  20   20
#python3 workspaceExors.py     40002

#python3 netFieldEvoSSExors.py 40003   0.1   -1.e-3  1.e-8  1.0  20   10
#python3 workspaceExors.py     40003

#python3 netFieldEvoSSExors.py 40004   0.1   -1.e-3  1.e-9  1.0  40   10
#python3 workspaceExors.py     40004

#python3 netFieldEvoSSExors.py 40005   0.1   -1.e-3  1.e-10  1.0  40   10
#python3 workspaceExors.py     40005

#python3 netFieldEvoSSExors.py 40006   0.1   -1.e-3  1.e-9  1.0  6   20
#python3 workspaceExors.py     40006

#python3 netFieldEvoSSExors.py 40007   0.1   -1.e-3  1.e-10  1.0  6   20
#python3 workspaceExors.py     40007



#python3 netFieldEvoSSExors.py 50000   1.0   -1.e-3  1.e-8  1.0  20   20
#python3 workspaceExors.py     50000

#python3 netFieldEvoSSExors.py 50001   1.0   -1.e-3  1.e-9  1.0  20   20
#python3 workspaceExors.py     50001

#python3 netFieldEvoSSExors.py 50002   1.0   -1.e-3  1.e-10  1.0  20   20
#python3 workspaceExors.py     50002

#python3 netFieldEvoSSExors.py 50003   1.0   -1.e-3  1.e-8  1.0  20   10
#python3 workspaceExors.py     50003

#python3 netFieldEvoSSExors.py 50004   1.0   -1.e-3  1.e-9  1.0  40   10
#python3 workspaceExors.py     50004

#python3 netFieldEvoSSExors.py 50005   1.0   -1.e-3  1.e-10  1.0  40   10
#python3 workspaceExors.py     50005

python3 netFieldEvoSSExors.py 50006   1.0   -1.e-3  1.e-9  1.0  6   20
python3 workspaceExors.py     50006

#python3 netFieldEvoSSExors.py 50007   1.0   -1.e-3  1.e-10  1.0  6   20
#python3 workspaceExors.py     50007
