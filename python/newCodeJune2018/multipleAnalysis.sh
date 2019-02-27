#!/bin/bash

module purge
module load gcc
module load python/3.6.1

# means still needs to be run
## means done


python3 workspaceExors.py ${1}00 
python3 workspaceExors.py 1001 
python3 workspaceExors.py 1001 
python3 workspaceExors.py 1001 
python3 workspaceExors.py 1001 

