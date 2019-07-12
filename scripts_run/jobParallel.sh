#!/bin/bash
# Set the name of the job
#SBATCH -J netFieldEvo
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks 24
#SBATCH -o ../output/outfile_%j.txt
#SBATCH -e ../output/outfile_%j.txt
#SBATCH --qos=normal
#SBATCH -A ucb-general

# The following commands will be executed when this script is run.
module purge
module load python/2.7.11
module load intel
module load impi
cd ../python
mpirun -n 4 python netFieldEvoSSExors.py 5000 5001 5002 5003


# End of example job shell script
#
