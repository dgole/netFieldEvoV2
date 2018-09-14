#!/bin/bash
# Set the name of the job
#SBATCH -J netFieldEvo
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks 24
#SBATCH -o ../output/outfile_6000.txt
#SBATCH -e ../output/errfile_6000.txt
#SBATCH --qos=normal
#SBATCH -A ucb-general

# The following commands will be executed when this script is run.
module purge
module load python/2.7.11
cd ../python
python netFieldEvoSSExors.py 6000


# End of example job shell script
#
