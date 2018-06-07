#!/bin/bash
#SBATCH --job-name=myresults

#SBATCH --ntasks=8
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=19900
#SBATCH --nodes=1

module load scipy/1.0.0-python-3.6.4
module load gsl/2.4
module load astropy/2.0.3-python-3.6.4
module load numpy/1.14.1-python-3.6.4
module load openmpi/3.0.0
module load gcc/6.4.0
module load matplotlib/2.2.2-python-3.6.4
module load mpi4py/3.0.0-python-3.6.4

PROGNAME="/home/jseiler/self_consistent_SAGE/output/myresults.py"
NUMPROC=8

mpirun -np $NUMPROC python3 $PROGNAME
