#!/bin/bash
#SBATCH --job-name=corrfunc

#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=20gb
#SBATCH --nodes=1

module load scipy/1.0.0-python-3.6.4
module load gsl/2.4
module load astropy/2.0.3-python-3.6.4
module load numpy/1.14.1-python-3.6.4
module load openmpi/3.0.0
module load gcc/6.4.0
module load matplotlib/2.2.2-python-3.6.4
module load mpi4py/3.0.0-python-3.6.4
module load h5py/2.7.1-python-3.6.4-serial
module load ipython/5.5.0-python-3.6.4
module load git/2.16.0

PROGNAME="/home/jseiler/self_consistent_SAGE/output/corrfunc_checking.py"

python3 $PROGNAME
