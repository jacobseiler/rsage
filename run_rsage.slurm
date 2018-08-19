#!/bin/bash
#SBATCH --job-name=shifted_SFR

#SBATCH --ntasks=32
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --nodes=1

module purge
module load scipy/1.0.0-python-3.6.4
module load gsl/2.4
module load astropy/2.0.3-python-3.6.4
module load numpy/1.14.1-python-3.6.4
module load openmpi/3.0.0
module load gcc/6.4.0

NUMPROC=32

SAGE_ini="/home/jseiler/tmp/rsage/ini_files/kali_SAGE.ini"
cifog_ini="/home/jseiler/tmp/rsage/ini_files/kali_cifog.ini"

path_to_exe="/home/jseiler/rsage/bin/rsage"

run_prefix="test_RSAGE" # Attaches a prefix to the ini and log files to make it easier to discriminate between runs. 
path_to_log="/home/jseiler/rsage/logs/${run_prefix}.log"

echo "Running ${path_to_exe} on ${NUMPROC} processors with SAGE ini file ${SAGE_ini} and cifog ini file ${cifog_ini}"

mpirun -np ${NUMPROC} ${path_to_exe} ${SAGE_ini} ${cifog_ini} > ${path_to_log} 
