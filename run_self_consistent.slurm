#!/bin/bash
#SBATCH --job-name=self_consistent_SN_fej

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

ini_dir="/home/jseiler/self_consistent_SAGE/ini_files/SN_fej/"
log_dir="/home/jseiler/self_consistent_SAGE/logs/SN_fej/"
run_prefix="SN_fej_0.7_0.0" # Attaches a prefix to the ini and log files to make it easier to discriminate between runs. 

LowSnap=27
HighSnap=97

NUMPROC=32

for ((snap_idx = LowSnap; snap_idx < HighSnap+ 1; ++snap_idx))
{

  echo "Running Snapshot ${snap_idx}"
  SAGE_ini="${ini_dir}/${run_prefix}_SAGE_snap${snap_idx}.ini"
  cifog_ini="${ini_dir}/${run_prefix}_cifog_snap${snap_idx}.ini"

  # Run SAGE
  echo "Running SAGE"
  mpirun -np ${NUMPROC} /home/jseiler/self_consistent_SAGE/sage/sage ${SAGE_ini} > ${log_dir}/${run_prefix}_SAGE_${snap_idx}.log
  exit_code=$? 
  if [ $exit_code -ne 0 ]; then
    echo "SAGE exited with errorcode $exit_code"
    exit $exit_code 
  fi

  # Run cifog with -r -s for 1 snapshot
  if ((snap_idx == LowSnap)) ; then
    cifog_run_flag="-s"
  else
    cifog_run_flag="-r -s"
  fi   

  echo "Running cifog"
  mpirun -np ${NUMPROC} /home/jseiler/grid-model/cifog ${cifog_ini} ${cifog_run_flag} > ${log_dir}/${run_prefix}_cifog_${snap_idx}.log
  exit_code=$? 
  if [ $exit_code -ne 0 ]; then
    echo "cifog exited with errorcode $exit_code"
    exit $exit_code 
  fi
    
  # Generate reionization redshift grid and increment SAGE+cifog ini
  echo "Running Python script"
  python3 /home/jseiler/self_consistent_SAGE/misc_functions.py -f ${SAGE_ini} -c ${cifog_ini} -r 1 -i ${ini_dir} -x ${run_prefix}
  exit_code=$? 
  if [ $exit_code -ne 0 ]; then
    echo "The Python cleanup script exited with errorcode $exit_code"
    exit $exit_code 
  fi

  # Generate reionization modifiers + filter masses
  if ((snap_idx == LowSnap)) ; then
    filter_mass_run_flag=1
  else
    filter_mass_run_flag=0
  fi   

  echo "Generating filter mass" 
  mpirun -np ${NUMPROC} /home/jseiler/self_consistent_SAGE/create_filter_mass ${SAGE_ini}${ini_tag} ${snap_idx} ${filter_mass_run_flag} > ${log_dir}/${run_prefix}_filtermass_${snap_idx}.log 
  exit_code=$? 
  if [ $exit_code -ne 0 ]; then
    echo "create_filter_mass exited with errorcode $exit_code"
    exit $exit_code 
  fi

} 
