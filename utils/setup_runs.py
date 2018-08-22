"""
Creates directories and ``.ini`` files for ``RSAGE`` runs. By specifying lists of
variables in ``SAGE_fields_update`` and ``cifog_fields_update``, you are able
to create a unique combination of variables to update for each run. 

The script also gives the option of making ``slurm`` files with the paths
correctly specified and also submit them (**BE CAREFUL WHEN SUBMITTING**).

All ``.ini`` and ``slurm`` files are created using template files.  Examples
files are included in the base repo.
"""

#!/usr/bin/env python
from __future__ import print_function

import numpy as np
import sys
import os
from shutil import copyfile
import subprocess

sys.path.append('./output/')
import ReadScripts
import AllVars

script_dir = os.path.dirname(os.path.realpath(__file__))


def check_input_parameters(SAGE_fields_update, cifog_fields_update,
                           run_directories, base_SAGE_ini, base_cifog_ini,
                           base_slurm_file):
    """
    Checks that the update fields all have the same number of inputs.  Also
    checks that the template ``.ini`` and ``.slurm`` files exist. 
 
    Parameters
    ----------

    SAGE_fields_update, cifog_fields_update : Dictionaries
        Fields that will be updated and their new value.

    run_directories : List of strings, length equal to number of runs 
        Path to the base ``RSAGE`` directory for each run where all the model
        output will be placed.

    base_SAGE_ini, base_cifog_ini : Strings
        Paths to the template SAGE and cifog ``.ini`` files.

    base_slurm_file : String
        Path to the template ``.slurm`` file.

    Returns
    ----------
    
    None.

    Errors 
    ----------
        ValueError:
            Raised if any of the fields in ``SAGE_fields_update`` or
            ``cifog_fields_update`` have different lengths.

            Raised if any of the fields in ``SAGE_fields_update`` or
            ``cifog_fields_update`` have different lengths to that of
            ``run_directories``.

            Raised in any of the template files do not exist.    
    """

    # Compare every field within ``SAGE_fields_update`` and
    # ``cifog_fields_update`` to ensure they have the same lengths.
    for count, my_dict in enumerate([SAGE_fields_update, cifog_fields_update]):
        for key_1 in SAGE_fields_update.keys(): 
            for key_2 in SAGE_fields_update.keys():
 
                if len(SAGE_fields_update[key_1]) != len(SAGE_fields_update[key_2]):
                    if count == 0:
                        which_dict = "SAGE"
                    else:
                        which_dict = "cifog"
 
                    print("For the {0} update dictionary, Key {1} did not have "
                          "the same length as key {2}".format(which_dict,
                                                              key_1, 
                                                              key_2))
                    raise ValueError

            if len(SAGE_fields_update[key_1]) != len(run_directories):
                print("For the {0} update dictionary, Key {1} did not have "
                      "the same length as the number of run directories" \
                      .format(which_dict, key_1))
                raise ValueError

    # Check all the template files exist.
    for my_file in [base_SAGE_ini, base_cifog_ini, base_slurm_file]: 
        if not os.path.isfile(my_file):

            print("File {0} does not exist.".format(my_file))
            raise ValueError


def create_directories(run_directory):
    """
    Creates the directories to house all the ``RSAGE`` outputs.   
 
    Parameters
    ----------

    run_directory: String 
        Path to the base ``RSAGE`` directory. 

    Returns
    ----------
    
    None.
    """

    # First create the base directory.
    base_dir = run_directory 
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
        print("Created directory {0}".format(base_dir))

    # Directory that will contain all the SAGE output.
    gal_dir = "{0}/galaxies".format(base_dir)
    if not os.path.exists(gal_dir):
        os.makedirs(gal_dir)
        print("Created directory {0}".format(gal_dir))

    # Directory that will contain all the grids. 
    grids_dir = "{0}/grids".format(base_dir)
    if not os.path.exists(grids_dir):
        os.makedirs(grids_dir)
        print("Created directory {0}".format(grids_dir))
       
        # If the grids directory didn't exist, there's no way these will.

        dirs = ["grids/nion", "grids/cifog",
                "grids/cifog/reionization_modifiers",
                "ini_files", "slurm_files", "log_files"]
        for directory in dirs:
            dir_ = "{0}/{1}".format(base_dir, directory)
            os.makedirs(dir_)
            print("Created directory {0}".format(dir_))

   
def update_ini_files(base_SAGE_ini, base_cifog_ini,
                     SAGE_fields_update, cifog_fields_update,
                     run_directory):
    """
    Using template ini files for ``SAGE`` and ``cifog``, creates new ones with 
    the directory paths and field names updated. 

    Parameters
    ----------

    base_SAGE_ini, base_cifog_ini : Strings
        Paths to the template SAGE and cifog ini files.

    SAGE_fields_update, cifog_fields_update : Dictionaries
        Fields that will be updated and their new value.

    run_directory : String
        Path to the base ``RSAGE`` directory.

    Returns
    ----------

    SAGE_fname, cifog_fname : Strings
        Names of the newly created ``SAGE`` and ``cifog`` ini files.
    """

    SAGE_params = ReadScripts.read_SAGE_ini(base_SAGE_ini)
    cifog_params, cifog_headers = ReadScripts.read_cifog_ini(base_cifog_ini)

    # These are paths and don't depend on `FileNameGalaxies`. 
    SAGE_params["OutputDir"] = "{0}/galaxies".format(run_directory)
    SAGE_params["GridOutputDir"] = "{0}/grids/nion".format(run_directory)
    SAGE_params["PhotoionDir"] = "{0}/grids/cifog".format(run_directory)

    # Now go through the parameters and update them.
    for name in SAGE_fields_update:
        SAGE_params[name] = SAGE_fields_update[name] 

    for name in cifog_fields_update:
        cifog_params[name] = cifog_fields_update[name] 

    # The unique identifier amongst each run will be `FileNameGalaxies`. 
    prefix_tag = SAGE_params["FileNameGalaxies"]

    # Now update all the fields for this specific run.
    SAGE_params["PhotoionName"] = "{0}_photHI".format(prefix_tag)
    SAGE_params["ReionRedshiftName"] = "{0}_reionization_redshift" \
                                       .format(prefix_tag)

    # The name for the ionizing photon file depends upon the escape fraction
    # prescription chosen.
    nion_fname = get_nion_fname(SAGE_params)

    cifog_params["inputNionFile"] = "{0}/grids/nion/{1}" \
                                    .format(run_directory, nion_fname)
    cifog_params["output_XHII_file"] = "{0}/grids/cifog/{1}_XHII" \
                                       .format(run_directory,
                                               prefix_tag)
    cifog_params["output_photHI_file"] = "{0}/grids/cifog/{1}_photHI" \
                                         .format(run_directory,
                                                 prefix_tag)
    cifog_params["output_restart_file"] = "{0}/grids/cifog/{1}_restart" \
                                          .format(run_directory,
                                                  prefix_tag)

    # Write out the new ini files, using `FileNameGalaxies` as the tag.
    SAGE_fname = "{0}/ini_files/{1}_SAGE.ini".format(run_directory,
                                                     prefix_tag) 

    cifog_fname = "{0}/ini_files/{1}_cifog.ini".format(run_directory,
                                                       prefix_tag) 

    with open (SAGE_fname, "w+") as f:
        for name in SAGE_params.keys():
            string = "{0} {1}\n".format(name, SAGE_params[name])
            f.write(string)

    with open (cifog_fname, "w+") as f:
        for name in cifog_params.keys():

            if name in cifog_headers.keys():
                string = "{0}\n".format(cifog_headers[name])
                f.write(string)

            string = "{0} = {1}\n".format(name, cifog_params[name])
            f.write(string)

    return SAGE_fname, cifog_fname


def get_nion_fname(SAGE_params):
    """
    Using the ``fescPrescription`` specified in the ``SAGE.ini`` file,
    determines the name of the output ionizing photon file. 

    Parameters
    ----------

    SAGE_params: Dictionary
        Dictionary containing the ``SAGE.ini`` file parameters. 
    
    Returns
    ----------

    nion_fname: String
        Tag for the ionizing photon files that ``SAGE`` creates. 
    """

    fesc_prescription = SAGE_params["fescPrescription"]
    
    if fesc_prescription == 0:
        nion_fname = "{0}_fesc{1:.2f}_HaloPartCut{2}_nionHI" \
                      .format(SAGE_params["FileNameGalaxies"],
                              SAGE_params["beta"],
                              SAGE_params["HaloPartCut"])

    elif fesc_prescription == 1:
        nion_fname = "{0}_ejected_{1:.3f}_{2:.3f}_HaloPartCut{3}_nionHI" \
                     .format(SAGE_params["FileNameGalaxies"],
                             SAGE_params["alpha"],
                             SAGE_params["beta"],
                             SAGE_params["HaloPartCut"])

    elif fesc_prescription == 2:
        nion_fname = "{0}_quasar_{1:.2f}_{2:.2f}_{3:.2f}_HaloPartCut{4}_nionHI" \
                     .format(SAGE_params["FileNameGalaxies"],
                             SAGE_params["quasar_baseline"],
                             SAGE_params["quasar_boosted"],
                             SAGE_params["N_dyntime"],
                             SAGE_params["HaloPartCut"])

    elif fesc_prescription == 3 or fesc_prescription == 4:
        nion_fname = "{0}_AnneMH_{1:.3e}_{2:.2f}_{3:.3e}_{4:.2f}_HaloPartCut{5}_nionHI" \
                     .format(SAGE_params["FileNameGalaxies"],
                             SAGE_params["MH_low"],
                             SAGE_params["fesc_low"],
                             SAGE_params["MH_high"],
                             SAGE_params["fesc_high"],
                             SAGE_params["HaloPartCut"])

    elif fesc_prescription == 5: 
        nion_fname = "{0}_SFR_HaloPartCut{1}_nionHI" \
                     .format(SAGE_params["FileNameGalaxies"],
                             SAGE_params["HaloPartCut"])

    else:
        print("Select a valid fescPrescription, [0, 4].")
        raise ValueError

    return nion_fname


def make_ini_files(base_SAGE_ini, base_cifog_ini, 
                   SAGE_fields_update, cifog_fields_update,
                   run_directories):
    """
    Makes new ``SAGE`` and ``cifog`` ``.ini`` files for each run. 

    Parameters
    ----------

    base_SAGE_ini, base_cifog_ini : String
        Paths to the template ``SAGE`` and ``cifog`` ``.ini`` files.

    SAGE_fields_update, cifog_fields_update : Dictionaries of lists
        The ``SAGE`` and ``cifog`` parameter fields that will be updated for
        each run.  The values for each parameter field key are length Nx1 where
        N is the number of runs. 
        
    run_directories : List of strings, length equal to number of runs 
        Path to the base ``RSAGE`` directory for each run where all the model
        output will be placed.
 
    Returns
    ----------

    SAGE_ini_names, cifog_ini_names : List of strings, length equal to number
                                      of runs
        Paths to the ini file created for each run.
    """

    SAGE_ini_names = []
    cifog_ini_names = []

    # Now for each run, create a unique dictionary containing the fields for 
    # this run, update the ini files then create all the output directories. 
    for run_number in range(len(run_directories)):

        create_directories(run_directories[run_number])

        thisrun_SAGE_update = {}
        for name in SAGE_fields_update.keys():
            thisrun_SAGE_update[name] = SAGE_fields_update[name][run_number]

        thisrun_cifog_update = {}
        for name in cifog_fields_update.keys():
            thisrun_cifog_update[name] = cifog_fields_update[name][run_number]

        SAGE_fname, cifog_fname = update_ini_files(base_SAGE_ini, base_cifog_ini,
                                                   thisrun_SAGE_update, thisrun_cifog_update,
                                                   run_directories[run_number])        

        SAGE_ini_names.append(SAGE_fname)
        cifog_ini_names.append(cifog_fname)

    return SAGE_ini_names, cifog_ini_names


def make_slurm_files(base_slurm_file, SAGE_ini_names, cifog_ini_names, 
                     run_directories, Nproc): 
    """
    Makes ``slurm`` files for each run. 

    Parameters
    ----------

    base_slurm_file : String
        Path to the template slurm file.

    SAGE_ini_names, cifog_ini_names : List of strings, length equal to number
                                      of runs
        Paths to the ini file created for each run.

    run_directories : List of strings, length equal to number of runs 
        Path to the base ``RSAGE`` directory for each run where all the model
        output will be placed.

    Nproc : Integer
        Number of processors that each run will be executed with.
 
    Returns
    ----------

    slurm_names : List of strings, length equal to number
                                      of runs
        Paths to the ``slurm`` file created for each run.
    """
    slurm_names = []

    for run_number in range(len(SAGE_ini_names)):
        
        SAGE_params = ReadScripts.read_SAGE_ini(SAGE_ini_names[run_number])
        run_name = SAGE_params["FileNameGalaxies"]

        slurm_fname = "{0}/slurm_files/{1}.slurm".format(run_directories[run_number],
                                                         run_name) 

        tmp_slurm_fname = "{0}.tmp".format(base_slurm_file)
        copyfile(base_slurm_file, tmp_slurm_fname)

        # Want to replace lines in the slurm file. Set up the strings. 
        job_name = "#SBATCH --job-name={0}".format(run_name) 
        ntask = "#SBATCH --ntasks={0}".format(Nproc)
        NUMPROC = "NUMPROC={0}".format(Nproc)
        SAGE_ini = 'SAGE_ini="{0}"'.format(SAGE_ini_names[run_number])
        cifog_ini = 'cifog_ini="{0}"'.format(cifog_ini_names[run_number])
        run_prefix = 'run_prefix="{0}"'.format(run_name) 
        path_to_log = 'path_to_log="{0}/log_files/{1}.log"'.format(run_directories[run_number], run_name)

        # Replace strings at specific line numbers.
        line_numbers = [2, 4, 17, 19, 20, 24, 25]  
        string_names = [job_name, ntask, NUMPROC, SAGE_ini, cifog_ini, 
                        run_prefix, path_to_log]
        for line, name in zip(line_numbers, string_names): 
            command = "sed -i '{0}s@.*@{1}@' {2} ".format(line, name,
                                                          tmp_slurm_fname)
            subprocess.call(command, shell=True)

        # Finally move the temporary file to the final location.
        command = "mv {0} {1}".format(tmp_slurm_fname, slurm_fname)
        subprocess.call(command, shell=True)
        print("Created {0}".format(slurm_fname))

        slurm_names.append(slurm_fname)

    return slurm_names


def submit_slurm_jobs(slurm_names):
    """
    Submits the ``slurm`` jobs for each run. 

    Parameters
    ----------

    slurm_names : List of strings, length equal to number
                                      of runs
        Paths to the ``slurm`` file created for each run.
 
    Returns
    ----------

    None.
    """

    for slurm_fname in slurm_names:

        command = "sbatch {0}".format(slurm_fname)        
        subprocess.call(command, shell=True)
         
 
if __name__ == '__main__':

    #########################################################################
    # Specify here the SAGE parameters that you want to change for each run #
    #########################################################################

    fescPrescription = [5]
    #alpha = [0.30, 0.30, 0.30, 0.50, 0.50, 0.50]
    #beta = [0.0, 0.15, 0.20, 0.0, 0.15, 0.20]
    FileNameGalaxies = ["SFR"]

    SAGE_fields_update = { "fescPrescription" : fescPrescription,
                           #"alpha" : alpha,
                           # "beta" : beta,
                           "FileNameGalaxies" : FileNameGalaxies
                         }

    ##########################################################################
    # Specify here the cifog parameters that you want to change for each run #
    ##########################################################################

    cifog_fields_update = {}

    ################################################
    # Specify here the path directory for each run #
    ################################################

    run_directories = ["/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR"]

    #########################################################################
    # Specify here the path to the base ini files (shouldn't need to touch) #  
    #########################################################################

    base_SAGE_ini = "{0}/../ini_files/kali_SAGE.ini".format(script_dir)
    base_cifog_ini = "{0}/../ini_files/kali_cifog.ini".format(script_dir)
    base_slurm_file = "{0}/../run_rsage.slurm".format(script_dir)

    #####################
    # Misc for each run # 
    #####################
    
    Nproc = 32  # Number of processors to run on.

    #############
    # Functions # 
    #############

    check_input_parameters(SAGE_fields_update, cifog_fields_update,
                           run_directories, base_SAGE_ini, base_cifog_ini,
                           base_slurm_file)

    SAGE_ini_names, cifog_ini_names = make_ini_files(base_SAGE_ini, base_cifog_ini, 
                                                     SAGE_fields_update, cifog_fields_update,
                                                     run_directories)

    slurm_names = make_slurm_files(base_slurm_file, SAGE_ini_names, 
                                   cifog_ini_names, run_directories, Nproc)

    ###########################################################################
    # CAREFUL CAREFUL CAREFUL CAREFUL CAREFUL CAREFUL CAREFUL CAREFUL CAREFUL # 
    ###########################################################################
    # This function will submit all the jobs onto the queue.
    #submit_slurm_jobs(slurm_names)
