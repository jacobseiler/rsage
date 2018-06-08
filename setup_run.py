"""
Creates directories and updates ini files with the correct file names.

The user specifies the base output directory and the .ini files for SAGE and
cifog.

The script then creates directories such as '<BaseDir>/galaxies'.  It then
reads the escape fraction prescription and constants, and determines what the
eventual name for the output nion files will be. Finally it updates the .ini
files with these nion filepaths.
"""

#!/usr/bin/env python
from __future__ import print_function

import numpy as np
import sys
import os
import argparse
import subprocess

sys.path.append('./output/')
import ReadScripts
import AllVars

def parse_input_arguments():
    """
    Parses the command line input arguments.

    If there has not been a runtime directory, SAGE .ini or cifog .ini file
    specified a RuntimeError will be raised. 

    Parameters
    ----------

    None.

    Returns
    ----------

    args: Dictionary.  Required.
        Dictionary of arguments from the ``argparse`` package.
        Dictionary is keyed by the argument name (e.g., args["SAGE_fname"]).
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", dest="run_directory", 
                        help="Path to the directory the output files will be "
                             "located.  Required.  Enter WITHOUT the final /")

    parser.add_argument("-f", "--SAGE_ini", dest="SAGE_fname", 
                        help="Location of the SAGE ini file.  Required.")

    parser.add_argument("-c", "--cifog_ini", dest="cifog_fname", 
                        help="Location of the cifog ini file.  Required.")

    parser.add_argument("-p", "--run_prefix", dest="prefix",
                        help="Prefix for the naming of files. Useful to run "
                            "multiple models with only slight variations. "
                            "Default: No prefixes.", default = None)

    args = parser.parse_args()

    if (args.SAGE_fname is None or args.run_directory is None or
        args.cifog_fname is None):
        parser.print_help()
        raise RuntimeError 

    return vars(args) 


def create_directories(args):
    """
    Creates the directories to house all the R-SAGE outputs.   
 
    Parameters
    ----------

    args: Dictionary.  Required.
        Dictionary containing the input parameters specified at runtime.

    Returns
    ----------
    
    None.
    """

    # First create the base directory.
    base_dir = args["run_directory"]
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
                "grids/nion/properties"]
        for directory in dirs:
            dir_ = "{0}/{1}".format(base_dir, directory)
            os.makedirs(dir_)
            print("Created directory {0}".format(dir_))

   
def update_ini_files(args):
    """
    Rewrites the ini files to point to the correct output files.  

    Parameters
    ----------

    args: Dictionary.  Required.
        Dictionary containing the input parameters specified at runtime.
    
    Returns
    ----------

    None.
    """

    SAGE_params, SAGE_params_names = ReadScripts.read_SAGE_ini(args["SAGE_fname"])
    cifog_params, cifog_params_names, cifog_headers = ReadScripts.read_cifog_ini(args["cifog_fname"])
 
    if args["prefix"] is None:
        prefix_tag = ""
    else:
        prefix_tag = args["prefix"]
 
    SAGE_params["OutputDir"] = "{0}/galaxies".format(args["run_directory"])
    SAGE_params["GridOutputDir"] = "{0}/grids/nion".format(args["run_directory"])
    SAGE_params["PhotoionDir"] = "{0}/grids/cifog".format(args["run_directory"])
    SAGE_params["PhotoionName"] = "{0}_photHI".format(prefix_tag)
    SAGE_params["ReionRedshiftName"] = "{0}_reionization_redshift" \
                                       .format(prefix_tag)

    nion_fname = get_nion_fname(SAGE_params) 
    cifog_params["inputNionFile"] = "{0}/grids/nion/{1}" \
                                    .format(args["run_directory"], nion_fname)
    cifog_params["output_XHII_file"] = "{0}/grids/cifog/{1}_XHII" \
                                       .format(args["run_directory"],
                                               prefix_tag)
    cifog_params["output_photHI_file"] = "{0}/grids/cifog/{1}_photHI" \
                                         .format(args["run_directory"],
                                                 prefix_tag)
    cifog_params["output_restart_file"] = "{0}/grids/cifog/{1}_restart" \
                                          .format(args["run_directory"],
                                                  prefix_tag)
    print("Nion_fname {0}".format(nion_fname))

    with open (args["SAGE_fname"], "w+") as f:
        for name in SAGE_params_names:
            string = "{0} {1}\n".format(name, SAGE_params[name][0])
            f.write(string)

    with open (args["cifog_fname"], "w+") as f:
        for name in cifog_params_names:
            if name in cifog_headers:
                header_string = "{0}".format(cifog_headers[name])
                f.write(header_string)
            string = "{0} = {1}\n".format(name, cifog_params[name][0])
            f.write(string)


def get_nion_fname(SAGE_params):
    """
    Using the fescPrescription specified in the SAGE.ini file, determines the
    name of the output nion files. 

    NOTE: fescPrescription == 1 is deprecated and is not allowed. If the .ini
    file uses this fescPrescription a ValueError will be raised.

    Parameters
    ----------

    SAGE_params: Dictionary.  Required.
        Dictionary containing the SAGE .ini file parameters. 
    
    Returns
    ----------

    nion_fname: String. Required.
        Base name of the eventual output nion files. 
    """

    fesc_prescription = SAGE_params["fescPrescription"]
    
    if fesc_prescription == 0:
        nion_fname = "{0}_fesc{1}_HaloPartCut{2}_nionHI" \
                      .format(SAGE_params["FileNameGalaxies"][0],
                              SAGE_params["fesc"][0],
                              SAGE_params["HaloPartCut"][0])

    elif fesc_prescription == 1:
        print("Using fesc_prescription of 1 is deprecated.")
        raise ValueError 

    elif fesc_prescription == 2:
        alpha, beta = determine_fesc_constants(SAGE_params)              
        nion_fname = "{0}_MH_{1:.3e}_{2:.2f}_{3:.3e}_{4:.2f}_HaloPartCut{5}_nionHI" \
                     .format(SAGE_params["FileNameGalaxies"][0],
                             SAGE_params["MH_low"][0],
                             SAGE_params["fesc_low"][0],
                             SAGE_params["MH_high"][0],
                             SAGE_params["fesc_high"][0],
                             SAGE_params["HaloPartCut"][0])

    elif fesc_prescription == 3:
        nion_fname = "{0}_ejected_{1:.3f}_{2:.3f}_HaloPartCut{3}_nionHI" \
                     .format(SAGE_params["FileNameGalaxies"][0],
                             SAGE_params["alpha"][0],
                             SAGE_params["beta"][0],
                             SAGE_params["HaloPartCut"][0])

    elif fesc_prescription == 4:
        nion_fname = "{0}_quasar_{1:.2f}_{2:.2f}_{3:.2f}_HaloPartCut{4}_nionHI" \
                     .format(SAGE_params["FileNameGalaxies"][0],
                             SAGE_params["quasar_baseline"][0],
                             SAGE_params["quasar_boosted"][0],
                             SAGE_params["N_dyntime"][0],
                             SAGE_params["HaloPartCut"][0])

    elif fesc_prescription == 5 or fesc_prescription == 6:
        nion_fname = "{0}_AnneMH_{1:.3e}_{2:.2f}_{3:.3e}_{4:.2f}_HaloPartCut{5}_nionHI" \
                     .format(SAGE_params["FileNameGalaxies"][0],
                             SAGE_params["MH_low"][0],
                             SAGE_params["fesc_low"][0],
                             SAGE_params["MH_high"][0],
                             SAGE_params["fesc_high"][0],
                             SAGE_params["HaloPartCut"][0])

    elif fesc_prescription == 7:
        nion_fname = "{0}_ejectedpower_{1:.3e}_{2:.2f}_{3:.3e}_{4:.2f}_HaloPartCut{5}_nionHI" \
                     .format(SAGE_params["FileNameGalaxies"][0],
                             SAGE_params["MH_low"][0],
                             SAGE_params["fesc_low"][0],
                             SAGE_params["MH_high"][0],
                             SAGE_params["fesc_high"][0],
                             SAGE_params["HaloPartCut"][0])
        
    elif fesc_prescription == 8:
        nion_fname = "{0}_mstar_{1:.3e}_{2:.3e}_{3:.2f}_{4:.2f}_HaloPartCut{5}_nionHI" \
                     .format(SAGE_params["FileNameGalaxies"][0],
                             SAGE_params["fesc_Mstar_low"][0],
                             SAGE_params["fesc_Mstar_high"][0],
                             SAGE_params["fesc_Mstar"][0],
                             SAGE_params["fesc_not_Mstar"][0],
                             SAGE_params["HaloPartCut"][0])

    elif fesc_prescription == 9:
        nion_fname = "{0}_ejectedSN_{1:.3f}_{2:.3f}_HaloPartCut{3}_nionHI" \
                     .format(SAGE_params["FileNameGalaxies"][0],
                             SAGE_params["alpha"][0],
                             SAGE_params["beta"][0],
                             SAGE_params["HaloPartCut"][0])

    elif fesc_prescription == 10:
        nion_fname = "{0}_ejectedQSO_{1:.3f}_{2:.3f}_HaloPartCut{3}_nionHI" \
                     .format(SAGE_params["FileNameGalaxies"][0],
                             SAGE_params["alpha"][0],
                             SAGE_params["beta"][0],
                             SAGE_params["HaloPartCut"][0])

    else:
        print("Select a valid fescPrescription (0 to 7 inclusive).")
        raise ValueError

    return nion_fname


def determine_fesc_constants(SAGE_params):
    """
    If the fescPrescription depends on Halo mass, the functional form is fesc =
    alpha*MH^(beta).  This function determines the values of alpha and beta
    depending upon the fixed points specified in the SAGE parameter file.

    If the fixed points have had their halo mass specified in log units a
    ValueError will be raised.

    Parameters
    ----------

    SAGE_params: Dictionary.  Required.
        Dictionary containing the SAGE .ini file parameters. 
    
    Returns
    ----------

    alpha, beta: Floats. Required.
        Constants for fesc equation. 
    """

    # The SAGE ini file specifies two fixed points.
    fesc_high = SAGE_params["fesc_high"][0]
    MH_high = SAGE_params["MH_low"][0]

    fesc_low = SAGE_params["fesc_low"][0]    
    MH_low = SAGE_params["MH_high"][0]

    # The values for halo mass should be in non-log units. Do a quick check.
    if (MH_high < 1e6 or MH_low < 1e6):
        print("If using fescPrescription == 2 (fesc depends on halo mass) the "
              "fixed points need to have their halo mass specified in Msun, NOT "
              "LOG MSUN.")
        raise ValueError 
    

    log_A = (np.log10(fesc_high) - np.log10(fesc_low)*np.log10(MH_high) / np.log10(MH_low)) \
            * pow(1 - np.log10(MH_high) / np.log10(MH_low), -1)

    B = (np.log10(fesc_low) - log_A) / np.log10(MH_low);
    A = pow(10, log_A);

    alpha = A;
    beta = B;

    return alpha, beta


if __name__ == '__main__':

    args = parse_input_arguments()
    
    create_directories(args)
    update_ini_files(args)
