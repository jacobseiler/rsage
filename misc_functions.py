#!/usr/bin/env python
from __future__ import print_function

import sys
import numpy as np
import os
import argparse

sys.path.append('/home/jseiler/self_consistent_SAGE/output/')
import ReadScripts
import AllVars

import subprocess
def parse_input_arguments():
    """
    Parses the command line input arguments.

    If there has not been a SAGE .ini or cifog .ini file specified a 
    RuntimeError will be raised. 

    Parameters
    ----------

    None.

    Returns
    ----------

    args: Dictionary. Required.
        Dictionary of arguments from the ``argparse`` package.
        Dictionary is keyed by the argument name (e.g., args["SAGE_fname"]).
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--SAGE_ini", dest="SAGE_fname", 
                        help="Location of the SAGE ini file.  Required.")
    parser.add_argument("-c", "--cifog_ini", dest="cifog_fname", 
                        help="Location of the cifog ini file.")
    parser.add_argument("-p", "--precision", dest="precision", 
                        help="Precision of the grid files. 0 for int, 1 for " 
                        "float, 2 for double. Default: 2", default=2, 
                        type=int) 
    parser.add_argument("-r", "--reionredshift", dest="reionredshift", 
                        help="Set to 1 to generate the reionization redshift "
                        "grid for the snapshot range specified. Default: 0.", default = 0, 
                        type=int)
    parser.add_argument("-i", "--ini_dir", dest="ini_dir", 
                        help="Specifies the directory to create two new .ini "
                        "files that are identical to the specified .ini files "
                        "but with value of HighSnap and stopSnapshot "
                        "incremented by 1.  Default: None", 
                        default = None)
    parser.add_argument("-x", "--run_prefix", dest="prefix",
                        help="Prefix for the naming of files. Useful to run "
                            "multiple models with only slight variations. "
                            "Default: No prefixes.", default = None)

    args = parser.parse_args()

    if (args.SAGE_fname is None):
        parser.print_help()
        raise RuntimeError 

    if (args.ini_dir is not None and args.cifog_fname is None):
        print("To increment the ini files both a SAGE and cifog ini files must "
              "be specified.")
        parser.print_help() 
        raise RuntimeError

    return vars(args) 


def create_redshift_grid(args):
    """
    Creates a grid containing the redshift each grid cell was ionized at.

    Parameters
    ----------

    args: Dictionary. Required.
        Dictionary containing the input parameters specified at runtime.

    Returns
    ----------

    None. The new reionization_redshift grid is written to the directory
    specified by ``SAGE_params["PhotoIonDir"]``. 
    """

    SAGE_params, SAGE_params_names = ReadScripts.read_SAGE_ini(args["SAGE_fname"])

    SnapList = np.arange(SAGE_params["LowSnap"] + 1, 
                         SAGE_params["HighSnap"] + 2) 

    AllVars.Set_Params_Kali()

    GridSize = SAGE_params["GridSize"][0]
    xHII_base = "{0}/{1}_XHII".format(SAGE_params["PhotoionDir"][0],
                                      args["prefix"])     
    redshift_output_base = "{0}/{1}".format(SAGE_params["PhotoionDir"][0],
                                                SAGE_params["ReionRedshiftName"][0]) 

    reionization_redshift_grid = np.full((pow(GridSize, 3)), -1.0)

    for snapshot in SnapList:
        xHII_fname = "{0}_{1:03d}".format(xHII_base, snapshot)
        xHII_grid = ReadScripts.read_binary_grid(xHII_fname, GridSize,
                                                 args["precision"], False) 
        
        w_ionized_snap = np.where((xHII_grid > 0.9))[0] # Indices of the ionized cells.
        w_not_ionized = np.where((reionization_redshift_grid < -0.5))[0] # Indices of the cells that have yet to be ionized. 
    
        w_to_update = np.intersect1d(w_ionized_snap, w_not_ionized) 

        print("")
        print("For snapshot {0} there were {1} cells already ionized, {2} cells ionized in this snapshot resulting in {3} cells we are about to ionized".format(snapshot, pow(GridSize, 3) - len(w_not_ionized), len(w_ionized_snap), len(w_to_update))) 
        print("{0:.2f}% of cells will be ionized after this snapshot.".format((pow(GridSize, 3) - len(w_not_ionized) + len(w_to_update)) / pow(GridSize, 3) * 100.0))
        print("")
        reionization_redshift_grid[w_to_update] = AllVars.SnapZ[snapshot]
        
    fname_out = "{0}".format(redshift_output_base) 
    reionization_redshift_grid.tofile(fname_out)

    print("Reionization redshift grid saved to file {0}".format(fname_out))


def increment_ini(args):
    """
    Rewrites a new SAGE and cifog .ini file with HighSnap and stopSnapshot
    incremented by 1. 

    This is done to run the next iteration of R-SAGE.

    Parameters
    ----------

    args: Dictionary.  Required.
        Dictionary containing the input parameters specified at runtime.

    Returns
    ----------

    None. The new ini files are written to the directory specified by
    ``args["ini_dir"]``.
    """

    if args["prefix"] is None:
        prefix_tag = ""
    else:
        prefix_tag = args["prefix"]

    cifog_params, cifog_params_names, cifog_headers = ReadScripts.read_cifog_ini(args["cifog_fname"])
    SAGE_params, SAGE_params_names = ReadScripts.read_SAGE_ini(args["SAGE_fname"])

    concat_properties(SAGE_params)

    SAGE_params['HighSnap'] = [SAGE_params['HighSnap'][0] + 1] 
    SAGE_params['ReionSnap'] = [SAGE_params['ReionSnap'][0] + 1] 

    cifog_params["stopSnapshot"] = [cifog_params["stopSnapshot"][0] + 1]
 
    fname_SAGE = "{0}/{2}_SAGE_snap{1}.ini".format(args["ini_dir"],
                                                   SAGE_params["HighSnap"][0],
                                                   prefix_tag)
    fname_cifog = "{0}/{2}_cifog_snap{1}.ini".format(args["ini_dir"],
                                                     cifog_params["stopSnapshot"][0],
                                                     prefix_tag)
   
    with open (fname_SAGE, "w+") as f:
        for name in SAGE_params_names:
            string = "{0} {1}\n".format(name, SAGE_params[name][0])
            f.write(string)

    with open (fname_cifog, "w+") as f:
        for name in cifog_params_names:
            if name in cifog_headers:
                header_string = "{0}".format(cifog_headers[name])
                f.write(header_string)
            string = "{0} = {1}\n".format(name, cifog_params[name][0])
            f.write(string)

    print("Successfully created ini file {0}".format(fname_SAGE))
    print("Successfully created ini file {0}".format(fname_cifog))


def concat_properties(SAGE_params):

    command = "cat {0}/properties/{2}_misc_properties_{1:03d}_*" \
                 .format(SAGE_params["GridOutputDir"][0],
                 SAGE_params["HighSnap"][0],
                 SAGE_params["FileNameGalaxies"][0])

    fname = "{0}/properties/{2}_misc_properties_{1:03d}" \
                 .format(SAGE_params["GridOutputDir"][0],
                 SAGE_params["HighSnap"][0],
                 SAGE_params["FileNameGalaxies"][0])
    
    f = open(fname, "w")
    subprocess.call(command, stdout=f, shell=True)
    f.close() 

    command = "rm {0}/properties/misc_properties_{1:03d}_*" \
                 .format(SAGE_params["GridOutputDir"][0],
                 SAGE_params["HighSnap"][0])

    subprocess.call(command, shell=True)

if __name__ == '__main__':

    args = parse_input_arguments()
               
    if args["reionredshift"]:
        create_redshift_grid(args)

    if args["ini_dir"] is not None:
        increment_ini(args)

