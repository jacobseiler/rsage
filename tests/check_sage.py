#!/usr/bin/env python

"""
This file controls testing for the R-SAGE model.

As R-SAGE is run iteratively using a number of separate codes, there are
multiple tiers of testing.

- 0th Order -
Codes actually run and produce outputs that can be read.

For SAGE: This involves checking that the outputs can produce a stellar 
mass function. 
"""
 
from __future__ import print_function
import numpy as np
import argparse
import sys
import os
import pytest
from tqdm import tqdm

# Get the directory the testing happens in.
# Used a global variable for convenience as quite a few functions will use this.
test_dir = os.path.dirname(os.path.realpath(__file__))

scripts_dir = "{0}/../output/".format(test_dir)
sys.path.append(scripts_dir)

import AllVars
import ReadScripts
import PlotScripts


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

    parser.add_argument("-f", "--SAGE_file", dest="SAGE_fname", 
                        help="Location of the SAGE file to check.  Required.")

    parser.add_argument("-s", "--simulation", dest="simulation", 
                        help="Name of the simulation we're checking. "
                             "Default:Kali.  Accepted: Kali, MiniMill.",
                        default="Kali")

    args = parser.parse_args()

    if (args.SAGE_fname is None): 
        parser.print_help()
        raise RuntimeError 

    return vars(args) 


def check_sage_file(SAGE_fname, simulation="Kali"):
    """
    Goes through the specified SAGE file and does a number of checks.

    Parameters
    ----------

    SAGE_fname: String.  Required.
        Name of the SAGE file to check.

    simulation: String.  Default: Kali.
        Name of the simulation SAGE used. 

    Returns
    ----------

    None.  If a check fails, the program will exit.
    """

    if simulation == "Kali":
        AllVars.Set_Params_Kali()
    elif simulation == "MiniMill":
        AllVars.Set_Params_MiniMill()
    else:
        print("A valid simulation has not been selected when calling "
              "`check_sage_file.")
        raise ValueError

    Gals, Gals_Desc = ReadScripts.ReadGals_SAGE(SAGE_fname, None, 
                                                len(AllVars.SnapZ)) 

    NTrees = ReadScripts.Read_SAGE_header(SAGE_fname, None)

    check_centrals(Gals, NTrees)


def check_centrals(Gals, NTrees):
    """
    Ensures that there is only one central for each FoF Halo. 

    Parameters
    ----------

    Gals: numpy structured array.  Required.
        Galaxies that will be checked.
        See `ReadGals_SAGE` in `output/ReadScripts.py` for full data type
        description. 

    NTrees: Integer.  Required.
        Number of trees in the file. 

    Returns
    ----------

    None.  If a check fails, the program will exit.
    """

    for treenr in tqdm(range(NTrees)): 
        w_gal = np.where(Gals.TreeNr == treenr)[0]

        if len(w_gal) == 0:
            continue
   
        for snapshot in range(len(AllVars.SnapZ)):            
            w_alive = w_gal[np.where((Gals.GridHistory[w_gal, snapshot] != -1))[0]]
                
            FoFNr = Gals.GridFoFHaloNr[w_alive, snapshot]
            unique_FoFNr = np.unique(FoFNr)
            for fof in unique_FoFNr:
                w_fof = w_alive[np.where(FoFNr == fof)[0]]
                   
                unique, unique_count = np.unique(Gals.GridType[w_fof,
                                                 snapshot], return_counts = True)

                if (unique_count[0] == 0 and unique_count[0] > 1):
                    print("For Tree {0}, FoFHaloNr {1} (existing at snapshot "
                          "{2}) had {3} centrals.  There should only ever be " 
                          "1 central per FoF Halo.".format(treenr, fof,
                                                           snapshot,
                                                           unique_count[0]))
                    print("FoFNr: {0}".format(Gals.GridFoFHaloNr[w_fof,
                                              snapshot]))
                    print("Type: {0}".format(Gals.GridType[w_fof, snapshot]))
                    raise ValueError

if __name__ == "__main__":

    args = parse_input_arguments()
    check_sage_file(args["SAGE_fname"], args["simulation"])
