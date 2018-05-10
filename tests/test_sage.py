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

try: # Python2
    from urllib import urlretrieve
except ImportError: #Python3
    from urllib.request import urlretrieve
import subprocess

# Get the directory the testing happens in.
# Used a global variable for convenience as quite a few functions will use this.
test_dir = os.path.dirname(os.path.realpath(__file__))

scripts_dir = "{0}/../output/".format(test_dir)
sys.path.append(scripts_dir)

import AllVars
import ReadScripts


def get_trees():
    """
    Grabs the trees and galaxy output needed for the testing. 

    First checks the test directory to see if they trees are there.
    Otherwise downloads the test_RSAGE repo.
    
    Parameters
    ----------

    None.

    Returns
    ----------

    downloaded_repo: Boolean.  Required.
        Flag to check if we downloaded the test_RSAGE repo.  Useful for cleanup as we get 
        rid of all the unecessary files (e.g., ``README.rst``, ``LICENSE``, etc). 
    """

    print("")
    print("Checking to see if we need to download the test tree and output file.")

    tree_file = "{0}/trees_063_000.dat".format(test_dir)    
    if not os.path.isfile(tree_file):
        print("{0} does not exist, downloading the test_RSAGE repo and "
              "unzipping.".format(tree_file))
        master_branch = urlretrieve("https://github.com/jacobseiler/"
                                    "test_RSAGE/archive/master.zip")
        subprocess.call(["unzip", "-j", master_branch[0], "-d", test_dir])
        downloaded_repo = True 
    else:
        print("{0} exists so no need to download test_RSAGE repo."
              .format(tree_file))
        downloaded_repo = False

    print("Done")
    print("")

    return downloaded_repo

def check_sage_dirs(galaxy_name="test"):
    """
    """

    print("First checking that the required output directories are present.")
    print("")
    
    directory = "{0}/test_output/galaxies/".format(test_dir)
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:

        # Get rid of any old runs.
        print("Cleaning up old output files.")
        command = "rm {0}/test_output/galaxies/{1}_*".format(test_dir,
                                                             galaxy_name) 
        subprocess.call(command, shell=True)

    directory = "{0}/test_output/grids/".format(test_dir) 
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    print("Done.")
    print("")

def run_my_sage(ini_name="test_mini_millennium.ini"):

    print("Executing SAGE.")
    print("")

    path_to_sage = "{0}/../sage/sage".format(test_dir)
    path_to_ini = "{0}/test_ini_files/{1}".format(test_dir, ini_name)
    returncode = subprocess.call([path_to_sage, path_to_ini])

    if returncode != 0:
        print("SAGE exited with error code {0}".format(returncode))
        raise RuntimeError

    print("Done")
    print("")

def check_smf(galaxy_name="test"):

    print("")
    print("Now checking the stellar mass function for the final snapshot of "
          "mini-millennium.")

    AllVars.Set_Params_MiniMill()
    max_snap = len(AllVars.SnapZ) - 1

    # First check that the output of the test run can be read. 
   
    gal_name = "{0}/test_output/galaxies/{1}_z0.000".format(test_dir,
                                                            galaxy_name) 
    Gals, Gals_Desc = ReadScripts.ReadGals_SAGE(gal_name, 0, max_snap + 1)

    gal_name = "{0}/test_output/galaxies/{1}_MergedGalaxies".format(test_dir,
                                                                    galaxy_name) 
    Gals_Merged, _= ReadScripts.ReadGals_SAGE(gal_name, 0, max_snap + 1)
    Gals = ReadScripts.Join_Arrays(Gals, Gals_Merged, Gals_Desc)

    # Gals is now a recarray containing all galaxies at all snapshots. #

    w_gal = np.where((Gals.GridHistory[:, max_snap] != -1) & 
                     (Gals.GridStellarMass[:, max_snap] > 0.0))[0]
        
    position = Gals.GridHistory[w_gal, max_snap]
    w_wrong = np.where(position > pow(256,3))[0] 
    if (len(w_wrong) > 0):
        print("Found grid cells that had indices outside of the allowed "
              "values.")
        print("The indices must be between 0 and {0}".format(pow(256,3)))
        print("Galaxies {0} had indices {1}.".format(w_gal[w_wrong],
                                              position[w_wrong]))
        raise RuntimeError

    mass_test = np.log10(Gals.GridStellarMass[w_gal, max_snap] * 1.0e10 / AllVars.Hubble_h)   
    
    w_wrong = np.where(mass_test <= 0.0)[0] 
    if (len(w_wrong) > 0):
        print("The mass of the acceptable galaxies must be greater than 0.0.")
        print("Galaxies {0} had stellar mass {1}.".format(w_gal[w_wrong],
                                                          mass_test[w_wrong]))
        raise RuntimeError

    # Now let's check compare the mass of the test to the data. 

    mass_data = np.loadtxt("./mini_millennium_testmass.txt") 
   
    mass_difference = mass_test - mass_data
    w_wrong = np.where(mass_difference > 3e-3)[0]
    if (len(w_wrong) > 0):
        print("There must be no difference between the mass of the test run and"
              " the data in the test_RSAGE repository")
        print("Test Galaxies {0} had stellar mass {1} and data "
              "have stellar mass {2}".format(w_gal[w_wrong],
                                             mass_test[w_wrong],
                                             mass_data[w_wrong]))
        print("The difference is {0}".format(mass_difference[w_wrong]))
        raise RuntimeError
    

    print("")
    print("================")
    print("All tests passed")
    print("================")
    print("")

def test_run():
    """
    Wrapper to run all the tests.

    Parameters
    ----------

    None.

    Returns
    ----------

    None.
    """

    print("")
    print("Welcome to the RSAGE testing funhouse!")
    print("")

    downloaded_repo = get_trees() #  Download mini-millennium tree if we needed 

    ini_files = ["PhotonPrescription0_mini_millennium.ini", 
                 "PhotonPrescription1_mini_millennium.ini"]
    galaxy_names = ["PhotonPrescription0",
                    "PhotonPrescription1"] 

    for ini_file, galaxy_name in zip(ini_files, galaxy_names):

        check_sage_dirs(galaxy_name) # First check that directories for output are present. 
        run_my_sage(ini_file) #  Run my version of SAGE (not full R-SAGE yet).
        check_smf(galaxy_name) #  Attempt to make a stellar mass function.

    print("Done")
    print("")
    
    cleanup(downloaded_repo)

def cleanup(downloaded_repo):
    """
    Removes uncessary files from the directory. 

    Parameters
    ----------

    downloaded_repo: Boolean.  Required. 
        If we downloaded the entire test_RSAGE repo, need to get rid of files such as ``README.rst``.

    Returns
    ----------

    None.
    """

    print("Cleaning up files.")
   
    # If we downloaded the git repo, remove the uneeded files. 
    files = []   

    if downloaded_repo:
        subprocess.call(["rm", "README.rst", "LICENSE", ".gitignore"]) 

    # Remove all the output galaxy files.
   
    print("Done")
    print("")

if __name__ == "__main__":

    test_run()
