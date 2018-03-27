#!/usr/bin/env python

"""
This file controls testing for the R-SAGE model.

As R-SAGE is run iteratively using a number of separate codes, there are
multiple tiers of testing.j

- 0th Order -
Codes actually run and produce outputs that can be read.

For SAGE: This involves checking that the outputs can produce a stellar 
mass function. 
"""
 
from __future__ import print_function
import numpy as np
import argparse
import sys
import h5py
import os
import pytest

from urllib.request import urlretrieve
import subprocess

# Get the directory the testing happens in.
# Used a global variable for convenience as quite a few functions will use this.
test_dir = os.path.dirname(os.path.realpath(__file__))


def get_trees():
    """
    Grabs the trees needed for the testing. 

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
    print("Checking to see if we need to download the test tree file.")

    tree_file = "{0}/trees_063.0".format(test_dir)
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


def run_my_sage():
    """
    """

    print("Now running my version of SAGE (not full R-SAGE yet).")
    print("")

    print("First checking that the required output directories are present.")
    
    directory = "{0}/test_output/galaxies/".format(test_dir) 
    if not os.path.exists(directory):
        os.makedirs(directory)

    directory = "{0}/test_output/grids/".format(test_dir) 
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    print("Done.")
    print("Executing SAGE.")    

    subprocess.call(["../sage/sage",
                     "./test_ini_files/test_mini_millennium.ini"])

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

    downloaded_repo = get_trees() #  Download mini-millennium tree if we need to. 

    run_my_sage() #  Run my version of SAGE (not full R-SAGE yet).

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
   
    if downloaded_repo:
        subprocess.call(["rm", "README.rst", "LICENSE", ".gitignore"]) 

    print("Done")
    print("")

if __name__ == "__main__":

    test_run()
