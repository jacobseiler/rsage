#!/usr/bin/env python

"""
This file controls testing for the R-SAGE model.

As R-SAGE is run iteratively using a number of separate codes, there are
multiple tiers of testing.

For comparison, the results of the test models are compared to the data in the
'test_RSAGE' directory.

- 0th Order -
Codes actually run and produce outputs that can be read.

For SAGE we check the stellar masses of the galaxies.  Since the output format
of SAGE can change between versions, we store the 'correct' stellar mass as a
.txt file in the 'test_RSAGE' directory.  See `check_smf()` for list of
conditions that fail the test.
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

# This is the directory name that will contain all the **correct** data from
# the test repository.
test_datadir = "rsage_testdata"

scripts_dir = "{0}/../output/".format(test_dir)
sys.path.append(scripts_dir)

import AllVars
import ReadScripts


def get_trees():
    """
    Grabs the trees and galaxy output needed for the testing. 

    First checks the test directory to see if they trees are there.
    Otherwise downloads the ``test_datadir`` repo.
    
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

    treefile = "{0}/{1}/kali512/trees/kali512_000.dat".format(test_dir, test_datadir)
    if not os.path.isfile(treefile):
        print("{0} does not exist, downloading the {1}"
              "repo.".format(treefile, test_datadir)) 

        repo_url = "https://github.com/jacobseiler/{0}".format(test_datadir)
        command = "git clone {0} --depth 1".format(repo_url)
        subprocess.call(command, shell=True)

        downloaded_repo = True

    else:
        print("{0} exists so no need to download {1} repo."
              .format(treefile, test_datadir))
        downloaded_repo = False

    print("Done")
    print("")

    return downloaded_repo


def check_sage_dirs(galaxy_name="test"):
    """
    Ensures the directories for outputs exist. 

    Also removes old files that may be in the directory.

    Parameters
    ----------

    galaxy_name: String. Optional, default: 'test'. 
        Prefix name for the galaxies. 

    Returns
    ----------

    None.
    """

    print("First checking that the required output directories are present.")
    print("")
   
    # Check to see if the output directory exists.
    # If it does, remove any old output.
    # Otherwise, create the directory.  
    directory = "{0}/test_output/galaxies/".format(test_dir)
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        print("Cleaning up old output files.")
        command = "rm {0}/test_output/galaxies/{1}_*".format(test_dir,
                                                             galaxy_name) 
        subprocess.call(command, shell=True)

    # Check that the log directory exists.
    directory = "{0}/test_logs/".format(test_dir)
    if not os.path.exists(directory):
        os.makedirs(directory)

    print("Done.")
    print("")


def run_my_sage(ini_name):
    """
    Executes my version of SAGE. 

    This uses the base analytic prescription for reionization (i.e., no
    self-consistent). 

    Parameters
    ----------

    ini_name: String
        Name of the input ini file. 

    Returns
    ----------

    None.

    If the execution of SAGE fails (exit code is not 0), a RuntimeError is
    raised.
    """

    print("Executing SAGE.")
    print("")

    # Use paths relative to the file this script is in.
    path_to_sage = "{0}/../bin/rsage".format(test_dir)
    path_to_sage_ini = "{0}/test_ini_files/{1}_SAGE.ini".format(test_dir, ini_name)

    path_to_cifog_ini = "{0}/test_ini_files/{1}_cifog.ini".format(test_dir, ini_name)
    command = "{0} {1} {2}".format(path_to_sage, path_to_sage_ini, path_to_cifog_ini)

    returncode = subprocess.call(command, shell=True)

    if returncode != 0:
        print("SAGE exited with error code {0}".format(returncode))
        raise RuntimeError

    print("Done")
    print("")


def check_smf(Gals, galaxy_name, max_snap, sage_params, mass_tol=3.0e-3,
              update_data=0):
    """
    Checks the stellar mass of the galaxies.

    The galaxies of the model being tested are loaded in and compared to the
    'correct' masses from the 'test_RSAGE' repo. 

    Parameters
    ----------

    mass_tol: Float. Optional, default: 3.0e-3.
        The allowed difference between the stellar mass of the test galaxies
        and the stellar mass in the 'test_RSAGE' repo.  If any galaxy has a
        (absolute) difference larger than this, a RuntimeError is raised. 

    Returns
    ----------

    None.

    If any test galaxy has a GridPosition greater than 256^3 a RuntimeError is
    raised. 

    If test model produces a different number of galaxies to the 'test_RSAGE'
    galaxies a RuntimeError is raised. 

    If the mass of any test galaxy is negative a RuntimeError is raised.

    If the mass of any test galaxy differs by more than ``mass_tol`` dex 
    compared to the 'test_RSAGE' galaxies a RuntimeError is raised. 
    """

    GridSize = int(sage_params["GridSize"])
    Hubble_h = float(sage_params["Hubble_h"])
    LastSnapshot = int(sage_params["LastSnapShotNr"])

    print("")
    print("Now checking the stellar mass function for the final snapshot of "
          "Kali.")

    w_gal = np.where((Gals.GridHistory[:, max_snap] != -1) & 
                     (Gals.GridStellarMass[:, max_snap] > 0.0))[0]
        
    position = Gals.GridHistory[w_gal, max_snap]
    w_wrong = np.where(position > pow(GridSize, 3))[0] 
    if (len(w_wrong) > 0):
        print("Found grid cells that had indices outside of the allowed "
              "values.")
        print("The indices must be between 0 and {0}".format(pow(GridSize, 3)))
        print("Galaxies {1} had indices {1}.".format(w_gal[w_wrong],
                                              position[w_wrong]))
        raise RuntimeError

    mass_test = np.log10(Gals.GridStellarMass[w_gal, max_snap] * 1.0e10 / Hubble_h)   
    mass_test_nolog = Gals.GridStellarMass[w_gal, max_snap] * 1.0e10 / Hubble_h 
    
    w_wrong = np.where(mass_test_nolog <= 0.0)[0] 
    if (len(w_wrong) > 0):
        print("The mass of the acceptable galaxies must be greater than 0.0.")
        print("Galaxies {0} had stellar mass {1}.".format(w_gal[w_wrong],
                                                          mass_test_nolog[w_wrong]))
        raise RuntimeError

    if update_data:
        print("=======================================================")
        print("WARNING WARNING WARNING WARNING WARNING WARNING WARNING")
        print("=======================================================")
        print("=======================================================")
        print("WARNING WARNING WARNING WARNING WARNING WARNING WARNING")
        print("=======================================================")
        input("YOU ARE ABOUT TO OVERWRITE THE TEST DATA. IF THIS IS WHAT YOU"
              "WANT, PRESS ENTER OTHERWISE CTRL-C TO GET OUTTA HERE!")
        input("JUST CHECKING ONCE MORE!")

        fname = "{0}/{1}/kali512/data/{2}_testmass.txt".format(test_dir, test_datadir,
                                                               galaxy_name)
        np.savetxt(fname, mass_test)
        print("Saved mass data as {0}".format(fname))

        print("All test data updated. Exiting checking now (because it'll "
              "obviously be correct.)")
        return

    # Now let's check compare the mass of the test to the 'test_RSAGE' repo. 
    mass_repo = np.loadtxt("{0}/{1}/kali512/data/{2}_testmass.txt".format(test_dir,
                                                                          test_datadir,
                                                                          galaxy_name))

    if len(mass_test) != len(mass_repo):
        print("For the test data we had {0} galaxies at Snapshot {2} with > 0 "
              "Stellar Mass.  This is compared to {1} galaxies for the "
              "CORRECT mass data.".format(len(mass_test), len(mass_data),
                                          LastSnapshot))
        print("First 10 galaxies for the test data {0}".format(mass_test[0:10]))
        print("First 10 galaxies for the repo data {0}".format(mass_repo[0:10]))

        raise RuntimeError
 
    mass_difference = mass_test - mass_repo
    w_wrong = np.where(abs(mass_difference) > mass_tol)[0]
    if (len(w_wrong) > 0):
        print("There must be no difference between the mass of the test run and"
              " the data in the test_RSAGE repository")
        print("Test Galaxies {0} had stellar mass {1} and repo "
              "have stellar mass {2}".format(w_gal[w_wrong],
                                             mass_test[w_wrong],
                                             mass_repo[w_wrong]))
        print("The difference is {0}".format(mass_difference[w_wrong]))
        print("Out of a total {0} galaxies, {1} had a mass difference greater "
              "than the tolerance level of {2}".format(len(mass_difference),
                                                       len(w_wrong),
                                                       mass_tol))
        print("The average mass difference of these galaxies is "
              "{0} with an overall average mass difference of {1}" \
              .format(np.mean(abs(mass_difference[w_wrong])),
                      np.mean(abs(mass_difference))))

        max_diff = np.argmax(abs(mass_difference))
        
        print("The largest mass difference is {0} which corresponds to a {1}% " 
              "shift compared to the test data (mass {2})." \
              .format(mass_difference[max_diff],
                      abs(mass_difference[max_diff]/mass_test[max_diff]*100.0),
                      mass_test[max_diff]))
        raise RuntimeError


def load_gals(max_snap, galaxy_name="test"):

    # First check that the output of the test run can be read. 

    gal_name = "{0}/test_output/galaxies/{1}_z5.829".format(test_dir,
                                                            galaxy_name) 
    Gals, Gals_Desc = ReadScripts.ReadGals_SAGE(gal_name, 0, max_snap + 1)

    gal_name = "{0}/test_output/galaxies/{1}_MergedGalaxies".format(test_dir,
                                                                    galaxy_name) 
    Gals_Merged, _= ReadScripts.ReadGals_SAGE(gal_name, 0, max_snap + 1)
    Gals = ReadScripts.Join_Arrays(Gals, Gals_Merged, Gals_Desc)

    # Gals is now a recarray containing all galaxies at all snapshots. 

    return Gals


def read_grids(SAGE_params, cifog_params, snapshot, reshape=False):
    """
    Reads the grids for a specific ``RSAGE`` run.

    Parameters
    ----------

    SAGE_params : Dictionary
        Dictionary keyed by the ``SAGE`` parameter field names and containing 
        the values from the ``.ini`` file. See ``read_SAGE_ini`` in
        ``output/ReadScripts.py`` for full details.

    cifog_params : Dictionary
        Dictionary keyed by the ``cifog`` parameter field names and containing 
        the values from the ``.ini`` file. See ``read_cifog_ini`` in
        ``output/ReadScripts.py`` for full details.

    reshape : Boolean, default False
        Controls whether the grids should be recast to NxNxN arrays (True) or kept as 1D
        arrays (False).
    """

    print("Reading the nionHI, XHII and photHI grids for the test run.")

    RunPrefix = SAGE_params["RunPrefix"]
    RunPrefix = SAGE_params["RunPrefix"]
    OutputDir = SAGE_params["OutputDir"]
    GridSize = int(SAGE_params["GridSize"])

    # Here the precision is 1 for float, 2 for double. XHII and photHI are
    # hardcoded to be in double.

    nion_prefix = get_nion_prefix(SAGE_params)
    nion_path = "{0}/grids/nion/{1}_{2}_nionHI_{3:03d}".format(OutputDir, RunPrefix,
                                                       nion_prefix, snapshot)
    nion_precision = int(cifog_params["nionFilesAreInDoublePrecision"]) + 1
    nion_grid = ReadScripts.read_binary_grid(nion_path, GridSize,
                                             nion_precision, reshape=reshape)

    XHII_path = "{0}/grids/cifog/{1}_XHII_{2:03d}".format(OutputDir, RunPrefix,
                                                          snapshot)
    XHII_precision = 2
    XHII_grid = ReadScripts.read_binary_grid(XHII_path, GridSize,
                                             XHII_precision, reshape=reshape)

    photHI_path = "{0}/grids/cifog/{1}_photHI_{2:03d}".format(OutputDir,
                                                              RunPrefix,
                                                              snapshot)

    photHI_precision = 2
    photHI_grid = ReadScripts.read_binary_grid(photHI_path, GridSize,
                                             photHI_precision, reshape=reshape)


    return nion_grid, XHII_grid, photHI_grid


def get_nion_prefix(SAGE_params):

    fescPrescription = int(SAGE_params["fescPrescription"])
    HaloPartCut = int(SAGE_params["HaloPartCut"])

    if fescPrescription == 0:
        beta = float(SAGE_params["beta"])
        nion_prefix = "fesc{0:.2f}_HaloPartCut{1}".format(beta,
                                                          HaloPartCut)

    elif fescPrescription == 1:
        alpha = float(SAGE_params["alpha"])
        beta = float(SAGE_params["beta"])

        nion_prefix = "ejected_{0:.3f}_{1:.3f}_HaloPartCut{2}".format(alpha,
                                                                      beta, 
                                                                      HaloPartCut)

    elif fescPrescription == 2:
        baseline = float(SAGE_params["quasar_baseline"])
        boosted = float(SAGE_params["quasar_boosted"])
        N_dyntime = float(SAGE_params["N_dyntime"])

        nion_prefix = "quasar_{0:.2f}_{1:.2f}_{2:.2f}_HaloPartCut{3}".format(baseline,
                                                                             boosted,
                                                                             N_dyntime,
                                                                             HaloPartCut)

    elif fescPrescription == 3 or fescPrescription == 4:
        MH_low = float(SAGE_params["MH_low"])
        MH_high = float(SAGE_params["MH_high"])
        fesc_low = float(SAGE_params["fesc_low"])
        fesc_high = float(SAGE_params["fesc_high"])

        nion_prefix = "AnneMH_{0:.3e}_{1:.2f}_{2:.3e}_{3:.2f}_HaloPartCut{4}".format(MH_low, fesc_low, MH_high,
                                              fesc_high)

    else:
        alpha = float(SAGE_params["alpha"])
        beta = float(SAGE_params["beta"])
        delta = float(SAGE_params["delta"])

        nion_prefix = "SFR_{0:.3f}_{1:.3f}_{2:.3f}_HaloPartCut{3}".format(alpha,
                                                                          beta,
                                                                          delta,
                                                                          HaloPartCut)

    return nion_prefix


def check_grids(nion_grid, XHII_grid, photHI_grid, SAGE_params, cifog_params,
                snapshot, update_data=0):

    print("")
    print("Now checking that the nionHI, XHII and photHI grids match the test "
          "data.")

    RunPrefix = SAGE_params["RunPrefix"]
    nion_precision = int(cifog_params["nionFilesAreInDoublePrecision"]) + 1
    GridSize = int(SAGE_params["GridSize"])

    tags = ["nionHI", "XHII", "photHI"]
    grids = [nion_grid, XHII_grid, photHI_grid]
    precisions = [nion_precision, 2, 2]  # XHII and photHI hard code as double.

    if update_data:
        print("=======================================================")
        print("WARNING WARNING WARNING WARNING WARNING WARNING WARNING")
        print("=======================================================")
        print("=======================================================")
        print("WARNING WARNING WARNING WARNING WARNING WARNING WARNING")
        print("=======================================================")
        input("YOU ARE ABOUT TO OVERWRITE THE GRID TEST DATA. IF THIS IS WHAT "
              "YOU WANT, PRESS ENTER OTHERWISE CTRL-C TO GET OUTTA HERE!")
        input("JUST CHECKING ONCE MORE!")

        for (grid, tag, precision) in zip(grids, tags, precions):

            fname = "{0}/{1}/kali512/data/{2}_test_{3}_{4:03d}"\
                    .format(test_dir, test_datadir, RunPrefix, tag, snapshot)

            # We are passed a 3D array, need to save it as 1D binary.
            np.savetxt(fname, grid)
            print("Saved {0} grid data as {1}".format(tag, fname))

        print("All grid test data updated. Exiting checking now (because it'll "
              "obviously be correct.)")
        return

    # Now let's compare the grids to the test data.

    for (grid, tag, precision) in zip(grids, tags, precisions):
        
        fname = "{0}/{1}/kali512/data/{2}_test_{3}_{4:03d}"\
                .format(test_dir, test_datadir, RunPrefix, tag, snapshot)

        test_grid = ReadScripts.read_binary_grid(fname, GridSize, precision, reshape=False)

        diff = np.abs(grid-test_grid)       
        if len(diff[diff > 1e-8]) > 0: 

            print("Found that the {0} grids disagreed.".format(tag))
            print("The SAGE dictionary is {0}".format(SAGE_params))
            print("The cifog dictionary is {0}".format(cifog_params))

            print("The mean value of the run grid is {0:.6e} compared to the mean "
                  "value of the test grid is {1:.6e}".format(np.mean(grid),
                                                             np.mean(test_grid)))
            print("The non-zero difference values are {0}".format(diff[diff > 0]))

            tol = 0.1 
            print("Checking if any of these values have a fractional difference "
                  "greater than {0}".format(tol))
            greater_0 = np.where(test_grid > 0.0)[0]
            fractional = diff[greater_0]/test_grid[greater_0]
            wrong_vals = greater_0[np.where(fractional > tol)[0]]
            if len(wrong_vals) > 0:
                print("Run values {0}".format(grid[wrong_vals]))
                print("Correct values {0}".format(test_grid[wrong_vals]))
                print("Diff {0}".format(diff[wrong_vals]))
                print("Fractional diff {0}".format(diff[wrong_vals]/test_grid[wrong_vals]))
                raise ValueError

    print("Grids are checked and all correct!")
    print("")

    return


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

    downloaded_repo = get_trees()  # Download Kali tree if we needed 

    # We have multiple test parameter specs we want to test.
    # Need all the names of the ini files and the name of the galaxies they
    # produce. 

    ini_files = ["kali512"]
   
    for ini_file in ini_files:

        # First read the ini file to get the runtime parameters.
        path_to_sage_ini = "{0}/test_ini_files/{1}_SAGE.ini".format(test_dir, ini_file)
        path_to_cifog_ini = "{0}/test_ini_files/{1}_cifog.ini".format(test_dir, ini_file)

        SAGE_params = ReadScripts.read_SAGE_ini(path_to_sage_ini)
        cifog_params, cifog_headers = ReadScripts.read_cifog_ini(path_to_cifog_ini) 

        run_prefix = SAGE_params["RunPrefix"]
        max_snap = int(SAGE_params["LastSnapShotNr"])

        # Make sure all the directories we need are present.
        # SAGE itself will take care of the directories for the results. 
        check_sage_dirs(run_prefix)

        # Then run SAGE.
        run_my_sage(ini_file)

        print("")
        print("SAGE run, now reading in the Galaxies.")
        print("")

        # Read the results and check the stellar mass function.
        Gals = load_gals(max_snap, run_prefix)
        check_smf(Gals, run_prefix, max_snap, SAGE_params)

        # Read the grids and check.
        snapshot = 61
        nion_grid, XHII_grid, photHI_grid = read_grids(SAGE_params,
                                                       cifog_params,
                                                       snapshot)
        check_grids(nion_grid, XHII_grid, photHI_grid, SAGE_params,
                    cifog_params, snapshot)
    
    print("Done")
    print("")
    
    #cleanup(downloaded_repo)

    print("")
    print("================")
    print("All tests passed")
    print("================")
    print("")


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

    # Remove all the output galaxy files.
   
    print("Done")
    print("")


if __name__ == "__main__":

    test_run()
