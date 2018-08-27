#!/usr/bin/env python
"""
This file contains the functions for calculating the reionization data used to 
make plots with ``ReionPlots.py``.  It is called from ``paper_plots.py`` with a
dictionary containing which exact plots we require. 

You should not need to touch this file.  Please refer to the ``paper_plots.py``
documentation for full information on how to use this plotting pipeline. 

Author: Jacob Seiler
Version: 0.1
"""

from __future__ import print_function

import numpy as np
import os

import ReadScripts as rs
import PlotScripts as ps
import CollectiveStats as collective
import GalaxyData as gd 
import ReionPlots as reionplot

def plot_reion_properties(rank, size, comm, reion_ini_files, gal_ini_files,
                          model_tags, reion_plots, output_dir):
    """    
    Wrapper function to handle reading in of data + calculating reionization 
    properties, then calling the specified plotting routines.

    Parameters
    ----------

    reion_ini_files, gal_ini_files : List of strings 
        ``.ini`` file corresponding to each model that we're plotting.  We need
        both the galaxy (``SAGE``) and reionization (``cifog``) ``.ini`` files. 

    model_tags : List of strings
        String that will appear on the legend of the plot for each model.

    reion_plots : Dictionary
        Controls which of the plots we will make.  Keys are the name of each
        plot (e.g., ``history``) and the value specifies if we are plotting it.

    output_dir : String
        Directory where the plots are saved. If this directory does not exist,
        it is created beforehand. 

    Returns
    ---------

    None.
    """

    # First calculate all the properties and statistics we need.
    reion_data = generate_data(rank, size, comm, reion_ini_files,
                               gal_ini_files, reion_plots)

    # Check to see if the output directory exists.
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print("Made output directory {0}".format(output_dir))

    # Set label sizes, colours etc.
    ps.Set_Params_Plot()

    # Then find out what we need and plot em!
    if reion_plots["history"]:
        reionplot.plot_history(rank, comm,
                               reion_data["z_array_reion_allmodels"],
                               reion_data["lookback_array_reion_allmodels"],
                               reion_data["mass_frac_allmodels"],
                               reion_data["cosmology_allmodels"],
                               reion_data["t_bigbang_allmodels"],
                               model_tags, output_dir, "history")


def generate_data(rank, size, comm, reion_ini_files, gal_ini_files,
                  reion_plots):
    """    
    Reads in the galaxy data for calculate all the require properties for each
    models.

    Parameters
    ----------

    reion_ini_files, gal_ini_files : List of strings 
        ``.ini`` file corresponding to each model that we're plotting.  We need
        both the galaxy (``SAGE``) and reionization (``cifog``) ``.ini`` files. 

    reion_plots : Dictionary
        Controls which of the plots we will make.  Keys are the name of each
        plot (e.g., ``reion``) and the value specifies if we are plotting it. If
        we're not plotting a property we don't need to calculate stuff for it! 

    Returns
    ---------

    reion_data : Dictionary
        All of the calculated properties required to create the reionization
        plots.
    """

    # ============================= #
    # General stuff for each model.
    # ============================= #

    # Unlike GalaxyData where we use the whole redshift range, here we only use
    # the range that covers reionization.
    z_array_reion_allmodels = []
    lookback_array_reion_allmodels = []

    cosmology_allmodels = []
    t_bigbang_allmodels = []

    volume_frac_allmodels = []
    mass_frac_allmodels = []

    nion_allmodels = []

    # All outer arrays set up, time to read in the data!
    for model_number, (reion_ini_file, gal_ini_file) in \
            enumerate(zip(reion_ini_files, gal_ini_files)):

        # Read in the parameters and set some initial variables.
        cifog_params, _ = rs.read_cifog_ini(reion_ini_file)
        SAGE_params = rs.read_SAGE_ini(gal_ini_file)

        cosmology, t_bigbang = gd.set_cosmology(float(SAGE_params["Hubble_h"]),
                                                     float(SAGE_params["Omega"]))

        cosmology_allmodels.append(cosmology)
        t_bigbang_allmodels.append(t_bigbang)

        first_snap = int(SAGE_params["LowSnap"])
        last_snap = int(SAGE_params["LastSnapShotNr"])
        GridSize = int(SAGE_params["GridSize"])
        model_volume = pow(float(SAGE_params["BoxSize"]) / \
                           float(SAGE_params["Hubble_h"]),3)
        model_hubble_h = float(SAGE_params["Hubble_h"])
        model_halopartcut = int(SAGE_params["HaloPartCut"])

        # Load the redshift file and calculate the lookback times. 
        z_array_full, lookback_array_full = gd.load_redshifts(SAGE_params["FileWithSnapList"],
                                                              cosmology, t_bigbang)
        z_array_reion = np.array(z_array_full[first_snap:last_snap+1])
        lookback_array_reion = np.array(lookback_array_full[first_snap:last_snap+1])

        z_array_reion_allmodels.append(z_array_reion)
        lookback_array_reion_allmodels.append(lookback_array_reion)

        XHII_fbase = cifog_params["output_XHII_file"]
        density_fbase = cifog_params["inputIgmDensityFile"]
        nion_fbase = cifog_params["inputNionFile"]

        # cifog uses 0 for floating point and 1 for double precision.
        # I use 0 for integer, 1 for floating point and 2 for double precision.
        density_precision = int(cifog_params["densityFilesAreInDoublePrecision"])
        density_precision += 1

        nion_precision = int(cifog_params["nionFilesAreInDoublePrecision"])
        nion_precision += 1

        XHII_precision = 2  # XHII is assumed to have double precision.

        # ============ #
        # Array setup. #
        # ============ #

        ## NOTE NOTE NOTE NOTE NOTE NOTE ##
        # These are all HI fractions # 
        volume_frac_allmodels.append(np.zeros(last_snap - first_snap + 1))
        mass_frac_allmodels.append(np.zeros(last_snap - first_snap + 1))

        nion_allmodels.append(np.zeros(last_snap - first_snap + 1))

        # All arrays done, now loop over snapshots and read in.
        for count, snapnum in enumerate(range(first_snap + rank,
                                              last_snap + 1, size)):

            # cifog numbering is weird and is shifted by +1.
            # E.g., nion file 027 is used to produce XHII file 028.
            cifog_snapnum = snapnum + 1

            if reion_plots["history"]:
                XHII_path = "{0}_{1:03d}".format(XHII_fbase, cifog_snapnum)
                XHII = rs.read_binary_grid(XHII_path, GridSize, XHII_precision) 

                density_path = "{0}{1:03d}.dens.dat".format(density_fbase, snapnum)
                density = rs.read_binary_grid(density_path, GridSize,
                                              density_precision)

                # Be aware we track everything using the neutral HI fraction.
                # For the mass fraction, weight it by the density and normalize.
                volume_frac = 1.0 - np.mean(XHII)  # Nothing special for volume frac.
                mass_frac = 1.0 - np.sum(XHII * density / np.sum(density))

                volume_frac_allmodels[model_number][count] = volume_frac
                mass_frac_allmodels[model_number][count] = mass_frac

            if reion_plots["nion"]:
                nion_path = "{0}_{1:03d}".format(nion_fbase, snapnum)
                nion = rs.read_binary_grid(nion_path, GridSize, nion_precision)

                nion_allmodels[model_number][count] = np.sum(nion)

        # Snapshot Loop.
    # Model Number Loop.

    reion_data = {"z_array_reion_allmodels" : z_array_reion_allmodels,
                  "lookback_array_reion_allmodels" : lookback_array_reion_allmodels,
                  "cosmology_allmodels" : cosmology_allmodels,
                  "t_bigbang_allmodels" : t_bigbang_allmodels,    
                  "volume_frac_allmodels" : volume_frac_allmodels,
                  "mass_frac_allmodels" : mass_frac_allmodels,
                  "nion_allmodels" : nion_allmodels}

    return reion_data
