#!/usr/bin/env python
"""
This file contains the functions for calculating the galaxy data used to make
plots with ``GalaxyPlots.py``.  It is called from ``paper_plots.py`` with a
dictionary containing which exact plots we require. 

You should not need to touch this file.  Please refer to the ``paper_plots.py``
documentation for full information on how to use this plotting pipeline. 

Author: Jacob Seiler
Version: 0.2
"""

from __future__ import print_function

import numpy as np
import os

from astropy import units as u
from astropy import cosmology
from scipy import stats

import ReadScripts as rs
import PlotScripts as ps
import CollectiveStats as collective
import GalaxyPlots as galplot


def set_cosmology(Hubble_h, Omega_m, Omega_b):
    """    
    Sets a flat, LCDM cosmology.

    Parameters
    ----------

    Hubble_h : float
        Value of Hubble little h (i.e., between 0 and 1).


    Omega_m : float
        Value of critical matter density.

    Returns
    ---------

    cosmo : Class ``astropy.cosmology`` 
        ``Astropy`` class containing the cosmology for this model.

    t_bigbang : float
        The lookback time to the Big Bang in megayears. 
    """

    cosmo = cosmology.FlatLambdaCDM(H0 = Hubble_h*100, Om0 = Omega_m,
                                    Ob0 = Omega_b) 
    t_bigbang = (cosmo.lookback_time(100000).value)*1.0e3 # Lookback time to the Big Bang in Myr.

    return cosmo, t_bigbang


def load_redshifts(fname, cosmology, t_bigbang):
    """    
    Loads the redshifts for a specified model and calculates the corresponding
    lookback times.

    Parameters
    ----------

    fname : String 
        Path to the redshift file.

    cosmology : astropy.cosmology class
        ``Astropy`` class containing the cosmology for this model.  Set using
        ``set_cosmology()``.

    t_bigbang : float
        The lookback time to the Big Bang in megayears. Set using
        ``set_cosmology()``. 

    Returns
    ---------

    z : List of floats 
        Redshift at each snapshots of the specified model.

    t_bigbang : List of floats 
        Corresponding lookback time to each snapshot of the specified model.
    """

    with open(fname, "r") as f:

        a = np.loadtxt(f)

    z = 1.0/a - 1.0

    lookback = (t_bigbang - cosmology.lookback_time(z).value*1.0e3)  # In Myr.

    return z, lookback 


def do_2D_binning(data_x, data_y, mean_curr, std_curr, N_curr, bins):
    """    
    Updates the bin values (mean, standard deviation and number of data points)
    within 2D histograms.  That is, given x-data and y-data, updates y-data
    values binned on the x-data.

    Parameters
    ----------

    data_x, data_y : Numpy-arrays of floats 
        Data that we are using to update the histograms.        

    mean_curr, std_curr, N_curr : Numpy-arrays of floats
        Current mean, standard deviation and number of data points within each
        histogram bin.

    bins : Numpy-array of floats
        The bins we are binning the y-data on.  Defined in units/properties of
        the x-data.

    Returns
    ---------

    mean_curr, std_curr, N_curr : Numpy-arrays of floats
        The updated mean, standard deviation and number of data points within
        each histogram bin.
    """

    # First bin the y-data based on the binned x-data.
    snap_mean, _, _ = stats.binned_statistic(data_x, data_y, statistic='mean',
                                             bins=bins)
    snap_std, _, _ = stats.binned_statistic(data_x, data_y, statistic=np.std,
                                            bins=bins)
    snap_N, _, _ = stats.binned_statistic(data_x, data_y, statistic='count',
                                          bins=bins)

    # Prevent NaNs from messing things up.
    snap_mean[snap_N == 0] = 0.0
    snap_std[snap_N == 0] = 0.0

    # Update the running totals for these arrays. 
    mean_curr, std_curr = collective.update_cum_stats(mean_curr, std_curr,
                                                      N_curr, snap_mean,
                                                      snap_std, snap_N)
    N_curr += snap_N 

    return mean_curr, std_curr, N_curr


def plot_galaxy_properties(rank, size, comm, ini_files, model_tags, 
                           galaxy_plots, output_dir, output_format):
    """    
    Wrapper function to handle reading in of data + calculating galaxy
    properties, then calling the specified plotting routines.

    Parameters
    ----------

    rank : Integer
        This processor rank.

    size : Integer
        The total number of processors executing the pipeline.

    comm : Class ``mpi4py.MPI.Intracomm``
        The ``mpi4py`` communicator.

    ini_files : List of strings 
        ``.ini`` file corresponding to each model that we're plotting.

    model_tags : List of strings
        String that will appear on the legend of the plot for each model.

    galaxy_plots : Dictionary
        Controls which of the plots we will make.  Keys are the name of each
        plot (e.g., ``SMF``) and the value specifies if we are plotting it.

    output_dir : String
        Directory where the plots are saved. If this directory does not exist,
        it is created beforehand. 

    output_format : String
        The format of the saved figures.

    Returns
    ---------

    None. All figures are saved to the ``output_dir`` in format
    ``output_format``. 
    """

    # Check to see if the output directory exists.
    if rank == 0:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print("Made output directory {0}".format(output_dir))

    # First calculate all the properties and statistics we need.
    galaxy_data = generate_data(rank, size, comm, ini_files, galaxy_plots)

    # Set label sizes, colours etc.
    ps.Set_Params_Plot()

    # Then find what plots we need and plot em!
    if galaxy_plots["nion"]:
        galplot.plot_nion(rank, comm,
                          galaxy_data["z_array_full_allmodels"],
                          galaxy_data["lookback_array_full_allmodels"],
                          galaxy_data["sum_nion_allmodels"],
                          galaxy_data["cosmology_allmodels"],
                          galaxy_data["t_bigbang_allmodels"],
                          model_tags, output_dir, "nion", output_format)

    if galaxy_plots["mstar_fesc"]:
        galplot.plot_mstar_fesc(rank, comm,
                                galaxy_data["mstar_bins"],
                                galaxy_data["mstar_bin_width"],
                                galaxy_data["z_array_full_allmodels"],
                                galaxy_data["mean_mstar_fesc_allmodels"],
                                galaxy_data["std_mstar_fesc_allmodels"],
                                galaxy_data["N_mstar_fesc_allmodels"],
                                model_tags, output_dir, "mstar_fesc",
                                output_format,
                                plot_snaps_for_models=galaxy_plots["plot_snaps_for_models"],
                                plot_models_at_snaps=galaxy_plots["plot_models_at_snaps"])

    if galaxy_plots["SMF"]:
        galplot.plot_SMF(rank, comm,
                         galaxy_data["mstar_bins"],
                         galaxy_data["mstar_bin_width"],
                         galaxy_data["z_array_full_allmodels"],
                         galaxy_data["SMF_allmodels"],
                         galaxy_data["cosmology_allmodels"],
                         model_tags, output_dir, "SMF", output_format)

    if galaxy_plots["mstar_fej"]:
        galplot.plot_mstar_fej(rank, comm,
                               galaxy_data["mstar_bins"],
                               galaxy_data["mstar_bin_width"],
                               galaxy_data["z_array_full_allmodels"],
                               galaxy_data["mean_mstar_fej_allmodels"],
                               galaxy_data["std_mstar_fej_allmodels"],
                               galaxy_data["N_mstar_fej_allmodels"],
                               model_tags, output_dir, "mstar", output_format,
                               plot_snaps_for_models=galaxy_plots["plot_snaps_for_models"],
                               plot_models_at_snaps=galaxy_plots["plot_models_at_snaps"])

    if galaxy_plots["mstar_SFR"]:
        galplot.plot_mstar_SFR(rank, comm,
                               galaxy_data["mstar_bins"],
                               galaxy_data["mstar_bin_width"],
                               galaxy_data["z_array_full_allmodels"],
                               galaxy_data["mean_mstar_SFR_allmodels"],
                               galaxy_data["std_mstar_SFR_allmodels"],
                               galaxy_data["N_mstar_SFR_allmodels"],
                               model_tags, output_dir, "mstar_SFR",
                               output_format,
                               plot_snaps_for_models=galaxy_plots["plot_snaps_for_models"],
                               plot_models_at_snaps=galaxy_plots["plot_models_at_snaps"])


def generate_data(rank, size, comm, ini_files, galaxy_plots):
    """    
    Reads in the galaxy data for calculate all the require properties for each
    models.

    Parameters
    ----------

    rank : Integer
        This processor rank.

    size : Integer
        The total number of processors executing the pipeline.

    comm : Class ``mpi4py.MPI.Intracomm``
        The ``mpi4py`` communicator.

    ini_files : List of strings 
        ``.ini`` file corresponding to each model that we're plotting.

    galaxy_plots : Dictionary
        Controls which of the plots we will make.  Keys are the name of each
        plot (e.g., ``SMF``) and the value specifies if we are plotting it. If
        we're not plotting a property we don't need to calculate stuff for it! 

    Returns
    ---------

    galaxy_data : Dictionary
        All of the calculated properties required to create the plots.
    """

    # Binning parameters for stellar mass. 
    mstar_bin_low = 5.0
    mstar_bin_high = 12.0
    mstar_bin_width = 0.2
    mstar_Nbins = int((mstar_bin_high - mstar_bin_low) / mstar_bin_width)
    mstar_bins = np.arange(mstar_bin_low, 
                           mstar_bin_high + mstar_bin_width,
                           mstar_bin_width)

    # ======================================================================= #
    # We calculate values for all models and put them into lists that are     #
    # indexed by ``model_number``. So first we need to set up the outer-lists #
    # then we will append to these for each model.                            #
    # ======================================================================= #
    # General stuff for each model.
    z_array_full_allmodels = []
    lookback_array_full_allmodels = []

    z_array_reion_allmodels = []        
    lookback_array_reion_allmodels = []

    cosmology_allmodels = []
    t_bigbang_allmodels = []

    # These are the arrays for the number of ionizing photons at each snapshot.
    # Note: This is the ESCAPING ionizing photons. 
    sum_nion_allmodels = []

    # Escape fraction as a function of stellar mass (Mstar). 
    mean_mstar_fesc_allmodels = []
    std_mstar_fesc_allmodels = []
    N_mstar_fesc_allmodels = []

    # Stellar mass function.
    SMF_allmodels = []

    # Ejected fraction as a function of stellar mass (Mstar). 
    mean_mstar_fej_allmodels = []
    std_mstar_fej_allmodels = []
    N_mstar_fej_allmodels = []

    # Star formation rate as a function of stellar mass (Mstar). 
    mean_mstar_SFR_allmodels = []
    std_mstar_SFR_allmodels = []
    N_mstar_SFR_allmodels = []

    # All outer arrays set up, time to read in the data!
    for model_number, ini_file in enumerate(ini_files):

        # Read in the parameters and set some initial variables.
        SAGE_params = rs.read_SAGE_ini(ini_file)

        cosmology, t_bigbang = set_cosmology(float(SAGE_params["Hubble_h"]),
                                             float(SAGE_params["Omega"]),
                                             float(SAGE_params["BaryonFrac"]))

        cosmology_allmodels.append(cosmology)
        t_bigbang_allmodels.append(t_bigbang)

        first_snap = int(SAGE_params["LowSnap"])
        last_snap = int(SAGE_params["LastSnapShotNr"])
        GridSize = int(SAGE_params["GridSize"])
        model_hubble_h = float(SAGE_params["Hubble_h"])
        model_halopartcut = int(SAGE_params["HaloPartCut"])

        # Careful, volume is in Mpc^3.
        model_volume = pow(float(SAGE_params["BoxSize"]) / \
                           float(SAGE_params["Hubble_h"]),3)        

        # Load the redshift file and calculate the lookback times. 
        z_array_full, lookback_array_full = load_redshifts(SAGE_params["FileWithSnapList"],
                                                           cosmology, t_bigbang)
        z_array_full = np.array(z_array_full[0:last_snap+1])
        lookback_array_full = np.array(lookback_array_full[0:last_snap+1])

        z_array_full_allmodels.append(z_array_full)
        lookback_array_full_allmodels.append(lookback_array_full)

        # Specifically set the redshift range to cover reionization.
        z_array_reion = np.array(z_array_full[first_snap:last_snap+1])
        z_array_reion_allmodels.append(z_array_reion)
        lookback_array_reion = np.array(lookback_array_full[first_snap:last_snap+1])
        lookback_array_reion_allmodels.append(lookback_array_reion)

        # Set up names for the galaxies. 
        galaxy_name = "{0}/{1}_z{2:.3f}".format(SAGE_params["OutputDir"],
                                                SAGE_params["FileNameGalaxies"],
                                                z_array_reion[-1]) 

        merged_name = "{0}/{1}_MergedGalaxies".format(SAGE_params["OutputDir"],
                                                      SAGE_params["FileNameGalaxies"])

        # Initialize the ionizing photon array to 0.
        sum_nion_allmodels.append(np.zeros(len(z_array_full),
                                           dtype=np.float32))

        # Escape fraction as a function of stellar mass.
        mean_mstar_fesc_allmodels.append([]) 
        std_mstar_fesc_allmodels.append([])
        N_mstar_fesc_allmodels.append([]) 
        for snap_count in range(len(z_array_full)):
            mean_mstar_fesc_allmodels[model_number].append(np.zeros(mstar_Nbins,
                                                                    dtype=np.float32))
            std_mstar_fesc_allmodels[model_number].append(np.zeros(mstar_Nbins,
                                                                   dtype=np.float32))
            N_mstar_fesc_allmodels[model_number].append(np.zeros(mstar_Nbins,
                                                                 dtype=np.float32))

        # Stellar mass function. 
        SMF_allmodels.append([])
        for snap_count in range(len(z_array_full)):
            SMF_allmodels[model_number].append(np.zeros(mstar_Nbins,
                                               dtype=np.float32))

        # Ejected fraction as a function of stellar mass.
        mean_mstar_fej_allmodels.append([]) 
        std_mstar_fej_allmodels.append([])
        N_mstar_fej_allmodels.append([]) 
        for snap_count in range(len(z_array_full)):
            mean_mstar_fej_allmodels[model_number].append(np.zeros(mstar_Nbins,
                                                                   dtype=np.float32))
            std_mstar_fej_allmodels[model_number].append(np.zeros(mstar_Nbins,
                                                                  dtype=np.float32))
            N_mstar_fej_allmodels[model_number].append(np.zeros(mstar_Nbins,
                                                                dtype=np.float32))

        # Star formation rate as a function of stellar mass.
        mean_mstar_SFR_allmodels.append([]) 
        std_mstar_SFR_allmodels.append([])
        N_mstar_SFR_allmodels.append([]) 
        for snap_count in range(len(z_array_full)):
            mean_mstar_SFR_allmodels[model_number].append(np.zeros(mstar_Nbins,
                                                                   dtype=np.float32))
            std_mstar_SFR_allmodels[model_number].append(np.zeros(mstar_Nbins,
                                                                  dtype=np.float32))
            N_mstar_SFR_allmodels[model_number].append(np.zeros(mstar_Nbins,
                                                                dtype=np.float32))

        # Check to see if we're only using a subset of the files.
        if galaxy_plots["FirstFile"]:
            first_file = galaxy_plots["FirstFile"]
        else:
            first_file = int(SAGE_params["FirstFile"])

        if galaxy_plots["lastFile"]:
            last_file = galaxy_plots["FirstFile"]
        else:
            last_file = int(SAGE_params["LastFile"])

        # ========================================================= #
        # Now go through each file and calculate the stuff we need. #
        # ========================================================= #
        # Parallelize over number of files.
        for fnr in range(first_file + rank, last_file + 1, size):

            print("Rank {0}: Model {1} File {2}".format(rank, model_number,
                                                        fnr))

            # Read in both the galaxies, the merged ones and combine them into
            # a single array.
            GG, Gal_Desc = rs.ReadGals_SAGE(galaxy_name, fnr, len(z_array_full))
            G_Merged, _ = rs.ReadGals_SAGE(merged_name, fnr, len(z_array_full))
            G = rs.Join_Arrays(GG, G_Merged, Gal_Desc)

            # For each snapshot, calculate properties for galaxies that exist.
            for snap_count, snapnum in enumerate(range(len(z_array_full))):

                Gals_exist = np.where((G.GridHistory[:, snapnum] != -1) &
                                      (G.GridStellarMass[:, snapnum] > 0.0) &
                                      (G.LenHistory[:, snapnum] > model_halopartcut))[0]
                if len(Gals_exist) == 0:
                    continue

                sum_nion_allmodels[model_number][snap_count] += sum(G.GridNgamma_HI[Gals_exist, snapnum] * \
                                                                    G.Gridfesc[Gals_exist,snapnum])

                log_mass = np.log10(G.GridStellarMass[Gals_exist, snapnum] * 1.0e10 / model_hubble_h)
                fesc = G.Gridfesc[Gals_exist, snapnum]
                fej = G.EjectedFraction[Gals_exist, snapnum]
                SFR = G.GridSFR[Gals_exist, snapnum]

                # Calculate the mean fesc as a function of stellar mass.
                if galaxy_plots["mstar_fesc"]:
                    mean_mstar_fesc_allmodels[model_number][snap_count], \
                    std_mstar_fesc_allmodels[model_number][snap_count], \
                    N_mstar_fesc_allmodels[model_number][snap_count] = \
                        do_2D_binning(log_mass, fesc,
                                      mean_mstar_fesc_allmodels[model_number][snap_count],
                                      std_mstar_fesc_allmodels[model_number][snap_count],
                                      N_mstar_fesc_allmodels[model_number][snap_count],
                                      mstar_bins)

                # Calculate the mean ejected fraction as a function of stellar mass.
                if galaxy_plots["mstar_fej"]:
                    mean_mstar_fej_allmodels[model_number][snap_count], \
                    std_mstar_fej_allmodels[model_number][snap_count], \
                    N_mstar_fej_allmodels[model_number][snap_count] = \
                        do_2D_binning(log_mass, fej,
                                      mean_mstar_fej_allmodels[model_number][snap_count],
                                      std_mstar_fej_allmodels[model_number][snap_count],
                                      N_mstar_fej_allmodels[model_number][snap_count],
                                      mstar_bins)

                if galaxy_plots["mstar_SFR"]:
                    mean_mstar_SFR_allmodels[model_number][snap_count], \
                    std_mstar_SFR_allmodels[model_number][snap_count], \
                    N_mstar_SFR_allmodels[model_number][snap_count] = \
                        do_2D_binning(log_mass, SFR,
                                      mean_mstar_SFR_allmodels[model_number][snap_count],
                                      std_mstar_SFR_allmodels[model_number][snap_count],
                                      N_mstar_SFR_allmodels[model_number][snap_count],
                                      mstar_bins)
 
                SMF_thissnap = np.histogram(log_mass, bins=mstar_bins)
                SMF_allmodels[model_number][snap_count] += SMF_thissnap[0]
                # Snapshot loop.
            # File Loop.
        # Model Loop.

        # Ionizing emissitivty is scaled by the simulation volume (in Mpc^3).
        sum_nion_allmodels[model_number] /= model_volume
        
        # Stellar Mass Function is normalized by boxsize and bin width.
        SMF_allmodels[model_number] = np.divide(SMF_allmodels[model_number],
                                                model_volume * mstar_bin_width)

    # Everything has been calculated. Now construct a dictionary that contains
    # all the data (for easy passing) and return it. 
    galaxy_data = {"z_array_full_allmodels" : z_array_full_allmodels,
                   "lookback_array_full_allmodels" : lookback_array_full_allmodels,
                   "z_array_reion_allmodels" : z_array_reion_allmodels,
                   "lookback_array_reion_allmodels" : lookback_array_reion_allmodels,
                   "cosmology_allmodels" : cosmology_allmodels,
                   "t_bigbang_allmodels" : t_bigbang_allmodels,
                   "sum_nion_allmodels" : sum_nion_allmodels,
                   "mstar_bins" : mstar_bins, "mstar_bin_width" : mstar_bin_width,
                   "mean_mstar_fesc_allmodels" : mean_mstar_fesc_allmodels,
                   "std_mstar_fesc_allmodels" : std_mstar_fesc_allmodels,
                   "N_mstar_fesc_allmodels" : N_mstar_fesc_allmodels,
                   "SMF_allmodels" : SMF_allmodels,
                   "mean_mstar_fej_allmodels" : mean_mstar_fej_allmodels,
                   "std_mstar_fej_allmodels" : std_mstar_fej_allmodels,
                   "N_mstar_fej_allmodels" : N_mstar_fej_allmodels,
                   "mean_mstar_SFR_allmodels" : mean_mstar_SFR_allmodels,
                   "std_mstar_SFR_allmodels" : std_mstar_SFR_allmodels,
                   "N_mstar_SFR_allmodels" : N_mstar_SFR_allmodels}
    return galaxy_data
