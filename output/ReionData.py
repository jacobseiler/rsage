#!/usr/bin/env python
"""
This file contains the functions for calculating the reionization data used to 
make plots with ``ReionPlots.py``.  It is called from ``paper_plots.py`` with a
dictionary containing which exact plots we require. 

You should not need to touch this file.  Please refer to the ``paper_plots.py``
documentation for full information on how to use this plotting pipeline. 

Author: Jacob Seiler
Version: 0.2
"""

from __future__ import print_function

import numpy as np
from numpy.fft import fftn, ifftn
import scipy.integrate as integrate
import os

import AllVars as av
import ReadScripts as rs
import PlotScripts as ps
import CollectiveStats as collective
import GalaxyData as gd 
import ReionPlots as reionplot

def T_naught(z, h, OM, OB):
    """
    Calculates the 21cm brightness temperature for specified cosmology +
    redshift.

    Parameters
    ---------
    z : Float
        Redshift we're calculating at.

    h : Float
        The Hubble h parameter defined as H = 100*h.

    OM : Float
        Critical matter density.

    OB : Float
        Critical Baryon density.

    Returns
    -------
    T0: Float
        21cm brightness temperature in units of mK.

    Units
    -----
    All input parameters are unitless.
    Output brightness temperature has units of mK.
    """

    T0 = 28.5 * ((1.0+z)/10.0)**(0.5) * OB/0.042 * h/0.73 * (OM/0.24)**(-0.5)
    return T0


def calc_ps(XHII, density, boxsize):
    """
    Calculates the 21cm and XHI power spectrum. 

    Parameters
    ---------
    XHII: 3-Dimensional Array of floats. Required. 
        Grid that contains the fraction of ionized hydrogen (XHII) in each 
        cell. 

    density: 3-Dimensional Array of floats. Required. 
        Grid that contains the overdensity (rho/<rho>) of dark matter in each
        cell. 

    Returns
    -------
    kmid_bins: 1-Dimensional array of floats. 
        The middle of each wavenumber bin. 

    powerspec: 1-Dimensional array of floats.
        The 21cm power spectrum in each k bin.

    p_err: 1-Dimensiona array of floats.
        The error on the 21cm power spectrum in each k bin. 


    The same again except for just the HII power spectrum.

    Units
    -----
    XHII is unitless.
    density is the overdensity, rho/<rho>.
    The wavenumber bins (kmid_bins) is in units of h/Mpc.
    The 21cm power spectrum (and associated error) is in units of Mpc^3/h^3/(2*pi). 
    """

    GridSize = np.shape(XHII)[0]

    meanion = np.mean(XHII)
    Tb = (1.0 - XHII)*density
    modes = ifftn(Tb)

    kmid_bins, powerspec, p_err = av.modes_to_pspec(modes,
                                                    boxsize=boxsize)

    kmid_bins_XHII, pspec_XHII, p_err_XHII = av.modes_to_pspec(ifftn(XHII),
                                                               boxsize=boxsize)

    return (kmid_bins, powerspec, p_err,
            kmid_bins_XHII, pspec_XHII, p_err_XHII)


def determine_ps_fixed_XHI(rank, size, comm,
                           z_array_reion_allmodels, cosmology_allmodels,
                           mass_frac_allmodels, XHII_fbase_allmodels,
                           XHII_precision_allmodels, density_fbase_allmodels,
                           density_precision_allmodels, GridSize_allmodels,
                           boxsize_allmodels, first_snap_allmodels,
                           fixed_XHI_values):

    num_models = len(mass_frac_allmodels)
    num_fractions = len(fixed_XHI_values)

    # Determine which mass_frac indices correspond to the requested XHI values.
    snap_idx_target = []
    for model_number in range(num_models):
        snap_idx_target.append([])
        for val in fixed_XHI_values:
            idx = (np.abs(mass_frac_allmodels[model_number] - val)).argmin()
            snap_idx_target[model_number].append(idx)

    flat_snap_idx_target = [item for sublist in snap_idx_target for item in sublist]

    # Now every rank knows what snapshot indices correspond to the fixed HI
    # values. We run 'classic' MPI loop. 
    k = []
    P21 = []
    PHII = []
    for idx in range(rank, num_models*num_fractions, size): 
        # Determine what model this corresponds to.
        model_number = int(idx / num_fractions)
        model_boxsize = boxsize_allmodels[model_number]
        model_cosmo = cosmology_allmodels[model_number]

        # The `snap_idx_target` value will be the relative snapshot number.
        # Hence need to add the `first_snap` for this model to get absolute.
        snapnum = flat_snap_idx_target[idx]
        redshift = z_array_reion_allmodels[model_number][snapnum]

        snapnum += first_snap_allmodels[model_number]
        cifog_snapnum = snapnum + 1

        # Load the XHII and density fields and calculate!
        XHII_path = "{0}_{1:03d}".format(XHII_fbase_allmodels[model_number],
                                         cifog_snapnum)
        XHII = rs.read_binary_grid(XHII_path, GridSize_allmodels[model_number],
                                   XHII_precision_allmodels[model_number]) 

        density_path = "{0}{1:03d}.dens.dat".format(density_fbase_allmodels[model_number],
                                                    snapnum)
        density = rs.read_binary_grid(density_path,
                                      GridSize_allmodels[model_number],
                                      density_precision_allmodels[model_number])


        T0 = T_naught(redshift, model_cosmo.H(0).value/100.0,
                      model_cosmo.Om0, model_cosmo.Ob0)

        tmp_k, tmp_PowSpec, tmp_Error, \
        tmp_k_XHII, tmp_Pspec_HII, tmp_Error_XHII = calc_ps(XHII, density,
                                                             model_boxsize)

        k.append(tmp_k)
        P21.append(tmp_PowSpec * T0*T0 * tmp_k**3 * 2.0*np.pi)
        PHII.append(tmp_Pspec_HII)

    comm.Barrier()

    # Now at this point each rank has a subset of the power spectra.
    # What we want to do is go through each index again and pass all of these
    # back onto the master rank.
    if rank == 0:
        k_master = []
        P21_master = []
        PHII_master = []

        rank_count = 0
        model_count = -1
    
        for idx in range(0, num_models*num_fractions):
            # This is the idx within each rank.
            ps_array_idx = int(idx / size)

            # If we've reached the end of the number of fractions, go to next
            # model. 
            if idx % num_fractions == 0: 
                model_count += 1
                k_master.append([])
                P21_master.append([])
                PHII_master.append([])

            # For every non-zero rank, we need to wait to receive the data from
            # the other processers.
            if rank_count == 0:
                k_master[model_count].append(k[ps_array_idx])
                P21_master[model_count].append(P21[ps_array_idx])
                PHII_master[model_count].append(PHII[ps_array_idx])
            else:
                tag = int(rank_count*100 + ps_array_idx)

                k_this_idx = comm.recv(source = rank_count,
                                       tag = tag) 
                P21_this_idx = comm.recv(source = rank_count,
                                         tag = tag+1) 
                PHII_this_idx = comm.recv(source = rank_count,
                                          tag = tag+2) 

                k_master[model_count].append(k_this_idx)
                P21_master[model_count].append(P21_this_idx)
                PHII_master[model_count].append(PHII_this_idx)

            rank_count += 1

            if rank_count == size:
                rank_count = 0

        return k_master, P21_master, PHII_master

    else:
        # We generate a unique tag for each rank + idx combination and then
        # send this processor's data to the master rank.
        for idx in range(len(P21)):
            tag = int(rank*100 + idx)

            k_this_idx = k[idx]
            P21_this_idx = P21[idx]
            PHII_this_idx = PHII[idx]

            comm.send(k_this_idx, dest = 0, tag = tag)
            comm.send(P21_this_idx, dest = 0, tag = tag+1)
            comm.send(PHII_this_idx, dest = 0, tag = tag+2)


        return None, None, None


def calc_tau(z_array_reion_allmodels, cosmology_allmodels, helium_allmodels,
             mass_frac_allmodels):

    def integrand(z, h, OM):
        H = av.Hubble_Param(z, h, OM) / (av.pc_to_m * 1.0e6 / 1.0e3)
        return (((1 + z)**2) / H)

    tau = []
    for model_number in range(len(mass_frac_allmodels)):

        # Set up some things for the model cosmology etc.
        model_mass_frac = mass_frac_allmodels[model_number]
        model_helium = helium_allmodels[model_number]
        model_h = cosmology_allmodels[model_number].H(0).value/100.0
        model_OM = cosmology_allmodels[model_number].Om0
        model_OB = cosmology_allmodels[model_number].Ob0
        model_z = z_array_reion_allmodels[model_number]

        model_tau = np.zeros(len(model_mass_frac))

        # First determine optical depth for redshift 0 to 4.
        tau_04 = integrate.quad(integrand, 0, 4, args=(model_h, model_OM,))[0] 
        tau_04 *= (1 + 2*model_helium/(4 * (1-model_helium)))

        # Then determine optical depth from z = 4 to lowest z of model.
        tau_46 = integrate.quad(integrand, 4, model_z[-1], args=(model_h, model_OM,))[0]
        tau_46 *= (1 + model_helium/(4* (1-model_helium)))

        tau_06 = tau_04 + tau_46

        model_tau[-1] = tau_06

        for snapnum in np.arange(len(model_mass_frac) - 2, -1, -1):

            this_z = model_z[snapnum]
            prev_z = model_z[snapnum + 1]

            # Hubble Parameter in Mpc/s/Mpc. 
            H = av.Hubble_Param(this_z, model_h, model_OM) / (av.pc_to_m * 1.0e6 / 1.0e3)
            numerator = ((1 + this_z) **2) *  (1.0 - model_mass_frac[snapnum])
         
            model_tau[snapnum] = model_tau[snapnum+1] + (( numerator / H) * (this_z - prev_z) * (1 + model_helium/(4 * (1-model_helium)))) 

        model_tau *= av.n_HI(0, model_h, model_OB, model_helium) * av.c_in_ms * av.Sigmat

        tau.append(model_tau)

    return tau


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

    # Check to see if the output directory exists.
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print("Made output directory {0}".format(output_dir))

    # First calculate all the properties and statistics we need.
    reion_data = generate_data(rank, size, comm, reion_ini_files,
                               gal_ini_files, reion_plots)

    # Set label sizes, colours etc.
    ps.Set_Params_Plot()

    # Gather all the fractions onto the master process.
    # This will be used for many different plots. 
    master_mass_frac = collective.collect_hist_across_tasks(rank, comm, 
                                                            reion_data["mass_frac_allmodels"]) 
    master_mass_frac = comm.bcast(master_mass_frac, root = 0)

    # Then find out what we need and plot em!
    if reion_plots["history"] and rank == 0:
        print("Plotting the reionization history.")
        reionplot.plot_history(reion_data["z_array_reion_allmodels"],
                               reion_data["lookback_array_reion_allmodels"],
                               reion_data["cosmology_allmodels"],
                               reion_data["t_bigbang_allmodels"],
                               master_mass_frac,
                               reion_plots["duration_definition"],
                               model_tags, output_dir, "history")

    if reion_plots["nion"]:

        master_nion = collective.collect_hist_across_tasks(rank, comm, 
                                                           reion_data["nion_allmodels"])

        if rank == 0:
            print("Plotting the ionizing emissivity.")
            reionplot.plot_nion(reion_data["z_array_reion_allmodels"],
                                reion_data["lookback_array_reion_allmodels"],
                                reion_data["cosmology_allmodels"],
                                reion_data["t_bigbang_allmodels"],
                                master_nion, 
                                model_tags, output_dir, "nion")

    if reion_plots["ps_fixed_XHI"]:
        print("Rank {0} Calculating power spectra at fixed neutral "
              "fractions.".format(rank))
        k, P21, PHII = determine_ps_fixed_XHI(rank, size, comm,
                                              reion_data["z_array_reion_allmodels"],
                                              reion_data["cosmology_allmodels"],
                                              master_mass_frac, 
                                              reion_data["XHII_fbase_allmodels"],
                                              reion_data["XHII_precision_allmodels"],
                                              reion_data["density_fbase_allmodels"],
                                              reion_data["density_precision_allmodels"],
                                              reion_data["GridSize_allmodels"],
                                              reion_data["boxsize_allmodels"],
                                              reion_data["first_snap_allmodels"],
                                              reion_plots["fixed_XHI_values"])

        if rank == 0:
            print("Plotting PS at fixed neutral fraction.")
            reionplot.plot_ps_fixed_XHI(k, P21, PHII,
                                        reion_plots["fixed_XHI_values"],
                                        model_tags, output_dir, "ps_fixed_XHI")

    if reion_plots["duration_contours"] and rank == 0:
        print("Plotting the duration of reionization contours.")
        reionplot.plot_duration_contours(reion_data["z_array_reion_allmodels"],
                                         reion_data["lookback_array_reion_allmodels"],        
                                         reion_data["cosmology_allmodels"],
                                         reion_data["t_bigbang_allmodels"],
                                         master_mass_frac,
                                         reion_plots["duration_contours_limits"],
                                         reion_plots["duration_definition"],
                                         model_tags, output_dir,
                                         "duration_contours")


    if reion_plots["optical_depth"] and rank == 0:
        tau_allmodels = calc_tau(reion_data["z_array_reion_allmodels"],
                                 reion_data["cosmology_allmodels"],
                                 reion_data["helium_allmodels"],
                                 master_mass_frac)

        print("Plotting the optical depth.")
        reionplot.plot_tau(reion_data["z_array_reion_allmodels"],
                           reion_data["lookback_array_reion_allmodels"],        
                           reion_data["cosmology_allmodels"],
                           reion_data["t_bigbang_allmodels"],
                           tau_allmodels,
                           model_tags, output_dir, "optical_depth")

    if reion_plots["optical_depth"] and reion_plots["history"] and rank == 0:
        reionplot.plot_combined_history_tau(reion_data["z_array_reion_allmodels"],
                                            reion_data["lookback_array_reion_allmodels"],    
                                            reion_data["cosmology_allmodels"],
                                            reion_data["t_bigbang_allmodels"],
                                            master_mass_frac, tau_allmodels, 
                                            model_tags, output_dir,
                                            "history_tau")

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

    if rank == 0:
        print("Generating reionization data for a total of {0} "
              "models.".format(len(reion_ini_files)))

    # ============================= #
    # General stuff for each model.
    # ============================= #

    # Unlike GalaxyData where we use the whole redshift range, here we only use
    # the range that covers reionization.
    z_array_reion_allmodels = []
    lookback_array_reion_allmodels = []

    # Since we calculate some things after the main loop, need to remember the
    # file names/precisions for each model.
    XHII_fbase_allmodels = []
    XHII_precision_allmodels = []
    density_fbase_allmodels = []
    density_precision_allmodels = []
    first_snap_allmodels = []

    GridSize_allmodels = []
    boxsize_allmodels = []
    helium_allmodels = []

    cosmology_allmodels = []
    t_bigbang_allmodels = []

    volume_frac_allmodels = []
    mass_frac_allmodels = []

    nion_allmodels = []

    # All outer arrays set up, time to read in the data!
    for model_number, (reion_ini_file, gal_ini_file) in \
            enumerate(zip(reion_ini_files, gal_ini_files)):

        if rank == 0:
            print("Model {0}".format(model_number))

        # Read in the parameters and set some initial variables.
        cifog_params, _ = rs.read_cifog_ini(reion_ini_file)
        SAGE_params = rs.read_SAGE_ini(gal_ini_file)

        cosmology, t_bigbang = gd.set_cosmology(float(SAGE_params["Hubble_h"]),
                                                float(SAGE_params["Omega"]),
                                                float(cifog_params["omega_b"]))

        cosmology_allmodels.append(cosmology)
        t_bigbang_allmodels.append(t_bigbang)

        first_snap = int(SAGE_params["LowSnap"])
        first_snap_allmodels.append(first_snap)

        last_snap = int(SAGE_params["LastSnapShotNr"])

        GridSize = int(SAGE_params["GridSize"])
        GridSize_allmodels.append(GridSize)

        boxsize = float(SAGE_params["BoxSize"])
        boxsize_allmodels.append(boxsize)

        helium = float(cifog_params["Y"])
        helium_allmodels.append(helium)

        model_volume = pow(float(SAGE_params["BoxSize"]) / \
                           float(SAGE_params["Hubble_h"]),3)
        model_hubble_h = float(SAGE_params["Hubble_h"])
        model_halopartcut = int(SAGE_params["HaloPartCut"])

        # Load the redshift file and calculate the lookback times. 
        z_array_full, lookback_array_full = gd.load_redshifts(SAGE_params["FileWithSnapList"],
                                                              cosmology, t_bigbang)
        z_array_reion = np.array(z_array_full[first_snap:last_snap])
        lookback_array_reion = np.array(lookback_array_full[first_snap:last_snap])

        z_array_reion_allmodels.append(z_array_reion)
        lookback_array_reion_allmodels.append(lookback_array_reion)

        XHII_fbase = cifog_params["output_XHII_file"]
        XHII_fbase_allmodels.append(XHII_fbase)

        density_fbase = cifog_params["inputIgmDensityFile"]
        density_fbase_allmodels.append(density_fbase)

        nion_fbase = cifog_params["inputNionFile"]

        # cifog uses 0 for floating point and 1 for double precision.
        # I use 0 for integer, 1 for floating point and 2 for double precision.
        density_precision = int(cifog_params["densityFilesAreInDoublePrecision"])
        density_precision += 1
        density_precision_allmodels.append(density_precision)

        nion_precision = int(cifog_params["nionFilesAreInDoublePrecision"])
        nion_precision += 1

        XHII_precision = 2  # XHII is assumed to have double precision.
        XHII_precision_allmodels.append(XHII_precision)

        # ============ #
        # Array setup. #
        # ============ #

        ## NOTE NOTE NOTE NOTE NOTE NOTE ##
        # These are all HI fractions # 
        volume_frac_allmodels.append(np.zeros(last_snap - first_snap))
        mass_frac_allmodels.append(np.zeros(last_snap - first_snap))

        nion_allmodels.append(np.zeros(last_snap - first_snap))

        # All arrays done, now loop over snapshots and read in.
        for snapnum in range(first_snap + rank, last_snap, size):

            # Where this snapshot slices into the global arrays.
            snap_idx = snapnum - first_snap

            # cifog numbering is weird and is shifted by +1.
            # E.g., nion file 027 is used to produce XHII file 028.
            cifog_snapnum = snapnum + 1

            XHII_path = "{0}_{1:03d}".format(XHII_fbase, cifog_snapnum)
            XHII = rs.read_binary_grid(XHII_path, GridSize, XHII_precision) 

            density_path = "{0}{1:03d}.dens.dat".format(density_fbase, snapnum)
            density = rs.read_binary_grid(density_path, GridSize,
                                          density_precision)

            # Be aware we track everything using the neutral HI fraction.
            # For the mass fraction, weight it by the density and normalize.
            volume_frac = 1.0 - np.mean(XHII)  # Nothing special for volume frac.
            mass_frac = 1.0 - np.sum(XHII * density / np.sum(density))

            volume_frac_allmodels[model_number][snap_idx] = volume_frac
            mass_frac_allmodels[model_number][snap_idx] = mass_frac

            if reion_plots["nion"]:
                nion_path = "{0}_{1:03d}".format(nion_fbase, snapnum)
                nion = rs.read_binary_grid(nion_path, GridSize, nion_precision)

                nion_allmodels[model_number][snap_idx] = np.sum(nion)

        # Snapshot Loop.

        # Ionizing emissitivty is scaled by the simulation volume (in Mpc^3).
        nion_allmodels[model_number] /= model_volume

    # Model Number Loop.

    reion_data = {"z_array_reion_allmodels" : z_array_reion_allmodels,
                  "lookback_array_reion_allmodels" : lookback_array_reion_allmodels,
                  "cosmology_allmodels" : cosmology_allmodels,
                  "t_bigbang_allmodels" : t_bigbang_allmodels,    
                  "volume_frac_allmodels" : volume_frac_allmodels,
                  "mass_frac_allmodels" : mass_frac_allmodels,
                  "nion_allmodels" : nion_allmodels,
                  "XHII_fbase_allmodels" : XHII_fbase_allmodels,
                  "XHII_precision_allmodels" : XHII_precision_allmodels,
                  "density_fbase_allmodels" : density_fbase_allmodels,
                  "density_precision_allmodels" : density_precision_allmodels,
                  "GridSize_allmodels" : GridSize_allmodels,
                  "boxsize_allmodels" : boxsize_allmodels,
                  "helium_allmodels" : helium_allmodels,
                  "first_snap_allmodels" : first_snap_allmodels}

    return reion_data
