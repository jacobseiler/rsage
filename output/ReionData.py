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
import time
import random

import AllVars as av
import ReadScripts as rs
import PlotScripts as ps
import CollectiveStats as collective
import GalaxyData as gd 
import ReionPlots as reionplot


def calc_duration(z_array_reion_allmodels, lookback_array_reion_allmodels,
                  mass_frac_allmodels, duration_definition):
    """
    Determines the duration of reionization.

    Parameters
    ----------

    z_array_reion_allmodels : 2D nested list of floats. Outer length is number
                              of models, inner is number of snapshots.
        The redshift at each snapshot for each model.

    lookback_array_reion_allmodels : 2D nested list of floats. Dimensions
                                     identical to ``z_array_reion_allmodels``. 
        The lookback time at each snapshot for each model. Units are Myr.

    mass_frac_allmodels : 2D nested list of floats. Dimensions equal to
                          ``z_array_reion_allmodels``.
        The mass weighted neutral fraction at each snapshot for each model.

    duration_definition : List of floats with length 3.
        The neutral fractions that define reionization.  The first element is
        the start, second is the mid-point and third is the end of
        reionization.

    Returns
    ---------

    duration_z, duration_t : 2D nested list of floats. Outer length is number
                             of models, inner is 3.
        The redshift and lookback time corresponding to each element of the
        ``duration_definition`` list.

    reion_completed : List of integers.
        Flag to denote whether reionization has been completedd by the final
        snapshot.
    """

    duration_z = []
    duration_t = []
    reion_completed = []

    # We need to be careful here. For low values of fesc, reionization
    # won't actually complete.  Hence we need to check `duration_z` and see
    # those models in which reionization is 'completed' at the last snapshot.
    for model_number in range(len(mass_frac_allmodels)):
        mass_frac_thismodel = mass_frac_allmodels[model_number]

        duration_z.append([]) 
        duration_t.append([])
        for val in duration_definition:
            idx = (np.abs(mass_frac_thismodel - val)).argmin()
            duration_z[model_number].append(z_array_reion_allmodels[model_number][idx]) 
            duration_t[model_number].append(lookback_array_reion_allmodels[model_number][idx]) 
            if (val == duration_definition[-1]) and \
               (idx == len(mass_frac_thismodel)-1):
                reion_completed.append(0)
            elif(val == duration_definition[-1]):
                reion_completed.append(1)

    return duration_z, duration_t, reion_completed


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

    T0 = 28.5 * ((1.0+z)/10.0)**(0.5) * OB/0.042 * h/0.73 * (0.24/OM)**(0.5)
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
    """    
    Calculates the 21cm and HII power spectra at fixed HI fractions.

    Parameters
    ----------

    rank : Integer
        This processor rank.

    size : Integer
        The total number of processors executing the pipeline.

    comm : Class ``mpi4py.MPI.Intracomm``
        The ``mpi4py`` communicator.

    z_array_reion_allmodels : 2D nested list of floats. Outer length is number
                              of models, inner is number of snapshots.
        The redshift at each snapshot for each model.

    cosmology_allmodels : List of class ``astropy.cosmology``. Length is number
                          of models.
        ``Astropy`` class containing the cosmology for each model.

    mass_frac_allmodels : 2D nested list of floats. Dimensions equal to
                          ``z_array_reion_allmodels``.
        The mass weighted neutral fraction at each snapshot for each model.

    XHII_fbase_allmodels : List of strings. Length is number of models.
        The base filename for the ionization fields for each model.

    XHII_precision_allmodels : List of integers. Length is number of models.
        The precision of the ionization fields for each model.
        1 : Float, 2 : Double.

    density_fbase_allmodels : List of strings. Length is number of models.
        The base filename for the density fields for each model.

    density_precision_allmodels : List of integers. Length is number of models.
        The precision of the density fields for each model.
        1 : Float, 2 : Double.

    GridSize_allmodels : List of integers. Length is number of models.
        The number of grid cells (along a box size) for each model.

    boxsize_allmodels : List of integers. Length is number of models.
        The simulation box size for each model (units are Mpc/h).

    first_snap_allmodels : List of integers. Length is number of models.
        The snapshot where ``cifog`` starts calculations for each model.

    fixed_XHI_values : List of floats.
        The neutral hydrogen fractions we're calculating the power spectra at.
        Defined by the user in ``paper_plots.py``.

    Returns
    ---------

    For all non-zero ranks, the returns are:
        None, None, None.

    k_master, P21_master, PHII_master : 3D nested lists of floats. Outer length
                                        is number of models, next is number of
                                        snapshots in the model and final is the
                                        number of wavenumber bins.
        The wavenumber bins, 21cm and HII power spectra for each model for the
        neutral hydrogen fractions `fixed_XHI_values`.
    """

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
        model_cosmo = cosmology_allmodels[model_number]
        model_boxsize = boxsize_allmodels[model_number] # Mpc/h.

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

        # Be careful, passing the boxsize at Mpc/h.
        tmp_k, tmp_PowSpec, tmp_Error, \
        tmp_k_XHII, tmp_Pspec_HII, tmp_Error_XHII = calc_ps(XHII, density,
                                                             model_boxsize)

        k.append(tmp_k)
        P21.append(tmp_PowSpec * T0*T0 * tmp_k**3 * 4.0*np.pi)
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
    """
    Calculates the Thomson integrated optical depth.

    Parameters
    ----------

    z_array_reion_allmodels : 2D nested list of floats. Outer length is number
                              of models, inner is number of snapshots.
        The redshift at each snapshot for each model.

    cosmology_allmodels : List of class ``astropy.cosmology``. Length is number
                          of models.
        ``Astropy`` class containing the cosmology for each model.

    helium_allmodels : Nested list of floats. Length is number of models,
        The helium fraction for each model.

    mass_frac_allmodels : 2D nested list of floats. Dimensions equal to
                          ``z_array_reion_allmodels``.
        The mass weighted neutral fraction at each snapshot for each model.

    Returns
    ---------

    tau : 2D nested list of floats. Dimensions equal to
          ``z_array_reion_allmodels``.
        The Thomson optical depth at each snapshot for each model.
    """

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

        # Then loop down through snapshots (low z to high z) and calculate tau.
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


def gather_ps(rank, size, comm, k_allmodels, P21_allmodels, PHII_allmodels,
              first_snap_allmodels, last_snap_allmodels):
    """
    Gathers the power spectra calculated on each processor onto the root rank.
    Each rank calculates the spectra of only a subset of snapshots so here we
    gather the spectra of all models and snapshots onto root rank. 

    Parameters
    ----------

    rank : Integer
        This processor rank.

    size : Integer
        The total number of processors executing the pipeline.

    comm : Class ``mpi4py.MPI.Intracomm``
        The ``mpi4py`` communicator.

    k_allmodels : 3D nested list of floats. Outer length is number of models,
                  next is number of snapshots processed by this processor and
                  final is the number of wavenumber bins.
        Wavenumber the spectra are binned on (units of Mpc/h). Processor
        unique.

    P21_allmodels, PHII_allmodels : 3D nested lists of floats. Dimensions are
                                    identical to ``k_allmodels``.
        The 21cm and HII power spectra for each model at each snapshot.
        Processor unique.

    first_snap_allmodels, last_snap_allmodels : List of integers. Length is
                                                number of models.
        The first and last snapshot that defines the snapshot range that
        ``cifog`` was run on. 

    Returns
    ---------

    For all non-zero ranks, the returns are:
        None, None, None.

    k_master, P21_master, PHII_master : 3D nested lists of floats. Dimensions
                                        are identical to `k_allmodels` except
                                        the 2nd axis length is the snapshot
                                        range for that particular model (not a
                                        subset).
        The wavenumber bins, 21cm and HII power spectra for each model for all
        snapshots. 
    """

    def generate_tag(rank):
        tag = int(rank*100)

        return tag

    # Rank 0 will gather the wavenumber bins/power spectra from all other
    # ranks. 
    if rank == 0:
        k_master = []
        P21_master = []
        PHII_master = []

        # Go through each model. 
        for model_number in range(len(k_allmodels)):

            k_master.append([])
            P21_master.append([])
            PHII_master.append([])

            model_k = k_allmodels[model_number]
            model_P21 = P21_allmodels[model_number]
            model_PHII = PHII_allmodels[model_number]

            num_snaps = last_snap_allmodels[model_number] - \
                        first_snap_allmodels[model_number]
            rank_count = 0
            my_count = 0

            # Then go through each snapshot.
            # In the main data loop (``generate_data()``) the snapshots are
            # scatter sequentially. Hence when we gather, we get snap0 from
            # rank 0, snap1 from rank 1 etc. So we increase rank_count for each
            # snapshot and then reset it when we reach `size`.
            for snap_idx in range(num_snaps):

                if rank_count == 0:
                    this_k = model_k[my_count] 
                    this_P21 = model_P21[my_count] 
                    this_PHII = model_PHII[my_count] 
                    my_count += 1
                else:
                    # Each rank will use a unique tag.
                    tag = generate_tag(rank_count) 

                    # Then the tag is offset for each data array. 
                    this_k = comm.recv(source = rank_count,
                                       tag = tag)
                    this_P21 = comm.recv(source = rank_count,
                                         tag = tag+1)
                    this_PHII = comm.recv(source = rank_count,
                                          tag = tag+2)

                # Now we have the data, append it to the master.
                k_master[model_number].append(this_k)
                P21_master[model_number].append(this_P21)
                PHII_master[model_number].append(this_PHII)

                rank_count += 1
                if rank_count == size:
                    rank_count = 0

            # Snapshot Loop.
        # Model Loop.

        return k_master, P21_master, PHII_master

    else:

        # For all other ranks, go through the power spectra it calculated and
        # send it back to the root rank.
        for model_number in range(len(k_allmodels)):
            for idx in range(len(P21_allmodels[model_number])):

                tag = generate_tag(rank) 

                k_this_idx = k_allmodels[model_number][idx]
                P21_this_idx = P21_allmodels[model_number][idx]
                PHII_this_idx = PHII_allmodels[model_number][idx]

                comm.send(k_this_idx, dest = 0, tag = tag)
                comm.send(P21_this_idx, dest = 0, tag = tag+1)
                comm.send(PHII_this_idx, dest = 0, tag = tag+2)

        # Non-zero ranks return junk.
        return None, None, None


def calculate_bubble_MC(ionized_cells, Ncell, output_file):
# Calculate the size of ionized/neutral bubbles by choosing an ionized/neutral cell and measuring the distance to a cell of the opposite phase.

## Input ##
# z is the redshift we are doing the MC walk at.
# ionized_cells is the array that contains the ionization state of a cell.

## Output ##
# The output file will be of the form '*OutputDir*/MC/*output_tag*_*z*.dat'

    def MC_walk(cells, indices, phase, N, Ncell, output_file):
        # Function for the MC walk to calculate bubble size.
    
        ## Input ##
        # cells is the array that contains the ionization field (0 or 1).
        # indices is a length 3 array containing the x,y,z indices of the cells that are neutral (phase = 0) or ionized (phase = 1)
        # phase is 0 if we are starting with neutral or 1 if we are starting with ionized cell.
        # N is the number of random points we select.
    
        ## Output ##
        # A .dat file containing the results of the MC walk (number of steps until it encounters a phase transition).

        radii = []

        for j in range(0,int(N)):
                        
            #if (j%20000 == 0):
            #    print("On walk number {0}".format(j))
    
			## Select a random direction to walk through ##
            direction = random.randint(1,6)

            if direction == 1:
                x = 1
                y = 0
                z = 0
            elif direction == 2:
                x = -1 
                y = 0
                z = 0
            elif direction == 3:
                x = 0
                y = 1
                z = 0
            elif direction == 4:
                x = 0
                y = -1
                z = 0
            elif direction == 5:
                x = 0
                y = 0
                z = 1
            else:
                x = 0
                y = 0
                z = -1

			# Pick the x,y,z coordinates of a random ionized/neutral cell
            random_index = random.randint(0,len(indices[0])-1)
            walk_x = indices[0][random_index]
            walk_y = indices[1][random_index]
            walk_z = indices[2][random_index]

            R = 0
            phase_transition = phase

            while (phase_transition == phase): # While we haven't changed phase yet.
                R += 1 # Keep increasing the radius until we do.
                phase_transition = cells[(walk_x + R*x) % Ncell, (walk_y + R*y) % Ncell, (walk_z + R*z) % Ncell]
                if (phase_transition > 0.8):
                    phase_transition = 1
                else:
                    phase_transition = 0	
                if (R >= Ncell): # If the radius has gone beyond the number of available cells, 
                    phase_transition = (phase + 1) % 2 # We force the change.

            radii.append(R)
			
        np.savetxt(output_file, radii, delimiter = ',')
        #print("MC file saved as {0}".format(outfile))

    #print("")
    #print("Calculating bubble size using MC walk.")	

    N = 1e5

    start_time = time.time()

    ionized_indices= np.where(ionized_cells > 0.8) # Calculate the array indices corresponding to neutral/ionized cells.
    unionized_indices = np.where(ionized_cells < 0.2)

    MC_walk(ionized_cells, ionized_indices, 1, N, Ncell, output_file)

    #print("MC took {0} seconds.".format(time.time() - start_time))


def plot_reion_properties(rank, size, comm, reion_ini_files, gal_ini_files,
                          model_tags, reion_plots, output_dir, output_format):
    """    
    Wrapper function to handle reading in of data + calculating reionization 
    properties, then calling the specified plotting routines.

    Parameters
    ----------

    rank : Integer
        This processor rank.

    size : Integer
        The total number of processors executing the pipeline.

    comm : Class ``mpi4py.MPI.Intracomm``
        The ``mpi4py`` communicator.

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
    reion_data = generate_data(rank, size, comm, reion_ini_files,
                               gal_ini_files, reion_plots, output_dir,
                               model_tags)

    # Gather all the fractions onto the master process.
    # This will be used for many different plots. 
    master_mass_frac = collective.collect_hist_across_tasks(rank, comm, 
                                                            reion_data["mass_frac_allmodels"]) 
    master_mass_frac = comm.bcast(master_mass_frac, root = 0)

    # Then find out what we need and plot em!
    if reion_plots["history"] and rank == 0:

        
        duration_z, duration_t, reion_completed = \
            calc_duration(reion_data["z_array_reion_allmodels"],
                          reion_data["lookback_array_reion_allmodels"],
                          master_mass_frac, reion_plots["duration_definition"])

        for model_number in range(len(master_mass_frac)):
            print("Model {0}: Start {1:.3f} \tMid {2:.3f}\tEnd {3:.23}\t"
                  "dz {4:.3f}\tdt {5:.2f}Myr\tReion Completed {6}" \
                  .format(model_number, duration_z[model_number][0],
                          duration_z[model_number][1], duration_z[model_number][-1],
                          duration_z[model_number][0]-duration_z[model_number][-1],
                          duration_t[model_number][-1]-duration_t[model_number][0],
                          reion_completed[model_number]))

        print("Plotting the reionization history.")
        reionplot.plot_history(reion_data["z_array_reion_allmodels"],
                               reion_data["lookback_array_reion_allmodels"],
                               reion_data["cosmology_allmodels"],
                               reion_data["t_bigbang_allmodels"],
                               master_mass_frac,
                               model_tags, output_dir, "history",
                               output_format)



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
                                reion_data["nion_factor_allmodels"], 
                                model_tags, output_dir, "nion", output_format)

    if reion_plots["ps_fixed_XHI"]:
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
                                        model_tags, output_dir, "ps_fixed_XHI",
                                        output_format)

    if reion_plots["contours"] and rank == 0:
        # tau is used for multiple plots. So check if we need to calculate it.
        try:
            tau_allmodels
        except NameError:
            tau_allmodels = calc_tau(reion_data["z_array_reion_allmodels"],
                                     reion_data["cosmology_allmodels"],
                                     reion_data["helium_allmodels"],
                                     master_mass_frac)

        # For the contours, only plot the optical depth at the highest z.
        tau_highz = []
        for model_number in range(len(tau_allmodels)):
            tau_highz.append(tau_allmodels[model_number][0])

        duration_z, duration_t, reion_completed = \
            calc_duration(reion_data["z_array_reion_allmodels"],
                          reion_data["lookback_array_reion_allmodels"],
                          master_mass_frac, reion_plots["duration_definition"])

        print("Plotting contours of constant tau.")
        reionplot.plot_tau_contours(tau_highz, reion_completed,
                                    reion_plots["alpha_beta_limits"],
                                    output_dir, "tau_contours", output_format)

        print("Plotting contours of constant reionization duration.")
        reionplot.plot_duration_contours(duration_t, reion_completed,
                                         reion_plots["alpha_beta_limits"],
                                         output_dir, "duration_contours",
                                         output_format)

    if reion_plots["optical_depth"] and rank == 0:
        # tau is used for multiple plots. So check if we need to calculate it.
        try:
            tau_allmodels
        except NameError:
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
                           model_tags, output_dir, "optical_depth",
                           output_format)

    if reion_plots["optical_depth"] and reion_plots["history"] and rank == 0:
        reionplot.plot_combined_history_tau(reion_data["z_array_reion_allmodels"],
                                            reion_data["lookback_array_reion_allmodels"],    
                                            reion_data["cosmology_allmodels"],
                                            reion_data["t_bigbang_allmodels"],
                                            master_mass_frac, tau_allmodels, 
                                            model_tags, output_dir,
                                            "history_tau", output_format)

    if reion_plots["ps_scales"]:
        k, P21, PHII = gather_ps(rank, size, comm,
                                 reion_data["k_allmodels"],
                                 reion_data["P21_allmodels"],
                                 reion_data["PHII_allmodels"],
                                 reion_data["first_snap_allmodels"],
                                 reion_data["last_snap_allmodels"])

        if rank == 0:
            reionplot.plot_ps_scales(k, P21, PHII, master_mass_frac,
                                     reion_plots["fixed_XHI_values"],
                                     reion_plots["small_scale_def"],
                                     reion_plots["large_scale_def"],
                                     model_tags, output_dir, "ps_scales",
                                     output_format)

    if reion_plots["slices_fixed_XHI"] and rank == 0:
        print("Plotting slices at fixed XHI fractions.")
        reionplot.plot_slices_XHI(reion_data["z_array_reion_allmodels"],
                                  reion_data["cosmology_allmodels"],
                                  master_mass_frac, 
                                  reion_data["XHII_fbase_allmodels"],
                                  reion_data["XHII_precision_allmodels"],
                                  reion_data["GridSize_allmodels"],
                                  reion_data["boxsize_allmodels"],
                                  reion_data["first_snap_allmodels"],
                                  reion_plots["fixed_XHI_values"],
                                  reion_plots["cut_slice"],
                                  reion_plots["cut_thickness"],
                                  model_tags, output_dir, "slices_XHI",
                                  output_format)


    if reion_plots["bubble_size"] and rank == 0:
        print("Determining bubble sizes at fixed XHI.")
        reionplot.determine_bubble_size(reion_data["z_array_reion_allmodels"],
                                        master_mass_frac,
                                        reion_data["first_snap_allmodels"],
                                        reion_data["GridSize_allmodels"],
                                        reion_data["boxsize_allmodels"],
                                        reion_plots["fixed_XHI_values"],
                                        model_tags, output_dir)


def generate_data(rank, size, comm, reion_ini_files, gal_ini_files,
                  reion_plots, output_dir, model_tags):
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

    reion_ini_files, gal_ini_files : List of strings 
        ``.ini`` file corresponding to each model that we're plotting.  We need
        both the galaxy (``SAGE``) and reionization (``cifog``) ``.ini`` files. 

    reion_plots : Dictionary
        Controls which of the plots we will make.  Keys are the name of each
        plot (e.g., ``reion``) and the value specifies if we are plotting it. If
        we're not plotting a property we don't need to calculate stuff for it! 

    output_dir : String
        Directory where the plots are saved. Used to save MC data.

    model_tags : List of strings
        String that will appear on the legend of the plot for each model. Used
        to save MC data with a unique name.

    Returns
    ---------

    reion_data : Dictionary
        All of the calculated properties required to create the reionization
        plots.
    """

    if rank == 0:
        print("Generating reionization data for a total of {0} "
              "models.".format(len(reion_ini_files)))

    # ======================================================================= #
    # We calculate values for all models and put them into lists that are     #
    # indexed by ``model_number``. So first we need to set up the outer-lists #
    # then we will append to these for each model.                            #
    # ======================================================================= #
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
    last_snap_allmodels = []

    GridSize_allmodels = []
    boxsize_allmodels = []
    helium_allmodels = []
    nion_factor_allmodels = []

    cosmology_allmodels = []
    t_bigbang_allmodels = []

    # Be careful, we use neutral fraction values here.
    volume_frac_allmodels = []
    mass_frac_allmodels = []

    # These are the nion grids used in cifog, so these are **actual** escaping
    # ionizing photons.
    nion_allmodels = []

    # Power Spectra.
    k_allmodels = []
    P21_allmodels = []
    PHII_allmodels = []

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
        last_snap_allmodels.append(last_snap)

        GridSize = int(SAGE_params["GridSize"])
        GridSize_allmodels.append(GridSize)

        # Careful, cifog uses Mpc/h.
        boxsize = float(SAGE_params["BoxSize"])
        boxsize_allmodels.append(boxsize)
        # However we use the volume as Mpc^3.
        model_volume = pow(float(SAGE_params["BoxSize"]) / \
                           float(SAGE_params["Hubble_h"]),3)

        helium = float(cifog_params["Y"])
        helium_allmodels.append(helium)

        nion_factor = float(cifog_params["nion_factor"])
        nion_factor_allmodels.append(nion_factor) 

        model_hubble_h = float(SAGE_params["Hubble_h"])
        model_halopartcut = int(SAGE_params["HaloPartCut"])

        # Load the redshift file and calculate the lookback times. 
        z_array_full, lookback_array_full = gd.load_redshifts(SAGE_params["FileWithSnapList"],
                                                              cosmology, t_bigbang)
        z_array_reion = np.array(z_array_full[first_snap:last_snap])
        lookback_array_reion = np.array(lookback_array_full[first_snap:last_snap])

        z_array_reion_allmodels.append(z_array_reion)
        lookback_array_reion_allmodels.append(lookback_array_reion)

        # Determine the base file names for the ionization, ionizing photons
        # and density fields. 
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

        # The ionization fields are assumed to have double precision.
        XHII_precision = 2
        XHII_precision_allmodels.append(XHII_precision)

        # Now it's time to set up all the arrays for this model number. 
        volume_frac_allmodels.append(np.zeros(last_snap - first_snap))
        mass_frac_allmodels.append(np.zeros(last_snap - first_snap))

        nion_allmodels.append(np.zeros(last_snap - first_snap))

        k_allmodels.append([])
        P21_allmodels.append([])
        PHII_allmodels.append([])

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
            volume_frac = 1.0 - np.mean(XHII)
            mass_frac = 1.0 - np.sum(XHII * density / np.sum(density))

            volume_frac_allmodels[model_number][snap_idx] = volume_frac
            mass_frac_allmodels[model_number][snap_idx] = mass_frac

            # Only need ionizing photons if we're plotting it.
            if reion_plots["nion"]:
                nion_path = "{0}_{1:03d}".format(nion_fbase, snapnum)
                nion = rs.read_binary_grid(nion_path, GridSize, nion_precision)

                nion_allmodels[model_number][snap_idx] = np.sum(nion)

            # If we're plotting the power spectra in scale space need to
            # calculate them at every single snapshot.
            if reion_plots["ps_scales"]:
                T0 = T_naught(z_array_reion[snap_idx], cosmology.H(0).value/100.0,
                              cosmology.Om0, cosmology.Ob0)

                # Be aware, using boxsize in Mpc/h.
                tmp_k, tmp_PowSpec, tmp_Error, \
                tmp_k_XHII, tmp_Pspec_HII, tmp_Error_XHII = calc_ps(XHII, density,
                                                                    boxsize)

                k_allmodels[model_number].append(tmp_k)
                P21_allmodels[model_number].append(tmp_PowSpec * T0*T0 * tmp_k**3 * 4.0*np.pi)
                PHII_allmodels[model_number].append(tmp_Pspec_HII)
           
            if reion_plots["bubble_size"] and \
               (mass_frac < 0.95 and mass_frac > 0.05):

                # Only calculate the MC if the file doesn't exist.
                MC_path = "{0}/MC/{1}_z_{2:.3f}.txt".format(output_dir,
                                                            model_tags[model_number],
                                                            z_array_reion[snap_idx]) 
                if (os.path.exists(infile) == False):
                    calculate_bubble_MC(XHII, GridSize, MC_path)
            
        # Snapshot Loop.

        # Ionizing emissitivty is scaled by the simulation volume (in Mpc^3).
        nion_allmodels[model_number] /= model_volume

    # Model Number Loop.

    # Everything has been calculated. Now construct a dictionary that contains
    # all the data (for easy passing) and return it. 
    reion_data = {"z_array_reion_allmodels" : z_array_reion_allmodels,
                  "lookback_array_reion_allmodels" : lookback_array_reion_allmodels,
                  "cosmology_allmodels" : cosmology_allmodels,
                  "t_bigbang_allmodels" : t_bigbang_allmodels,    
                  "volume_frac_allmodels" : volume_frac_allmodels,
                  "mass_frac_allmodels" : mass_frac_allmodels,
                  "nion_allmodels" : nion_allmodels,
                  "k_allmodels" : k_allmodels,
                  "P21_allmodels" : P21_allmodels,
                  "PHII_allmodels" : PHII_allmodels,
                  "XHII_fbase_allmodels" : XHII_fbase_allmodels,
                  "XHII_precision_allmodels" : XHII_precision_allmodels,
                  "density_fbase_allmodels" : density_fbase_allmodels,
                  "density_precision_allmodels" : density_precision_allmodels,
                  "GridSize_allmodels" : GridSize_allmodels,
                  "boxsize_allmodels" : boxsize_allmodels,
                  "helium_allmodels" : helium_allmodels,
                  "nion_factor_allmodels" : nion_factor_allmodels,
                  "first_snap_allmodels" : first_snap_allmodels,
                  "last_snap_allmodels" : last_snap_allmodels}

    return reion_data
