"""
Author: Jacob Seiler
"""
import numpy as np

from astropy import units as u
from astropy import cosmology

from scipy import stats

import ReadScripts as rs
import PlotScripts as ps
import CollectiveStats as collective
import GalaxyPlots as galplot

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

output_format = "png"


def set_cosmology(Hubble_h, Omega_m):
    
    cosmo = cosmology.FlatLambdaCDM(H0 = Hubble_h*100, Om0 = Omega_m) 
    t_BigBang = (cosmo.lookback_time(100000).value)*1.0e3 # Lookback time to the Big Bang in Myr.

    return cosmo, t_BigBang


def load_redshifts(fname, cosmology, t_BigBang):

    with open(fname, "r") as f:

        a = np.loadtxt(f)

    z = 1.0/a - 1.0

    lookback = (t_BigBang - cosmology.lookback_time(z).value*1.0e3)  # In Myr.

    return z, lookback 


def get_first_gridpos(GridHistory):

    GridPos = np.full(len(GridHistory), -1)

    for snapnum in range(len(GridHistory[0,:])):

        w_exist = np.where(GridHistory[:, snapnum] != -1)[0]
        w_not_updated = np.where(GridPos == -1)[0]

        w_to_update = w_exist[np.isin(w_exist, w_not_updated)]

        GridPos[w_to_update] = GridHistory[w_to_update, snapnum]

    return GridPos 


def get_central_inds(GridType):

    central_inds = []

    for gal_idx in range(len(GridType)):

        w_sat = np.where(GridType[gal_idx, :] == 1)[0]
        if len(w_sat) == 0: 
            central_inds.append(gal_idx)

    return central_inds 


def calculate_zreion_SFR(Galaxies, mean_zreion_SFR_thismodel, 
                         std_zreion_SFR_thismodel, N_zreion_SFR_thismodel,
                         zreion_inds, z_array_reion, first_snap, last_snap,
                         mass_cut):

    # Now within each zreion bin we want to calculate the mean/std
    # galaxy SFR across ALL snapshots.  So let's go through each bin
    # and first find all the galaxies within that bin. 
    for zreion_idx in range(max(zreion_inds)-1): 
        Gals_in_bin_inds = np.where(zreion_inds == zreion_idx)[0]
        Gals_in_bin = Galaxies[Gals_in_bin_inds]

        mean_SFR_this_bin = np.zeros(len(z_array_reion),dtype=np.float32)
        std_SFR_this_bin = np.zeros(len(z_array_reion),dtype=np.float32)
        N_SFR_this_bin = np.zeros(len(z_array_reion),dtype=np.float32)

        # Now that we have each galaxy in this bin, calculate the mean
        # SFR across all snapshots.
        for snap_count, snapnum in enumerate(np.arange(first_snap, 
                                                       last_snap + 1)): 
            Gals_exist = np.where((Gals_in_bin.GridHistory[:, snapnum] != -1) &
                                  (Gals_in_bin.GridFoFMass[:, snapnum] < mass_cut))[0]

            mean_SFR_this_bin[snap_count] = np.mean(Gals_in_bin.GridSFR[Gals_exist, snapnum])
            std_SFR_this_bin[snap_count] = np.std(Gals_in_bin.GridSFR[Gals_exist, snapnum])
            N_SFR_this_bin[snap_count] = len(Gals_exist)

        # Update the running totals for this zreion bin.
        mean_zreion_SFR_thismodel[zreion_idx], \
        std_zreion_SFR_thismodel[zreion_idx] = \
        collective.update_cum_stats( 
            mean_zreion_SFR_thismodel[zreion_idx], \
            std_zreion_SFR_thismodel[zreion_idx], \
            N_zreion_SFR_thismodel[zreion_idx], \
            mean_SFR_this_bin, std_SFR_this_bin, N_SFR_this_bin) 

        N_zreion_SFR_thismodel[zreion_idx] += N_SFR_this_bin
    # Proceed to next zreion bin.

    return mean_zreion_SFR_thismodel, std_zreion_SFR_thismodel, \
           N_zreion_SFR_thismodel


def do_2D_binning(data_x, data_y, mean_curr, std_curr, N_curr, bins):

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


def plot_galaxy_properties(ini_files, model_tags, galaxy_plots, output_dir="."):


    # First calculate all the properties and statistics we need.
    galaxy_data = read_data(ini_files, galaxy_plots)

    # Then find what plots we need and plot em!
    if galaxy_plots["zreion_hist"]:
        galplot.plot_zreion_hist(rank, comm,
                                 galaxy_data["zreion_hist_total_allmodels"],
                                 galaxy_data["zreion_bins_allmodels"], 
                                 model_tags, output_dir, "zreion_hist")

    if galaxy_plots["zreion_sfr"]:
        galplot.plot_zreion_SFR(rank, comm,
                                galaxy_data["mean_zreion_SFR_allmodels"],
                                galaxy_data["std_zreion_SFR_allmodels"],
                                galaxy_data["N_zreion_SFR_allmodels"],
                                galaxy_data["zreion_bins_allmodels"],
                        galaxy_data["z_array_reion_allmodels"],
                        model_tags, output_dir, "zreion_SFR_1e9")

    if galaxy_plots["nion"]:
        galplot.plot_nion(rank, comm,
                          galaxy_data["z_array_full_allmodels"],
                          galaxy_data["lookback_array_full_allmodels"],
                          galaxy_data["sum_nion_allmodels"],
                          galaxy_data["cosmology_allmodels"],
                          galaxy_data["t_bigbang_allmodels"],
                          model_tags, output_dir, "nion")

    if galaxy_plots["mstar_fesc"]:
        galplot.plot_mstar_fesc(rank, comm,
                                galaxy_data["mstar_bins"],
                                galaxy_data["mstar_bin_width"],
                                galaxy_data["z_array_full_allmodels"],
                                galaxy_data["mean_mstar_fesc_allmodels"],
                                galaxy_data["std_mstar_fesc_allmodels"],
                                galaxy_data["N_mstar_fesc_allmodels"],
                                model_tags, output_dir, "mstar_fesc",
                                plot_snaps_for_models=galaxy_plots["plot_snaps_for_models"],
                                plot_models_at_snaps=galaxy_plots["plot_models_at_snaps"])

    if galaxy_plots["SMF"]:
        galplot.plot_SMF(rank, comm,
                         galaxy_data["mstar_bins"],
                         galaxy_data["mstar_bin_width"],
                         galaxy_data["z_array_full_allmodels"],
                         galaxy_data["SMF_allmodels"],
                         galaxy_data["cosmology_allmodels"],
                         model_tags, output_dir, "SMF_gridinfall")

    if galaxy_plots["mstar_fej"]:
        galplot.plot_mstar_fej(rank, comm,
                               galaxy_data["mstar_bins"],
                               galaxy_data["mstar_bin_width"],
                               galaxy_data["z_array_full_allmodels"],
                               galaxy_data["mean_mstar_fej_allmodels"],
                               galaxy_data["std_mstar_fej_allmodels"],
                               galaxy_data["N_mstar_fej_allmodels"],
                               model_tags, output_dir, "mstar",
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
                               plot_snaps_for_models=galaxy_plots["plot_snaps_for_models"],
                               plot_models_at_snaps=galaxy_plots["plot_models_at_snaps"])

    if galaxy_plots["mstar_infall"]:
        galplot.plot_mstar_infall(rank, comm,
                               galaxy_data["mstar_bins"],
                               galaxy_data["mstar_bin_width"],
                               galaxy_data["z_array_full_allmodels"],
                               galaxy_data["mean_mstar_infall_allmodels"],
                               galaxy_data["std_mstar_infall_allmodels"],
                               galaxy_data["N_mstar_infall_allmodels"],
                               model_tags, output_dir, "mstar_infall",
                               plot_snaps_for_models=galaxy_plots["plot_snaps_for_models"],
                               plot_models_at_snaps=galaxy_plots["plot_models_at_snaps"])


def read_data(ini_files, galaxy_plots):

    # General stuff for each model.
    z_array_full_allmodels = []
    lookback_array_full_allmodels = []

    z_array_reion_allmodels = []        
    lookback_array_reion_allmodels = []

    zreion_bins_allmodels = []
    zreion_hist_total_allmodels = []

    cosmology_allmodels = []
    t_bigbang_allmodels = []

    # These are the arrays for the SFR at each snapshot within each zreion bin.
    mean_zreion_SFR_allmodels = []
    std_zreion_SFR_allmodels = []
    N_zreion_SFR_allmodels = []

    # These are the arrays for the number of ionizing photons at each snapshot.
    # Note: This is the ESCAPING ionizing photons. 
    sum_nion_allmodels = []

    # Binning parameters for stellar mass. 
    mstar_bin_low = 5.0
    mstar_bin_high = 12.0
    mstar_bin_width = 0.2
    mstar_Nbins = int((mstar_bin_high - mstar_bin_low) / mstar_bin_width)
    mstar_bins = np.arange(mstar_bin_low, 
                           mstar_bin_high + mstar_bin_width,
                           mstar_bin_width)

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

    # Infall rate as a function of stellar mass (Mstar). 
    mean_mstar_infall_allmodels = []
    std_mstar_infall_allmodels = []
    N_mstar_infall_allmodels = []

    # All outer arrays set up, time to read in the data!
    for model_number, ini_file in enumerate(ini_files):

        # Read in the parameters and set some initial variables.
        SAGE_params = rs.read_SAGE_ini(ini_file)

        cosmology, t_bigbang = set_cosmology(float(SAGE_params["Hubble_h"]),
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

        mass_cut = 1.0e9/1.0e10 * float(SAGE_params["Hubble_h"])

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

        # Read the reionization redshift grid.
        zreion_name = "{0}/{1}".format(SAGE_params["PhotoionDir"],
                                       SAGE_params["ReionRedshiftName"])

        zreion = rs.read_binary_grid(zreion_name, GridSize, 2, reshape=False)

        # We want to bin the galaxies on zreion.  Do this by splitting them
        # (roughly) up into 50Myr intervals.
        zreion_bins = [z_array_reion[0]]
        time = 0.0
        for count in np.arange(1, len(lookback_array_reion)-1):
           time += lookback_array_reion[count] - lookback_array_reion[count-1] 
           if time > 50.0:
                zreion_bins.append(z_array_reion[count])
                time = 0.0

        # Want the final bin to be the final redshift. 
        zreion_bins.append(z_array_reion[-1])
        # The bins need to be monotonically increasing, so flip them.
        zreion_bins = np.flip(zreion_bins, axis=-1)
        # Add to the ``allmodels`` array.
        zreion_bins_allmodels.append(zreion_bins[:-1]) 
        # Initialize the array that will hold all the histograms.
        zreion_hist_total_allmodels.append(np.zeros(len(zreion_bins[:-1]),
                                                    dtype=np.float32))

        # For the SFR at each snapshot within each zreion bin, each zreion bin
        # will contain an array with len(NumberSnapshots) that has the
        # mean/std/N data points.
        mean_zreion_SFR_allmodels.append([])
        std_zreion_SFR_allmodels.append([])
        N_zreion_SFR_allmodels.append([])
        for bin_count in range(len(zreion_bins)-1):
            mean_zreion_SFR_allmodels[model_number].append(np.zeros(len(z_array_reion),
                                                                    dtype=np.float32))
            std_zreion_SFR_allmodels[model_number].append(np.zeros(len(z_array_reion),
                                                                   dtype=np.float32))
            N_zreion_SFR_allmodels[model_number].append(np.zeros(len(z_array_reion),
                                                                 dtype=np.float32))

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

        # Star formation rate as a function of stellar mass.
        mean_mstar_infall_allmodels.append([]) 
        std_mstar_infall_allmodels.append([])
        N_mstar_infall_allmodels.append([]) 
        for snap_count in range(len(z_array_full)):
            mean_mstar_infall_allmodels[model_number].append(np.zeros(mstar_Nbins,
                                                                   dtype=np.float32))
            std_mstar_infall_allmodels[model_number].append(np.zeros(mstar_Nbins,
                                                                  dtype=np.float32))
            N_mstar_infall_allmodels[model_number].append(np.zeros(mstar_Nbins,
                                                                dtype=np.float32))


        # ========================================================= #
        # Now go through each file and calculate the stuff we need. #
        # ========================================================= #
        for fnr in range(int(SAGE_params["FirstFile"]) + rank,
                         int(SAGE_params["LastFile"])+1, size):
#                         1, size):


            print("Rank {0}: Model {1} File {2}".format(rank, model_number,
                                                        fnr))

            GG, Gal_Desc = rs.ReadGals_SAGE(galaxy_name, fnr, len(z_array_full))
            G_Merged, _ = rs.ReadGals_SAGE(merged_name, fnr, len(z_array_full))
            G = rs.Join_Arrays(GG, G_Merged, Gal_Desc)

            # Let's only deal with central galaxies.                
            central_inds = get_central_inds(G.GridType)
            G_centrals = G[central_inds]

            # Find the first position of each galaxy (assume it doesn't move
            # over its history). 
            GridPos = get_first_gridpos(G_centrals.GridHistory) 
            zreion_gals = zreion[GridPos]

            # Do the binning. 
            zreion_hist = np.histogram(zreion_gals, bins=zreion_bins)

            # Now when we bin with SFR we want to bunch it together...             
            zreion_inds = np.digitize(zreion_gals, bins=zreion_bins)

            # Calculate the SFR within each zreion bin.
            if galaxy_plots["zreion_sfr"]: 
                mean_zreion_SFR_allmodels[model_number], \
                std_zreion_SFR_allmodels[model_number], \
                N_zreion_SFR_allmodels[model_number] = \
                    calculate_zreion_SFR(G_centrals, 
                                         mean_zreion_SFR_allmodels[model_number],
                                         std_zreion_SFR_allmodels[model_number],
                                         N_zreion_SFR_allmodels[model_number],
                                         zreion_inds, z_array_reion,
                                         first_snap, last_snap, mass_cut)

            # Histogram count of the zreion values. 
            zreion_hist_total_allmodels[model_number] += zreion_hist[0]

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
                infall = G.GridInfallRate[Gals_exist, snapnum]

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

                if galaxy_plots["mstar_infall"]:   
                    mean_mstar_infall_allmodels[model_number][snap_count], \
                    std_mstar_infall_allmodels[model_number][snap_count], \
                    N_mstar_infall_allmodels[model_number][snap_count] = \
                        do_2D_binning(log_mass, infall,
                                      mean_mstar_infall_allmodels[model_number][snap_count],
                                      std_mstar_infall_allmodels[model_number][snap_count],
                                      N_mstar_infall_allmodels[model_number][snap_count],
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
 
    galaxy_data = {"z_array_full_allmodels" : z_array_full_allmodels,
                   "lookback_array_full_allmodels" : lookback_array_full_allmodels,
                   "z_array_reion_allmodels" : z_array_reion_allmodels,
                   "lookback_array_reion_allmodels" : lookback_array_reion_allmodels,
                   "cosmology_allmodels" : cosmology_allmodels,
                   "t_bigbang_allmodels" : t_bigbang_allmodels,
                   "zreion_bins_allmodels" : zreion_bins_allmodels,
                   "zreion_hist_total_allmodels" : zreion_hist_total_allmodels,
                   "mean_zreion_SFR_allmodels" : mean_zreion_SFR_allmodels,
                   "std_zreion_SFR_allmodels" : std_zreion_SFR_allmodels,
                   "N_zreion_SFR_allmodels" : N_zreion_SFR_allmodels, 
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
                   "N_mstar_SFR_allmodels" : N_mstar_SFR_allmodels,
                   "mean_mstar_infall_allmodels" : mean_mstar_infall_allmodels,
                   "std_mstar_infall_allmodels" : std_mstar_infall_allmodels,
                   "N_mstar_infall_allmodels" : N_mstar_infall_allmodels}

    return galaxy_data

 
if __name__ == "__main__":

    ps.Set_Params_Plot()

    ini_file_model1="/fred/oz004/jseiler/kali/self_consistent_output/rsage_constant/ini_files/const_0.3_SAGE.ini"
    ini_file_model2="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHneg/ini_files/MHneg_1e8_1e12_0.90_0.05_SAGE.ini"   
    ini_file_model3="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHpos/ini_files/MHpos_1e8_1e12_0.01_0.50_SAGE.ini"
    ini_file_model4="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.40_beta0.05_SAGE.ini"
    ini_file_model5="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha1.00_beta1.00_delta1.00_SAGE.ini"
    ini_file_model6="/home/jseiler/rsage/ini_files/kali_noreion_SAGE.ini"

    #ini_file_infall="/home/jseiler/rsage/ini_files/kali_SAGE_gridinfall.ini"

    ini_files = [ini_file_model2,
                 ini_file_model3,
                 ini_file_model4,
                 ini_file_model5,
                 ini_file_model6]

    #ini_files = [ini_file_infall]

    model_tags = [r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}^{-1}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: SFR}$",
                  r"$\mathbf{No \: Reion}$"]

    zreion_hist = 0
    zreion_sfr = 0
    nion = 0
    mstar_fesc = 0 
    SMF = 0
    mstar_fej = 0
    mstar_SFR = 0
    mstar_infall = 0
    
    plot_snaps_for_models = [[33, 50, 76, 93]]
    plot_models_at_snaps = None


    galaxy_plots = {"zreion_hist" : zreion_hist,
                    "zreion_sfr" : zreion_sfr,
                    "nion" : nion,
                    "mstar_fesc" : mstar_fesc,
                    "SMF" : SMF,
                    "mstar_fej" : mstar_fej,
                    "mstar_SFR" : mstar_SFR,
                    "mstar_infall" : mstar_infall,
                    "plot_snaps_for_models" : plot_snaps_for_models,
                    "plot_models_at_snaps" : plot_models_at_snaps}

    plot_galaxy_properties(ini_files, model_tags, galaxy_plots)
