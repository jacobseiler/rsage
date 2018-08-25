"""
This file contains the plotting scripts for the Seiler et. al (2018) paper. We
also include a number of extra plots that can be turned on/off as desired.
Please refer to the ``README`` in this directory (``output``) for more
information on how to use this for your own data.
Author: Jacob Seiler
Version: 0.1
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


def plot_galaxy_properties(ini_files, model_tags, galaxy_plots, output_dir="."):
    """    
    Wrapper function to handle reading in of data + calculating galaxy
    properties, then calling the specified plotting routines.

    Parameters
    ----------

    ini_files : List of strings 
        ``.ini`` file corresponding to each model that we're plotting.

    model_tags : List of strings
        String that will appear on the legend of the plot for each model.

    galaxy_plots : Dictionary
        Controls which of the plots we will make.  Keys are the name of each
        plot (e.g., ``SMF``) and the value specifies if we are plotting it.

    output_dir : String, optional
        Directory where the plots are saved. If this directory does not exist,
        it is created beforehand. 

    Returns
    ---------

    None.
    """

    # First calculate all the properties and statistics we need.
    galaxy_data = generate_data(ini_files, galaxy_plots)

    # Check to see if the output directory exists.
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print("Made output directory {0}".format(output_dir))

    # Then find what plots we need and plot em!
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


def generate_data(ini_files, galaxy_plots):
    """    
    Reads in the galaxy data for calculate all the require properties for each
    models.

    Parameters
    ----------

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
    # We calculate values for all models and put them into lists that are
    # indexed by ``model_number``. So first we need to set up the outerl-lists
    # then we will append to these for each model. 
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
        # Parallelize over number of files.
        for fnr in range(int(SAGE_params["FirstFile"]) + rank,
                         int(SAGE_params["LastFile"])+1, size):
#                         1, size):

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
                   "N_mstar_SFR_allmodels" : N_mstar_SFR_allmodels,
                   "mean_mstar_infall_allmodels" : mean_mstar_infall_allmodels,
                   "std_mstar_infall_allmodels" : std_mstar_infall_allmodels,
                   "N_mstar_infall_allmodels" : N_mstar_infall_allmodels}

    return galaxy_data

 
if __name__ == "__main__":

    # Plotting is driven entirely through specifying the ``.ini`` files. 
    # For this reason, the directories in the ``.ini`` files **MUST** be
    # absolute paths, **NOT** relative. 
    ini_file_model1="/fred/oz004/jseiler/kali/self_consistent_output/rsage_constant/ini_files/const_0.3_SAGE.ini"
    ini_file_model2="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHneg/ini_files/MHneg_1e8_1e12_0.90_0.05_SAGE.ini"   
    ini_file_model3="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHpos/ini_files/MHpos_1e8_1e12_0.01_0.50_SAGE.ini"
    ini_file_model4="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.40_beta0.05_SAGE.ini"
    ini_file_model5="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha1.00_beta1.00_delta1.00_SAGE.ini"
    ini_file_model6="/home/jseiler/rsage/ini_files/kali_noreion_SAGE.ini"

    #ini_file_infall="/home/jseiler/rsage/ini_files/kali_SAGE_gridinfall.ini"

    # All ``.ini`` files included in this array will be plotted.
    ini_files = [ini_file_model2,
                 ini_file_model3,
                 ini_file_model4,
                 ini_file_model5,
                 ini_file_model6]

    #ini_files = [ini_file_infall]

    # These are the labels that will appear on the axis legends for each model.
    model_tags = [r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}^{-1}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: SFR}$",
                  r"$\mathbf{No \: Reion}$"]


    # ============================================================= #
    # Switches to control what plots to make. 0 to skip, 1 to plot. #
    # ============================================================= #
    nion = 0
    mstar_fesc = 0 
    SMF = 0
    mstar_fej = 0
    mstar_SFR = 0
    mstar_infall = 0

    # For some plots, there is the option of plotting one model at different
    # snapshots (within one panel) or plotting all models at one snapshot
    # (within one panel).  
    plot_snaps_for_models = [[33, 50, 76, 93]  # For each panel, plot one model 
                             [33, 50, 76, 93]  # at specified snapshots.
                             [33, 50, 76, 93]   
                             [33, 50, 76, 93]]  

    plot_models_at_snaps = None  # For each panel, plot one specified snapshot 
                                 # for all models.

    # ====================== #
    # Don't touch below here # 
    # ====================== #
    galaxy_plots = {"nion" : nion,
                    "mstar_fesc" : mstar_fesc,
                    "SMF" : SMF,
                    "mstar_fej" : mstar_fej,
                    "mstar_SFR" : mstar_SFR,
                    "mstar_infall" : mstar_infall,
                    "plot_snaps_for_models" : plot_snaps_for_models,
                    "plot_models_at_snaps" : plot_models_at_snaps}

    plot_galaxy_properties(ini_files, model_tags, galaxy_plots)
