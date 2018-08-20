
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from astropy import units as u
from astropy import cosmology


import ReadScripts as rs
import PlotScripts as ps
import CollectiveStats as collective

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


def get_last_gridpos(GridHistory):

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


def plot_zreion_hist(zreion_hist_allmodels, z_array_allmodels, model_tags,
                     output_dir, output_tag): 

    master_counts = collective.reduce_hist_across_tasks(rank, comm, 
                                                        zreion_hist_allmodels)

    if rank != 0:
        return 

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    for model_number in range(len(z_array_allmodels)):
        ax1.plot(z_array_allmodels[model_number], 
                 master_counts[model_number],
                 color = ps.colors[model_number],
                 ls = ps.linestyles[model_number],
                 label = model_tags[model_number])

    ax1.set_yscale("log")

    ax1.set_xlabel(r"$\mathbf{z_{reion}}$",
                   fontsize = ps.global_labelsize)
    ax1.set_ylabel(r"$\mathbf{Count}$",
                   fontsize = ps.global_labelsize)

    leg = ax1.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(ps.global_legendsize)

    outputFile1 = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    fig1.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig1)


def plot_zreion_SFR(mean_zreion_SFR_allmodels, std_zreion_SFR_allmodels,
                    N_zreion_SFR_allmodels, zreion_bins_allmodels,
                    zreion_array_allmodels,
                    model_tags, output_dir, output_tag):   

    master_mean, master_std, master_N, _ = \
        collective.collect_across_tasks(rank, comm, \
                                        mean_zreion_SFR_allmodels, \
                                        std_zreion_SFR_allmodels, \
                                        N_zreion_SFR_allmodels, \
                                        zreion_bins_allmodels, \
                                        zreion_bins_allmodels, \
                                        binned=True)

    if rank != 0:
        return 

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)


    for model_number in range(len(zreion_bins_allmodels)):
        print(master_mean[model_number][10]) 
        ax1.plot(zreion_array_allmodels[model_number], 
                 master_mean[model_number][10],
                 color = ps.colors[model_number],
                 ls = ps.linestyles[model_number],
                 label = model_tags[model_number])

        z_label = r"$\mathbf{z = " + \
                    str(round(zreion_array_allmodels[model_number][10],2)) + \
                   "}$"                
        ax1.text(12.5, 1e-2, z_label, fontsize=ps.global_labelsize-4) 

    ax1.set_yscale("log")

    ax1.set_xlabel(r"$\mathbf{z}$",
                   fontsize = ps.global_labelsize)
    ax1.set_ylabel(r"$\mathbf{Mean SFR}$",
                   fontsize = ps.global_labelsize)

    leg = ax1.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(ps.global_legendsize)

    outputFile1 = "{0}/{1}.{2}".format(output_dir, output_tag, output_format)
    fig1.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig1)


def plot_galaxy_properties(ini_files, model_tags, galaxy_plots, output_dir="."):


    galaxy_data = read_data(ini_files)

    if galaxy_plots["zreion_hist"]:
        plot_zreion_hist(galaxy_data["zreion_hist_total_allmodels"],
                         galaxy_data["zreion_bins_allmodels"], 
                         model_tags, output_dir, "zreion_hist")

    if galaxy_plots["zreion_sfr"]:
        plot_zreion_SFR(galaxy_data["mean_zreion_SFR_allmodels"],
                        galaxy_data["std_zreion_SFR_allmodels"],
                        galaxy_data["N_zreion_SFR_allmodels"],
                        galaxy_data["zreion_bins_allmodels"],
                        galaxy_data["z_array_reion_allmodels"],
                        model_tags, output_dir, "zreion_SFR")


def read_data(ini_files):

    z_array_reion_allmodels = []        
    lookback_array_reion_allmodels = []

    zreion_bins_allmodels = []
    zreion_hist_total_allmodels = []

    # These are the arrays for the SFR at each snapshot within each zreion bin.
    mean_zreion_SFR_allmodels = []
    std_zreion_SFR_allmodels = []
    N_zreion_SFR_allmodels = []

    for model_number, ini_file in enumerate(ini_files):
        SAGE_params = rs.read_SAGE_ini(ini_file)

        cosmology, t_BigBang = set_cosmology(float(SAGE_params["Hubble_h"]),
                                             float(SAGE_params["Omega"]))

        first_snap = int(SAGE_params["LowSnap"])
        last_snap = int(SAGE_params["LastSnapShotNr"])

        z_array_full, lookback_array_full = load_redshifts(SAGE_params["FileWithSnapList"],
                                                           cosmology, t_BigBang)
        z_array_full = np.array(z_array_full[0:last_snap+1])
        lookback_array_full = np.array(lookback_array_full[0:last_snap+1])


        z_array_reion = np.array(z_array_full[first_snap:last_snap+1])
        z_array_reion_allmodels.append(z_array_reion)

        lookback_array_reion = np.array(lookback_array_full[first_snap:last_snap+1])
        lookback_array_reion_allmodels.append(lookback_array_reion)

        GridSize = int(SAGE_params["GridSize"])
 
        galaxy_name = "{0}/{1}_z{2:.3f}".format(SAGE_params["OutputDir"],
                                                SAGE_params["FileNameGalaxies"],
                                                z_array_reion[-1]) 

        merged_name = "{0}/{1}_MergedGalaxies".format(SAGE_params["OutputDir"],
                                                      SAGE_params["FileNameGalaxies"])

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



        for fnr in range(int(SAGE_params["FirstFile"]) + rank,
                         int(SAGE_params["LastFile"])+1, size):
#                         1, size):

                


            GG, Gal_Desc = rs.ReadGals_SAGE(galaxy_name, fnr, len(z_array_full))
            G_Merged, _ = rs.ReadGals_SAGE(merged_name, fnr, len(z_array_full))
            G = rs.Join_Arrays(GG, G_Merged, Gal_Desc)
            
            central_inds = get_central_inds(G.GridType)
            G_centrals = G[central_inds]
 
            GridPos = get_last_gridpos(G_centrals.GridHistory) 
            zreion_gals = zreion[GridPos]

            # Do the binning. 
            zreion_hist = np.histogram(zreion_gals, bins=zreion_bins)

            # Now when we bin with SFR we want to bunch it together...
             
            zreion_inds = np.digitize(zreion_gals, bins=zreion_bins)

            # Now within each zreion bin we want to calculate the mean/std
            # galaxy SFR across ALL snapshots.  So let's go through each bin
            # and first find all the galaxies within that bin. 
            #for zreion_idx in [63]:
            for zreion_idx in range(max(zreion_inds)-1): 
                Gals_in_bin_inds = np.where(zreion_inds == zreion_idx)[0]
                Gals_in_bin = G_centrals[Gals_in_bin_inds]

                mean_SFR_this_bin = np.zeros(len(z_array_reion),dtype=np.float32)
                std_SFR_this_bin = np.zeros(len(z_array_reion),dtype=np.float32)
                N_SFR_this_bin = np.zeros(len(z_array_reion),dtype=np.float32)

                # Now that we have each galaxy in this bin, calculate the mean
                # SFR across all snapshots.
                for snap_count, snapnum in enumerate(np.arange(first_snap, 
                                                               len(z_array_reion))):
                    Gals_exist = np.where(Gals_in_bin.GridHistory[:, snapnum] != -1)[0]

                    mean_SFR_this_bin[snap_count] = np.mean(Gals_in_bin.GridSFR[Gals_exist, snapnum])
                    std_SFR_this_bin[snap_count] = np.std(Gals_in_bin.GridSFR[Gals_exist, snapnum])
                    N_SFR_this_bin[snap_count] = len(Gals_exist)

                # Update the running totals for this zreion bin.
                mean_zreion_SFR_allmodels[model_number][zreion_idx], \
                std_zreion_SFR_allmodels[model_number][zreion_idx] = \
                collective.update_cumulative_stats( 
                    mean_zreion_SFR_allmodels[model_number][zreion_idx], \
                    std_zreion_SFR_allmodels[model_number][zreion_idx], \
                    N_zreion_SFR_allmodels[model_number][zreion_idx], \
                    mean_SFR_this_bin, std_SFR_this_bin, N_SFR_this_bin) 

                N_zreion_SFR_allmodels[model_number][zreion_idx] += N_SFR_this_bin
            # Proceed to next zreion bin.

            zreion_hist_total_allmodels[model_number] += zreion_hist[0]

            print(len(mean_zreion_SFR_allmodels[model_number]))
            print(len(zreion_bins_allmodels[model_number]))

    galaxy_data = {"z_array_reion_allmodels" : z_array_reion_allmodels,
                   "lookback_array_reion_allmodels" : lookback_array_reion_allmodels,
                   "zreion_bins_allmodels" : zreion_bins_allmodels,
                   "zreion_hist_total_allmodels" : zreion_hist_total_allmodels,
                   "mean_zreion_SFR_allmodels" : mean_zreion_SFR_allmodels,
                   "std_zreion_SFR_allmodels" : std_zreion_SFR_allmodels,
                   "N_zreion_SFR_allmodels" : N_zreion_SFR_allmodels} 


    return galaxy_data
 
if __name__ == "__main__":

    ps.Set_Params_Plot()

    ini_file_model1="/fred/oz004/jseiler/kali/self_consistent_output/rsage_constant/ini_files/const_0.3_SAGE.ini"
    ini_file_model2="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.40_beta0.05_SAGE.ini"
    ini_file_model3="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHneg/ini_files/MHneg_1e8_1e12_0.99_0.05_SAGE.ini"
    ini_file_model4="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHpos/ini_files/MHpos_1e8_1e12_0.01_0.50_SAGE.ini"

    ini_files = [ini_file_model1]
                 #ini_file_model2,
                 #ini_file_model3,
                 #ini_file_model4]

    model_tags = [r"$\mathbf{f_\mathrm{esc} \: Constant}$", 
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}^{-1}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}}$"]    

    zreion_hist = 0
    zreion_sfr = 1

    galaxy_plots = {"zreion_hist" : zreion_hist,
                    "zreion_sfr" : zreion_sfr}

    plot_galaxy_properties(ini_files, model_tags, galaxy_plots)     
    
