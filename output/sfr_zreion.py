
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import ReadScripts as rs
import PlotScripts as ps

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

output_format = "png"

def load_redshifts(fname):

    with open(fname, "r") as f:

        a = np.loadtxt(f)

    z = 1.0/a - 1.0

    return z 


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


def reduce_hist_across_tasks(hists):

    master_hists = []

    for model_number in range(len(hists)):

        if rank == 0:
            counts_total = np.zeros_like(hists[model_number])
        else:
            counts_total = None

        comm.Reduce([hists[model_number], MPI.FLOAT], [counts_total, MPI.FLOAT], op = MPI.SUM, root = 0) # Sum all the stellar mass and pass to Rank 0.

        master_hists.append(counts_total)

    return master_hists


def plot_zreion_hist(zreion_hist_allmodels, z_array_allmodels, model_tags,
                     output_tag="zreion_hist"):

    master_counts = reduce_hist_across_tasks(zreion_hist_allmodels)

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

    outputFile1 = "./{0}.{1}".format(output_tag, output_format) 
    fig1.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig1)


if __name__ == "__main__":

    ps.Set_Params_Plot()

    ini_file_model1="/fred/oz004/jseiler/kali/self_consistent_output/rsage_constant/ini_files/const_0.3_SAGE.ini"
    ini_file_model2="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.40_beta0.05_SAGE.ini"
    ini_file_model3="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHneg/ini_files/MHneg_1e8_1e12_0.99_0.05_SAGE.ini"
    ini_file_model4="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHpos/ini_files/MHpos_1e8_1e12_0.01_0.50_SAGE.ini"

    ini_files = [ini_file_model1,
                 ini_file_model2,
                 ini_file_model3,
                 ini_file_model4]

    model_tags = [r"$\mathbf{f_\mathrm{esc} \: Constant}$", 
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}^{-1}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}}$"]    

    z_array_allmodels = []    

    zreion_bins_allmodels = []
    zreion_hist_total_allmodels = []


    for model_number, ini_file in enumerate(ini_files):
        SAGE_params = rs.read_SAGE_ini(ini_file)

        z_array_full = load_redshifts(SAGE_params["FileWithSnapList"])
        z_array_full = np.array(z_array_full[0:int(SAGE_params["LastSnapShotNr"])+1])
        z_array_reion = np.array(z_array_full[int(SAGE_params["LowSnap"]):int(SAGE_params["LastSnapShotNr"])+1])
        z_array_allmodels.append(z_array_reion)

        # We want to bin the galaxies on zreion.  The bins need to be
        # monotonically increasing so flip the array.
        zreion_bins = np.flip(z_array_reion, axis=-1)
        # However we want the first bin to be [0.0, LastSimulationRedshift)
        # so add it on.
        zreion_bins = np.insert(zreion_bins, 0, 0.0) 
        # Add to the ``allmodels`` array.
        zreion_bins_allmodels.append(zreion_bins[:-1]) 

        GridSize = int(SAGE_params["GridSize"])
 
        galaxy_name = "{0}/{1}_z{2:.3f}".format(SAGE_params["OutputDir"],
                                                SAGE_params["FileNameGalaxies"],
                                                z_array_reion[-1]) 

        merged_name = "{0}/{1}_MergedGalaxies".format(SAGE_params["OutputDir"],
                                                      SAGE_params["FileNameGalaxies"])

        zreion_name = "{0}/{1}".format(SAGE_params["PhotoionDir"],
                                       SAGE_params["ReionRedshiftName"])

        zreion = rs.read_binary_grid(zreion_name, GridSize, 2, reshape=False) 

        for fnr in range(int(SAGE_params["FirstFile"]) + rank,
                         int(SAGE_params["LastFile"])+1, size):
                         #1, size):

            GG, Gal_Desc = rs.ReadGals_SAGE(galaxy_name, fnr, len(z_array_full))
            G_Merged, _ = rs.ReadGals_SAGE(merged_name, fnr, len(z_array_full))
            G = rs.Join_Arrays(GG, G_Merged, Gal_Desc)
            
            central_inds = get_central_inds(G.GridType) 
            GridPos = get_last_gridpos(G.GridHistory[central_inds]) 

            zreion_gals = zreion[GridPos]

            # Do the binning. 
            zreion_hist = np.histogram(zreion_gals, bins=zreion_bins)

            try:
                zreion_hist_total_allmodels[model_number] += zreion_hist[0]
            except IndexError:
                zreion_hist_total_allmodels.append(zreion_hist[0])

    plot_zreion_hist(zreion_hist_total_allmodels, z_array_allmodels, model_tags)
    
