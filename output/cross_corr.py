import os
import matplotlib
matplotlib.use('Agg')

import pylab as plt

import h5py
from matplotlib.colors import LogNorm
import numpy as np
from numpy.fft import fftn, ifftn

import ReadScripts
import AllVars
import PlotScripts

output_format = ".png"


def plot_crosscoeff(kbins, cross_coeff, output_tag): 

    conversion = (2*np.pi)**3/(2*np.pi**2)

    fig1 = plt.figure(figsize=(8,8))
    ax1 = fig1.add_subplot(111)

    ax1.axhline(y=0, xmin=0, xmax=10, color='k', ls = '--')

    #pspec = cross_pspec*kmid_bins**3*conversion 

    #ax1.set_ylim([min(pspec[1:-1]), max(pspec[1:-1])])
    #ax1.set_ylim([min(pspec[1:-1]), max(pspec[1:-1])])
    ax1.set_xscale("log")
    #ax1.set_yscale("symlog")

    #ax1.set_ylim([-0.1, 0.1])
    #ax1.plot(kmid_bins[1:-1], pspec[1:-1]) 
    ax1.plot(kmid_bins[1:-1], crosscorreff[1:-1]) 

    ax1.set_xlabel(r"$h [Mpc^{-1}h]$",
                   size = PlotScripts.global_labelsize)
    ax1.set_ylabel(r"$r_{\delta, N_{halo}}$",
                   size = PlotScripts.global_labelsize)
    #ax1.set_ylabel(r"$Cross Power$")

    plt.tight_layout()
    outputFile = "./{0}{1}".format(output_tag,
                                    output_format)                                    

    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))

    plt.close()


def generate_galaxy_density(galaxy_path, GridSize, snap, grid_name, num_files=64, num_snaps=99): 

    fname_gal="{0}_z5.782".format(galaxy_path)
    fname_merged="{0}_MergedGalaxies".format(galaxy_path)

    grid_path="./{0}_galaxygrid".format(grid_name)

    if os.path.isfile(grid_path):
        grid = np.load(grid_path)
        return grid        

    grid = np.zeros(pow(GridSize,3)) 

    for fnr in range(num_files): 
        GG, Gal_Desc = ReadScripts.ReadGals_SAGE(fname_gal, fnr, num_snaps) # Read galaxies
        G_Merged, _ = ReadScripts.ReadGals_SAGE(fname_merged, fnr, num_snaps) # Also need the merged galaxies.
        G = ReadScripts.Join_Arrays(GG, G_Merged, Gal_Desc) # Then join them together for all galaxies.

        # Find those galaxies that existed at the specified snapshot and their
        # locations.
        w_gal = np.where(G.GridHistory[:, snap] != -1)
        grid_inds = G.GridHistory[w_gal, snap]

        # Find the common grid locations and then add them their density to the
        # grid. 
        unique_inds, unique_count = np.unique(grid_inds, return_counts=True)
        grid[unique_inds] += unique_count

    print("Total of {0} galaxies".format(sum(grid)))

    grid.tofile(grid_path)
    print("Saved to {0}".format(grid_path))

    return grid


def generate_1D_grids(all_pos, GridSize, halos=0):

    def generate_grid(pos, GridSize):
        grid =  pos * GridSize / AllVars.BoxSize
        grid = np.around(grid)
        grid[grid >= GridSize] = GridSize - 1  
        grid = grid.astype(int)

        return grid

    x_grid = generate_grid(all_pos[:,0], GridSize)
    y_grid = generate_grid(all_pos[:,1], GridSize)
    z_grid = generate_grid(all_pos[:,2], GridSize)

    return x_grid, y_grid, z_grid


def generate_halo_density(halo_path, grid_name, GridSize=256, snap=97,
                          num_files=64, halo_mode=2, grid_mode="xyz",
                          check_file=0):
  

    if check_file:
        grid_path="./{0}_halogrid".format(grid_name)

        if os.path.isfile(grid_path):
            grid = np.load(grid_path)
            return grid        

 
    grid = np.zeros(pow(GridSize,3)) 

    for fnr in range(num_files): 

        if halo_mode == 2:
            fname = "{0}_{1:03d}.dat".format(halo_path, fnr) 
            Halos = ReadScripts.read_trees(fname)

            w = np.where(Halos["SnapNum"] == snap)[0]
            x_grid, y_grid, z_grid = generate_1D_grids(Halos["Pos"][w],
                                                       GridSize=GridSize)

            print("There are a total of {0} Halos in file and {1} Halos at " 
                  "Snapshot {2}".format(len(Halos), len(w), snap))

        if halo_mode == 3:
            fname = "{0}_{1:03d}.catalog_subgroups_properties/subfind_{1:03d}.catalog_subgroups_properties.{2}" \
                    .format(halo_path, snap, fnr)
            Halos = ReadScripts.read_subfind_halos(fname)

            x_grid, y_grid, z_grid = generate_1D_grids(Halos["position_COM"],
                                                       GridSize=GridSize)                                                       

        if grid_mode == "xyz":
            grid_inds = (x_grid*GridSize+y_grid)*GridSize+z_grid # xyz
        elif grid_mode == "xzy":
            grid_inds = (x_grid*GridSize+z_grid)*GridSize+y_grid # xzy
        elif grid_mode == "yxz":
            grid_inds = (y_grid*GridSize+x_grid)*GridSize+z_grid # yxz
        elif grid_mode == "yzx":
            grid_inds = (y_grid*GridSize+z_grid)*GridSize+x_grid # yzx
        elif grid_mode == "zxy":
            grid_inds = (z_grid*GridSize+x_grid)*GridSize+y_grid # zxy
        elif grid_mode == "zyx":
            grid_inds = (z_grid*GridSize+y_grid)*GridSize+x_grid # zyx
       
        unique_inds, unique_count = np.unique(grid_inds, return_counts=True)
        unique_inds = unique_inds.astype(int)
        grid[unique_inds] += unique_count

    if check_file:
        grid.tofile(grid_path) 
        print("Saved to {0}".format(grid_path))

    return grid


def generate_pseudo_density(path, GridSize, snap, num_files=256,
                            grid_mode="xyz"):

    grid = np.zeros(pow(GridSize,3))
 
    for fnr in range(num_files):
        print("Gridding file {0}".format(fnr))
        fname = "{0}_{1:03d}/kali_linker.{2}.hdf5".format(path, snap, fnr)

        with h5py.File(fname, "r") as f_in:
            pos = f_in["PartType1"]["Coordinates"][:]

        x_grid, y_grid, z_grid = generate_1D_grids(pos,
                                                   GridSize=GridSize)

        if grid_mode == "xyz":
            grid_inds = (x_grid*GridSize+y_grid)*GridSize+z_grid # xyz
        elif grid_mode == "xzy":
            grid_inds = (x_grid*GridSize+z_grid)*GridSize+y_grid # xzy
        elif grid_mode == "yxz":
            grid_inds = (y_grid*GridSize+x_grid)*GridSize+z_grid # yxz
        elif grid_mode == "yzx":
            grid_inds = (y_grid*GridSize+z_grid)*GridSize+x_grid # yzx
        elif grid_mode == "zxy":
            grid_inds = (z_grid*GridSize+x_grid)*GridSize+y_grid # zxy
        elif grid_mode == "zyx":
            grid_inds = (z_grid*GridSize+y_grid)*GridSize+x_grid # zyx
       
        unique_inds, unique_count = np.unique(grid_inds, return_counts=True)
        unique_inds = unique_inds.astype(int)
        grid[unique_inds] += unique_count


    grid_path="./pseudo_density_grid_{0}".format(grid_mode)
    grid.tofile(grid_path) 
    print("Saved to {0}".format(grid_path))

    return grid

def get_density(density_path, galaxy_path, halo_path, GridSize=256,
                precision=2, snap=97, mode=0, grid_mode="xyz"):

    if mode == 0:
        fname = "{0}{1:03d}.dens.dat".format(density_path, snap)
        density = ReadScripts.read_binary_grid(fname, GridSize, precision, False)       
        return density

    elif mode == 1:
        galaxy_density = generate_galaxy_density(galaxy_path, GridSize, snap,
                                                 "flip_nion", num_files=64)
        return galaxy_density       

    elif mode == 2 or mode == 3:
        grid_name = "kali_subfind_{0}_{1}".format(grid_mode, snap)
        halo_density = generate_halo_density(halo_path, grid_name,
                                             GridSize=GridSize, halo_mode=mode, 
                                             snap=snap, grid_mode=grid_mode)

        return halo_density

    elif mode == 4:
        density = generate_pseudo_density(density_path, GridSize, snap,
                                          grid_mode=grid_mode)
        return density


def plot_density(density_grid, GridSize, output_tag, cut_slice=127,
                        cut_thickness=3):

    density_grid = np.reshape(density_grid, (GridSize, GridSize, GridSize))

    fig1 = plt.figure(figsize=(8,8))
    ax1 = fig1.add_subplot(111)

    density_grid = density_grid/np.mean(density_grid)
    #density_grid[density_grid < 1] = 1

    min_dens = np.amin(density_grid)
    max_dens = np.amax(density_grid)

    print("min {0} max {1}".format(min_dens, max_dens))
    #max_dens = 20 

    cmap = plt.cm.get_cmap("afmhot_r")
    #cmap.set_under(color='black')           
    #im = ax1.imshow(np.mean(density_grid[cut_slice:cut_slice+cut_thickness,:,:],axis=0), # Cut yz
    #im = ax1.imshow(np.mean(density_grid[:,cut_slice:cut_slice+cut_thickness,:],axis=1), # Cut xz
                    #vmin=1e-5, vmax=max_dens, 
    im = ax1.imshow(np.mean(density_grid[:,:,cut_slice:cut_slice+cut_thickness],axis=2), #Cut xy
                    interpolation='nearest', origin='low',
                    norm=LogNorm(vmin=0.1, vmax=max_dens), 
                    extent=[0,AllVars.BoxSize,0,AllVars.BoxSize], 
                    cmap=cmap)

    cax = fig1.add_axes()
    cbar = fig1.colorbar(im, cax=cax)
    #cbar.ax.set_ylabel(r'$\mathbf{\delta_{psuedo}}$', rotation = 90, size = PlotScripts.global_labelsize)
    cbar.ax.set_ylabel(r'$\mathbf{\delta_{full}}$', rotation = 90, size = PlotScripts.global_labelsize)

    ax1.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$',
                   fontsize=PlotScripts.global_labelsize)  
    ax1.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$',
                   fontsize=PlotScripts.global_labelsize)  

    ax1.set_xlim([0.0, AllVars.BoxSize]) 
    ax1.set_ylim([0.0, AllVars.BoxSize])

    plt.tight_layout()
    #outputFile = OutputDir + output_tag + '_z' + str(z) + output_format 
    outputFile = "./{0}{1}".format(output_tag, output_format) 
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))
            
    plt.close()

 
if __name__ == '__main__':

    AllVars.Set_Params_Kali()
    PlotScripts.Set_Params_Plot()

    density_path_densfield="/fred/oz004/jseiler/kali/density_fields/1024_subsampled_256/snap"
    density_path_mine="/fred/oz004/jseiler/kali/tmp/snap"
    pseudo_snap_path="/fred/oz004/jseiler/kali/pseudo_snapshots/groups"

    galaxy_path="/fred/oz004/jseiler/kali/self_consistent_output/flip_nion/galaxies/flip_nion"
    #galaxy_path="/fred/oz004/jseiler/kali/self_consistent_output/new_fej/galaxies/fej_alpha0.0_beta0.2"

    #nion_path="/fred/oz004/jseiler/kali/self_consistent_output/flip_nion/grids/nion/flip_nion_fesc0.20_HaloPartCut32_nionHI_097"
    nion_path="/fred/oz004/jseiler/kali/self_consistent_output/new_fej/grids/nion/fej_alpha0.0_beta0.2_ejected_0.000_0.200_HaloPartCut32_nionHI_097"

    #XHII_path="/fred/oz004/jseiler/kali/self_consistent_output/grids/cifog/fej_alpha0.4_beta0.0_XHII_098"
    output_tag="density_galaxy_crosscorr_flip_swap01"
    #output_tag="nion_dens_crosscorr_fiducial_F"

    halo_path = "/fred/oz004/jseiler/kali/shifted_trees/subgroup_trees"
    subfind_halo_path = "/fred/oz004/jseiler/kali/subfind_catalogs/subfind"

    GridSize = 256 
    #density = get_density(density_path, galaxy_path, halo_path, mode=1)
   
    
    dm_density_densfield = get_density(density_path_densfield, galaxy_path, halo_path, mode=0,
                                       snap=92)

    #dm_density_mine = get_density(density_path_mine, galaxy_path, halo_path, mode=0,
    #                              snap=92)
   
    #dm_density_pseudo = get_density(pseudo_snap_path, galaxy_path, halo_path,
    #                                GridSize=GridSize, mode=4, snap=92, 
    #                                grid_mode="xyz")

    #galaxy_density = get_density(density_path, galaxy_path, halo_path, mode=1)
    #subfind_halo_density = get_density("", "", 
    #                                   subfind_halo_path, mode=3, snap=92,
    #                                   grid_mode="zyx", GridSize=GridSize)
   
    halo_density = get_density("", "", 
                               halo_path, mode=2, snap=92,
                               grid_mode="zyx", GridSize=GridSize)
 
    #plot_density(dm_density_pseudo, 256, "density_cuts/dm_pseudo_zyx_snap92")
    #exit()

    '''
    max_dens = np.argmax(density)
    print("Maximum density is {0} and Nion at this spot is "
          "{1}".format(density[max_dens], nion[max_dens])) 

    max_nion = np.argmax(nion)
    print("Maximum Nion is {0} and dens at this spot is "
          "{1}".format(nion[max_nion], density[max_nion])) 
    '''

    dm_density_densfield = np.reshape(dm_density_densfield, (GridSize, GridSize, GridSize))    
    #dm_density_mine = np.reshape(dm_density_mine, (GridSize, GridSize, GridSize))    
    #galaxy_density = np.reshape(galaxy_density, (GridSize, GridSize, GridSize))    
    halo_density = np.reshape(halo_density, (GridSize, GridSize, GridSize))

    #dm_density_pseudo = np.reshape(dm_density_pseudo, (GridSize, GridSize, GridSize))
    #subfind_halo_density = np.reshape(subfind_halo_density, (GridSize, GridSize, GridSize))


    kmid_bins, cross_pspec, pspec1, pspec2 = AllVars.calc_cross_corr(dm_density_densfield,
                                                                     halo_density,
                                                                     AllVars.BoxSize)

    crosscorreff = cross_pspec / (pspec1* pspec2)**0.5
    print(crosscorreff)
    plot_crosscoeff(kmid_bins, crosscorreff,
                    "crosscorrs/shifted_crosscorr_zyx") 
