#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')

import os
import heapq
import h5py as h5
import numpy as np
import pylab as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from numpy import *
from random import sample, seed, randint
from os.path import getsize as getFileSize
import math
import random
import csv
from cycler import cycler
from io import StringIO
#np.set_printoptions(threshold=np.nan)
from collections import Counter
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import AxesGrid
from astropy import units as u
from astropy import cosmology

import matplotlib.ticker as mtick
import PlotScripts
import ReadScripts
import AllVars

from mpi4py import MPI
from tqdm import tqdm

from hmf import MassFunction
from hmf import cosmo
from astropy.cosmology import FlatLambdaCDM

import h5py 

groupdir = '/lustre/projects/p004_swin/bsmith/1.6Gpc/means/halo_1721673/dm_gadget/data/'
snapdir = '/lustre/projects/p134_swin/jseiler/simulations/1.6Gpc/means/halo_1721673/dm_gadget/data/'
treedir = "/lustre/projects/p134_swin/jseiler/subfind_britton/trees/britton/vertical"
num_files = 125
num_cores = 256
bin_width = 0.1


Halo_Desc_full = [
('Descendant',          np.int32),
('FirstProgenitor',     np.int32),
('NextProgenitor',      np.int32),
('FirstHaloInFOFgroup', np.int32),
('NextHaloInFOFgroup',  np.int32),
('Len',                 np.int32),
('M_mean200',           np.float32),
('Mvir',                np.float32),
('M_TopHat',            np.float32),
('Pos',                 (np.float32, 3)),
('Vel',                 (np.float32, 3)),
('VelDisp',             np.float32),
('Vmax',                np.float32),
('Spin',                (np.float32, 3)),
('MostBoundID',         np.int64),
('SnapNum',             np.int32),
('Filenr',              np.int32),
('SubHaloIndex',        np.int32),
('SubHalfMass',         np.float32)
 ]

names = [Halo_Desc_full[i][0] for i in range(len(Halo_Desc_full))]
formats = [Halo_Desc_full[i][1] for i in range(len(Halo_Desc_full))]
Halo_Desc = np.dtype({'names':names, 'formats':formats}, align=True)
    
def plot_tree_hmf():

    cosmol = AllVars.Set_Params_Britton()

    total_halos = 0
    total_trees = 0   
 
    for file_idx in range(num_files):
        fname = "{0}/subgroup_trees_{1:03d}.dat".format(treedir, file_idx)

        fin = open(fname, 'rb')  # Open the file

        trees_thisfile = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file        
        halos_thisfile = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of gals in file.

        total_halos += halos_thisfile
        total_trees += trees_thisfile
	
        fin.close()
   
    print("There are {0} total halos".format(total_halos)) 
    # Initialize the storage array
    Halos = np.empty(total_halos, dtype=Halo_Desc)
    offset = 0   

    for file_idx in range(num_files):

        fname = "{0}/subgroup_trees_{1:03d}.dat".format(treedir, file_idx)
        print("Reading for file {0}".format(fname)) 
        fin = open(fname, 'rb')  # Open the file

        trees_thisfile = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file        
        halos_thisfile = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of gals in file.

        GalsPerTree = np.fromfile(fin, np.dtype((np.int32, trees_thisfile)),1) # Read the number of gals in each tree

        Halos_tmp = np.fromfile(fin, Halo_Desc, halos_thisfile)  # Read in the galaxy structures


        Halos[offset:offset+halos_thisfile]=Halos_tmp[0:halos_thisfile].copy()

        del(Halos_tmp)
        offset = offset + halos_thisfile # Update the offset position for the global array

        fin.close()  # Close the file       
    
    print min(Halos['Len'])
    
    count = 0
    for i in range(len(Halos)):
        if (i % 100000 == 0):
            print i
        if Halos['Len'][i] == 1:
            count += 1
    print count
    print len(Halos['Len'] == 1)
    exit()   
    for snap_idx in range(20, 92):

        redshift = AllVars.SnapZ[snap_idx]
        # HMF from the trees ##

        w = np.where((Halos['SnapNum'] == snap_idx))[0]
        halo_mass = np.log10(Halos['Mvir'][w] * 1.0e10 / AllVars.Hubble_h)

        (halo_counts_trees, bin_edges, bin_middle) = AllVars.Calculate_Histogram(halo_mass, bin_width, 0, 6, 12)

        ## Theoretical HMF from hmf ## 

        my_cosmo = cosmo.Cosmology(cosmo_model = cosmol) # Update the hmf cosmology.
        britton_cosmo = FlatLambdaCDM(H0 = 69.5, Om0 = 0.285, Ob0 = 0.04845)
        hmf = MassFunction()
        hmf.update(cosmo_params = {"H0" : 69.5, "Om0" : 0.285}, Mmax = 11, Mmin = 6, z = redshift)

        massfunc = hmf.dndlog10m
        hmf_bins = np.linspace(6.0, 11.0, num = (11.0 - 6.0) / 0.01)

        # Reading the FoF Tab #

        mass_fof_tab = []
        npart_fof_tab = []
        for core_idx in range(num_cores):
            fname_fof_tab_ids = "{1}groups_{0:03d}/fof_tab_{0:03d}.{2}.hdf5".format(snap_idx, groupdir, core_idx)

            with h5py.File(fname_fof_tab_ids, "r") as file_fof_tab:
                try:
                    Ngroups_foftab = file_fof_tab['Group']['GroupMass'].shape[0]
                except KeyError:
                    pass
                else:
                    for group_idx in range(Ngroups_foftab):
                        mass_fof_tab.append(np.log10(file_fof_tab['Group']['GroupMass'][group_idx] * 1.0e10))
        

        (counts_fof_tab, bin_edges, bin_middle) = AllVars.Calculate_Histogram(mass_fof_tab, bin_width, 0, 6, 12)


        # Plotting #
        ax1 = plt.subplot(111)

        title = "z = {0:.3f}".format(redshift)
        ax1.set_title(title)

        label = "Trees"
        ax1.plot(bin_middle, np.multiply(halo_counts_trees / pow(AllVars.BoxSize,3) * pow(AllVars.Hubble_h, 3) / bin_width, bin_middle), label = label)

        label = "HMF"
        ax1.plot(hmf_bins - np.log10(AllVars.Hubble_h), massfunc * pow(AllVars.Hubble_h,3), label = label)

        label = "FoF Tab"
        ax1.plot(bin_middle, np.multiply(counts_fof_tab / pow(AllVars.BoxSize,3) * pow(AllVars.Hubble_h, 3) / bin_width, bin_middle), label = label)

        ax1.set_ylabel(r'$\left(\frac{dn}{d\log{M}}\right) \ [\mathrm{Mpc}^{-3}]$', fontsize = PlotScripts.global_fontsize)
        ax1.set_yscale('log', nonposy='clip')
    #    plt.axis([6, 11.5, 1e-6, 5e0])

        leg = ax1.legend(loc='upper right', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(PlotScripts.global_legendsize)

        outputFile = "./HMF_Trees/TreeHMF_z{0:.3f}.png".format(AllVars.SnapZ[snap_idx])
        plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print('Saved file to {0}'.format(outputFile))
        plt.close()


 
if __name__ == '__main__':

    PlotScripts.Set_Params_Plot()
    plot_tree_hmf()    

