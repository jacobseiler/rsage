#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')

import sys
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
subfind_dir = '/lustre/projects/p134_swin/jseiler/subfind_britton/'
snapdir = '/lustre/projects/p134_swin/jseiler/simulations/1.6Gpc/means/halo_1721673/dm_gadget/data/'
num_files = 27
num_cores = 256
bin_width = 0.1
N_sub_files = 16

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

SUBFIND_Halo_Desc_full = [
('id_MBP',              np.int64),
('M_vir',               np.float64),
('n_particles',         np.int16),
('position_COM',        (np.float32, 3)),
('position_MBP',        (np.float32, 3)),
('velocity_COM',        (np.float32, 3)),
('velocity_MBP',        (np.float32, 3)),
('R_vir',               np.float32),
('R_halo',              np.float32),
('R_max',               np.float32),
('V_max',               np.float32),
('sigma_v',             np.float32),
('spin',                (np.float32, 3)),
('q_triaxial',          np.float32),
('s_triaxial',          np.float32),
('shape_eigen_vectors', (np.float32, (3,3))),
('padding',             (np.int16, 2))
                 ] # Note that there are also a padding of 8 bytes following this array. 

names = [SUBFIND_Halo_Desc_full[i][0] for i in range(len(SUBFIND_Halo_Desc_full))]
formats = [SUBFIND_Halo_Desc_full[i][1] for i in range(len(SUBFIND_Halo_Desc_full))]
SUBFIND_Halo_Desc = np.dtype({'names':names, 'formats':formats}, align=True)

def get_tree_offsets(treedir, file_idx):

    fname = "{0}/subgroup_trees_{1:03d}.dat".format(treedir, file_idx)
    print("Reading for file {0}".format(fname)) 
    fin = open(fname, 'rb')  # Open the file

    trees_thisfile = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file        
    halos_thisfile = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of gals in file.

    HalosPerTree = np.fromfile(fin, np.dtype((np.int32, trees_thisfile)),1)[0] # Read the number of Halos in each tree

    return HalosPerTree 

def read_trees_smallarray(treedir, file_idx, simulation):

    if (simulation == 0 or simulation == 1):
        fname = "{0}/subgroup_trees_{1:03d}.dat".format(treedir, file_idx)
    else:
        fname = "{0}/lhalotree.bin.{1}".format(treedir, file_idx)

    print("Reading for file {0}".format(fname)) 
    fin = open(fname, 'rb')  # Open the file

    trees_thisfile = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file.
    halos_thisfile = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of halos in file.

    HalosPerTree = np.fromfile(fin, np.dtype((np.int32, trees_thisfile)),1)[0] # Read the number of halos in each tree.

    Halos = np.fromfile(fin, Halo_Desc, halos_thisfile)  # Read in the halos.

    fin.close() # Close the file  

    return Halos, HalosPerTree

def read_trees_onearray(treedir):

    total_halos = 0
    total_trees = 0

    for file_idx in range(num_files):
        fname = "{0}/subgroup_trees_{1:03d}.dat".format(treedir, file_idx)

        fin = open(fname, 'rb')  # Open the file

        trees_thisfile = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file        
        halos_thisfile = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of halos in file.

        total_halos += halos_thisfile
        total_trees += trees_thisfile
	
        fin.close()
   
    print("There are {0} total halos".format(total_halos)) 
 
    Halos = np.empty(total_halos, dtype=Halo_Desc)
    offset = 0   

    for file_idx in range(num_files):

        fname = "{0}/subgroup_trees_{1:03d}.dat".format(treedir, file_idx)
        print("Reading for file {0}".format(fname)) 
        fin = open(fname, 'rb')  # Open the file

        trees_thisfile = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file        
        halos_thisfile = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of gals in file.

        HalosPerTree = np.fromfile(fin, np.dtype((np.int32, trees_thisfile)),1) # Read the number of Halos in each tree

        Halos_tmp = np.fromfile(fin, Halo_Desc, halos_thisfile)  # Read in the galaxy structures


        Halos[offset:offset+halos_thisfile]=Halos_tmp[0:halos_thisfile].copy()

        del(Halos_tmp)
        offset = offset + halos_thisfile # Update the offset position for the global array

    fin.close() # Close the file  

    return Halos 
   

def plot_hmf_comparison(snaplow, snaphigh, treesdir):

    m_low = 6
    m_high = 12

    num_files = 27
    use_PartMass = 0 
    
    NB = int((m_high - m_low)/bin_width) 
    tree_count_Britton = np.zeros((snaphigh - snaplow + 1, NB), dtype = np.int64)
    tree_count_Tiamat = np.zeros((snaphigh - snaplow + 1, NB), dtype = np.int64)

    for file_idx in range(num_files):
       
        Halos_Britton, dummy = read_trees_smallarray(treesdir[0], file_idx, 0)
        Halos_Tiamat, dummy = read_trees_smallarray(treesdir[1], file_idx, 1)

        print("File {0}".format(file_idx))
        count = 0
        for snapshot_idx in range(snaplow, snaphigh + 1):
      
            print("Snapshot {0}".format(snapshot_idx)) 
            # First do Britton's Sim #

            cosmol = AllVars.Set_Params_Britton()
            w = np.where((Halos_Britton['SnapNum'] == snapshot_idx))[0]
            if (use_PartMass == 0):
                halo_mass = np.log10(Halos_Britton['Mvir'][w] * 1.0e10 / AllVars.Hubble_h)
            else:
                halo_mass = np.log10(Halos_Britton['Len'][w] * AllVars.PartMass / AllVars.Hubble_h)

            (halo_counts_trees, bin_edges, bin_middle) = AllVars.Calculate_Histogram(halo_mass, bin_width, 0, m_low, m_high)
            tree_count_Britton[count] += halo_counts_trees   
 
            cosmol = AllVars.Set_Params_Tiamat_extended()
            w = np.where((Halos_Tiamat['SnapNum'] == snapshot_idx))[0]
            if (use_PartMass == 0):
                halo_mass = np.log10(Halos_Tiamat['Mvir'][w] * 1.0e10 / AllVars.Hubble_h)
            else:
                halo_mass = np.log10(Halos_Tiamat['Len'][w] * AllVars.PartMass / AllVars.Hubble_h)

            (halo_counts_trees, bin_edges, bin_middle) = AllVars.Calculate_Histogram(halo_mass, bin_width, 0, m_low, m_high)
            tree_count_Tiamat[count] += halo_counts_trees
            count += 1
   
    count = 0
    for snapshot_idx in range(snaplow, snaphigh + 1): 
        ax1 = plt.subplot(111) 
        # First Britton "

        cosmol = AllVars.Set_Params_Britton()
        label = "Britton z = {0:.3f}".format(AllVars.SnapZ[snapshot_idx])
        ax1.plot(bin_middle, np.multiply(tree_count_Britton[count] / pow(AllVars.BoxSize,3) * pow(AllVars.Hubble_h, 3) / bin_width, bin_middle), label = label, color = 'r')
        ax1.axvline(np.log10(AllVars.PartMass * 32 / AllVars.Hubble_h), ymin = -1, ymax = 1e5, color = 'r', linewidth = PlotScripts.global_linewidth, linestyle = '--') 


        # Then Tiamat # 

        cosmol = AllVars.Set_Params_Tiamat_extended()
        label = "Tiamat z = {0:.3f}".format(AllVars.SnapZ[snapshot_idx])
        ax1.plot(bin_middle, np.multiply(tree_count_Tiamat[count] / pow(AllVars.BoxSize,3) * pow(AllVars.Hubble_h, 3) / bin_width, bin_middle), label = label, color = 'b')
        ax1.axvline(np.log10(AllVars.PartMass * 32 / AllVars.Hubble_h), ymin = -1, ymax = 1e5, color = 'b', linewidth = PlotScripts.global_linewidth, linestyle = '--') 

        ## 
        count += 1
        ax1.set_xlabel(r'$\log_{10}\ M_{H} \:[M_{\odot}]$', fontsize = PlotScripts.global_fontsize)
        ax1.set_ylabel(r'$\left(\frac{dn}{d\log{M}}\right) \ [\mathrm{Mpc}^{-3}]$', fontsize = PlotScripts.global_fontsize)
        ax1.set_ylim([1e-4, 1e2])
        ax1.set_yscale('log', nonposy='clip')

        simulation_tag = "Comparison"

        leg = ax1.legend(loc='upper right', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(PlotScripts.global_legendsize)

        if (use_PartMass == 0):
            outputFile = "./HMF_Trees/{1}_DivideLittleh_TreeHMF_z{0:.3f}.png".format(AllVars.SnapZ[snapshot_idx], simulation_tag)
        else:
            outputFile = "./HMF_Trees/{1}_DivideLittleh_npart_TreeHMF_z{0:.3f}.png".format(AllVars.SnapZ[snapshot_idx], simulation_tag)
        plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print('Saved file to {0}'.format(outputFile))
        plt.close()
def plot_tree_hmf(snaplow, snaphigh, cosmol, treedir, simulation, compare_values=0, use_PartMass=0):

    m_low = 6
    m_high = 12

    total_halos = 0
    total_trees = 0   

    if (use_PartMass == 0):
        print("USING MVIR")
    else:
        print("USING NPART * PARTMASS")

    if (simulation == 0 or simulation == 1):
        num_files = 27
    else:
        num_files = 125
 
    for file_idx in range(num_files):
        if (simulation == 0 or simulation == 1):
            fname = "{0}/subgroup_trees_{1:03d}.dat".format(treedir, file_idx)
        else:
            fname = "{0}/lhalotree.bin.{1}".format(treedir, file_idx)
        fin = open(fname, 'rb')  # Open the file

        trees_thisfile = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file        
        halos_thisfile = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of gals in file.

        total_halos += halos_thisfile
        total_trees += trees_thisfile
	
        fin.close()
   
    print("There are {0} total halos".format(total_halos)) 
  
    NB = int((m_high - m_low)/bin_width) 
    halo_counts_trees_total = np.zeros((snaphigh - snaplow + 1, NB), dtype = np.int64)
       
    for file_idx in range(num_files):

        if (simulation == 0 or simulation == 1):
            fname = "{0}/subgroup_trees_{1:03d}.dat".format(treedir, file_idx)
        else:
            fname = "{0}/lhalotree.bin.{1}".format(treedir, file_idx)
        
        print("Reading for file {0}".format(fname)) 
        fin = open(fname, 'rb')  # Open the file

        trees_thisfile = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file        
        halos_thisfile = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of gals in file.

        GalsPerTree = np.fromfile(fin, np.dtype((np.int32, trees_thisfile)),1) # Read the number of gals in each tree

        Halos = np.fromfile(fin, Halo_Desc, halos_thisfile)  # Read in the galaxy structures

        fin.close()  # Close the file       
       
        count = 0 
        for snapshot_idx in range(snaplow, snaphigh + 1):

            # HMF from the trees ##

            w = np.where((Halos['SnapNum'] == snapshot_idx))[0]
            if (use_PartMass == 0):
                halo_mass = np.log10(Halos['Mvir'][w] * 1.0e10 / AllVars.Hubble_h) 
            else:
                halo_mass = np.log10(Halos['Len'][w] * AllVars.PartMass / AllVars.Hubble_h) 

            (halo_counts_trees, bin_edges, bin_middle) = AllVars.Calculate_Histogram(halo_mass, bin_width, 0, m_low, m_high)
            halo_counts_trees_total[count] += halo_counts_trees
            count += 1

    count = 0
    if (simulation == 0):
        print("BRITTON SIM")
    elif (simulation == 1):
        print("TIAMAT")
    else:
        print("MANODEEP SIM")

    for snapshot_idx in range(snaplow, snaphigh+ 1):       
        if (compare_values == 1): 
            print("==================================================")
            print("===========SNAPSHOT {0}========================".format(snapshot_idx))
            print("Mass Range \t Trees Value \t HMF Value \t HMF Value / Trees Value")
        ## Theoretical HMF from hmf ## 
        redshift = AllVars.SnapZ[snapshot_idx]
        my_cosmo = cosmo.Cosmology(cosmo_model = cosmol) # Update the hmf cosmology.
        britton_cosmo = FlatLambdaCDM(H0 = AllVars.Hubble_h * 100.0, Om0 = AllVars.Omega_m, Ob0 = 0.06) # Ob0 is a dummy variable so just hard code it. 
        hmf = MassFunction()
        hmf.update(cosmo_params = {"H0" : AllVars.Hubble_h * 100.0, "Om0" : AllVars.Omega_m}, Mmax = m_high, Mmin = m_low, z = redshift, dlog10m = bin_width)
       
        massfunc = hmf.dndlog10m
        hmf_bins = np.linspace(m_low, m_high, num = (m_high - m_low) / bin_width) 

        # Reading the FoF Tab #

        if (simulation == 0):
            mass_fof_tab = []
            npart_fof_tab = []
            for core_idx in range(num_cores):
                fname_fof_tab_ids = "{1}groups_{0:03d}/fof_tab_{0:03d}.{2}.hdf5".format(snapshot_idx, groupdir, core_idx)

                with h5py.File(fname_fof_tab_ids, "r") as file_fof_tab:
                    try:
                        Ngroups_foftab = file_fof_tab['Group']['GroupLen'].shape[0]
                    except KeyError:
                        pass
                    else:
                        for group_idx in range(Ngroups_foftab):
                            if (use_PartMass == 0):
                                mass_fof_tab.append(np.log10(file_fof_tab['Group']['GroupMass'][group_idx] * 1.0e10 / AllVars.Hubble_h))
                            else:
                                mass_fof_tab.append(np.log10(file_fof_tab['Group']['GroupLen'][group_idx] * AllVars.PartMass / AllVars.Hubble_h))
            

            (counts_fof_tab, bin_edges, bin_middle) = AllVars.Calculate_Histogram(mass_fof_tab, bin_width, 0, 6, 12)

        # Reading the SUBFIND files # 
            fname_subfind_groups = "{0}catalogs/subfind_{1:03d}.catalog_groups_properties/subfind_{1:03d}.catalog_groups_properties.0".format(subfind_dir, snapshot_idx)

            with open(fname_subfind_groups, 'rb') as file_subfind_groups:
                file_number = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)
                n_files = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0]
                N_groups_thisfile = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0]
                N_groups_allfile = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)

            print("Snapshot {0} has the SUBFIND results split across {1} files".format(snapshot_idx, n_files))
            for i_file in range(n_files):
                fname_subfind_groups = "{0}catalogs/subfind_{1:03d}.catalog_groups_properties/subfind_{1:03d}.catalog_groups_properties.{2}".format(subfind_dir, snapshot_idx, i_file)

#                print("Reading from file {0}".format(fname_subfind_groups))
                with open(fname_subfind_groups, 'rb') as file_subfind_groups:
                    file_number = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)
                    n_files = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)
                    N_groups_thisfile = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0]
                    N_groups_allfile = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)

                    #print("This file contains {0} groups".format(N_groups_thisfile))                    
                    SUBFIND_Halos = np.fromfile(file_subfind_groups, SUBFIND_Halo_Desc, N_groups_thisfile)  # Read in the galaxy structures

                    if (use_PartMass == 0):
                        halo_mass_subfind = np.log10(SUBFIND_Halos['M_vir'] / AllVars.Hubble_h)
                    else:
                        halo_mass_subfind = np.log10(SUBFIND_Halos['n_particles'] * AllVars.PartMass / AllVars.Hubble_h)

                    if (i_file == 0):
                        (halo_counts_subfind, bin_edges, bin_middle) = AllVars.Calculate_Histogram(halo_mass_subfind, bin_width, 0, m_low, m_high)
                    else:
                        (halo_counts_subfind_tmp, bin_edges, bin_middle) = AllVars.Calculate_Histogram(halo_mass_subfind, bin_width, 0, m_low, m_high)
                        halo_counts_subfind += halo_counts_subfind_tmp
            print("The total number of SUBFIND groups for Snapshot {0} is {1}".format(snapshot_idx, N_groups_allfile))
            print("The total number of Groups from the trees for snapshot {0} is {1}".format(snapshot_idx, np.sum(halo_counts_trees_total[count]))) 
        # Plotting #
        ax1 = plt.subplot(111)

        title = "z = {0:.3f}".format(redshift)
        ax1.set_title(title)

        label = "Trees"
        ax1.plot(bin_middle, np.multiply(halo_counts_trees_total[count] / pow(AllVars.BoxSize,3) * pow(AllVars.Hubble_h, 3) / bin_width, bin_middle), label = label)

        label = "HMF"
        ax1.plot(hmf_bins, massfunc * pow(AllVars.Hubble_h,3), label = label)
        
        if (compare_values == 1): 
            for bin_idx in range(len(bin_middle)):
                mine = halo_counts_trees_total[count][bin_idx] / pow(AllVars.BoxSize,3) * pow(AllVars.Hubble_h, 3) / bin_width * bin_middle[bin_idx]
                theory = massfunc[bin_idx] * pow(AllVars.Hubble_h,3)
                print("{0:7.4f} - {1:7.4f} \t {2:7.4e} \t {3:7.4e} \t {4:7.4e}".format(bin_middle[bin_idx] - bin_width / 2.0, bin_middle[bin_idx] + bin_width / 2.0, mine, theory, theory / mine)) 

        if (simulation == 0):
            label = "FoF Tab"
            ax1.plot(bin_middle, np.multiply(counts_fof_tab / pow(AllVars.BoxSize,3) * pow(AllVars.Hubble_h, 3) / bin_width, bin_middle), label = label)

            label = "Subfind"
            ax1.plot(bin_middle, np.multiply(halo_counts_subfind / pow(AllVars.BoxSize,3) * pow(AllVars.Hubble_h, 3) / bin_width, bin_middle), label = label)

        ax1.axvline(np.log10(AllVars.PartMass * 32 / AllVars.Hubble_h), ymin = -1, ymax = 1e5, color = 'k', linewidth = PlotScripts.global_linewidth, linestyle = '--') 
        ax1.set_xlabel(r'$\log_{10}\ M_{H} \:[M_{\odot}]$', fontsize = PlotScripts.global_fontsize)
        ax1.set_ylabel(r'$\left(\frac{dn}{d\log{M}}\right) \ [\mathrm{Mpc}^{-3}]$', fontsize = PlotScripts.global_fontsize)
        ax1.set_yscale('log', nonposy='clip')
    #    plt.axis([6, 11.5, 1e-6, 5e0])

        if (simulation == 0):
            simulation_tag = "Britton_Shift"
        elif (simulation == 1):   
            simulation_tag = "Tiamat"
        else:
            simulation_tag = "Manodeep"

        leg = ax1.legend(loc='upper right', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(PlotScripts.global_legendsize)

        if (use_PartMass == 0):
            outputFile = "./HMF_Trees/{1}_DivideLittleh_TreeHMF_z{0:.3f}.png".format(AllVars.SnapZ[snapshot_idx], simulation_tag)
        else:
            outputFile = "./HMF_Trees/{1}_DivideLittleh_npart_TreeHMF_z{0:.3f}.png".format(AllVars.SnapZ[snapshot_idx], simulation_tag)
        plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print('Saved file to {0}'.format(outputFile))
        plt.close()

        if (simulation == 0):
            ax2 = plt.subplot(111)
            ax2.set_title(title)
            ax2.plot(bin_middle, np.divide(halo_counts_trees_total[count], halo_counts_subfind))

            ax2.set_xlabel(r'$\log_{10}\ M_{H} \:[M_{\odot}]$', fontsize = PlotScripts.global_fontsize)
            ax2.set_ylabel(r'$N_{Trees} / N_{SUBFIND}$', fontsize = PlotScripts.global_fontsize)
       
            ax2.set_ylim([0, 1])

            ax3 = ax2.twinx()
            ax3.plot(bin_middle, np.subtract(halo_counts_subfind, halo_counts_trees_total[count]), color = 'r')
            ax3.set_ylabel(r"$\log_{10} N_{SUBFIND} - N_{Trees}$", fontsize = PlotScripts.global_fontsize)
            ax3.tick_params('y', colors='r')

            ax3.set_yscale('log', nonposy='clip')


            plt.tight_layout() 
            if (use_PartMass == 0):
                outputFile = "./HMF_Trees/{1}_DivideLittleh_Ratio_TreeHMF_z{0:.3f}.png".format(AllVars.SnapZ[snapshot_idx], simulation_tag)
            else:
                outputFile = "./HMF_Trees/{1}_DivideLittleh_npart_Ratio_TreeHMF_z{0:.3f}.png".format(AllVars.SnapZ[snapshot_idx], simulation_tag)
            plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
            print('Saved file to {0}'.format(outputFile))
            plt.close()
        
        count += 1
 

def plot_progline(Halos, simulation, treedir=None):

    ## First I want to find a progenitor line that has an appreciable length. ##

    max_progline_global_trees = []

    max_snap = len(AllVars.SnapZ) - 1 
    max_snap = 91

    file_range_low = [0, 6, 12, 18, 23] 
    file_range_high = [5, 11, 17, 22, num_files -1] 

    #for file_idx in range(num_files):
    for meta_file_idx in range(len(file_range_low)):
        for file_idx in range(file_range_low[meta_file_idx], file_range_high[meta_file_idx] + 1): 
            Halos_ThisFile, HalosPerTree = read_trees_smallarray(treedir, file_idx, simulation) 
            ## I now have all the halos from this file, along with the number of halos for each tree for this file. ##

            print("This file has {0} Halos spread across {1} trees.".format(len(Halos_ThisFile), len(HalosPerTree)))
            offset = 0
            progline_global_tree_count = []
            tree_num = []

            w_test = np.where((Halos_ThisFile['SnapNum'] == max_snap)) 
            
    #        print("There are {0} halos at snapshot {2} and of those, the maximum FistProgenitor pointer is {1} and the maximum Descendant pointer is {3}".format(len(w_test), max(Halos_ThisFile['FirstProgenitor'][w_test]), max_snap, max(Halos_ThisFile['Descendant'][w_test]))) 

            for tree_idx in range(len(HalosPerTree)): # Now lets loop over each individual tree within this file.
       
                Halos_ThisTree = Halos_ThisFile[offset:offset+int(HalosPerTree[tree_idx])] # Get those Halos that are in this tree.
                offset += HalosPerTree[tree_idx] # Update the offset.
                w = np.where((Halos_ThisTree['SnapNum']== max_snap))[0] # Get those halos that are still alive near the final snapshot. 
           
                for progline_idx in range(len(w)): # Then lets start to climb the progenitor line of the halos (if any) that exist at snapshot 91. 
                    progline_count = 0 
                    current_halo = w[progline_idx]
                    while (Halos_ThisTree['FirstProgenitor'][current_halo] != -1): # Keep climbing the tree until we reach the top.
                        progline_count += 1 
                        #print("Halo {0} \t FirstProg {1} \t Current SnapNum {2} \t NextProg[NextProg] {3} \t NextProg SnapNum {4}".format(current_halo, Halos_ThisTree['FirstProgenitor'][current_halo], Halos_ThisTree['SnapNum'][current_halo], Halos_ThisTree['FirstProgenitor'][Halos_ThisTree['FirstProgenitor'][current_halo]], Halos_ThisTree['SnapNum'][Halos_ThisTree['FirstProgenitor'][current_halo]])) 
                        current_halo = Halos_ThisTree['FirstProgenitor'][current_halo] # Move to the next halo.
                    progline_global_tree_count.append(progline_count)
                    tree_num.append(tree_idx)
                        
            longest_halo_tree = tree_num[np.argmax(progline_global_tree_count)]
            print("For file {0} the halo with the longest prog line (that is alive at snapshot {3}) is at tree {1} with {2} links".format(file_idx, longest_halo_tree, max(progline_global_tree_count), max_snap)) 

            ## Now we have the longest progenitor line, let's plot some statistics for it. ##

            ax1 = plt.subplot(121)
            ax2 = plt.subplot(122)

            offset_progline = np.sum(HalosPerTree[0:longest_halo_tree])
            Halos_ProgLine = Halos_ThisFile[offset_progline:offset_progline+int(HalosPerTree[longest_halo_tree])] # This is the halos that for the tree with the longest progenitor line within the file.

            w = np.where((Halos_ProgLine['SnapNum']== max_snap))[0] # Get those halos that are still alive near the final snapshot. 

            count_thistree = []
            # Lets just do a stupid thing to get the halo w index that starts the longest progenitor chain
            for progline_idx in range(len(w)): # Then lets start to climb the progenitor line of the halos (if any) that exist at snapshot 91. 
                progline_count = 0 
                current_halo = w[progline_idx]
                while (Halos_ProgLine['FirstProgenitor'][current_halo] != -1): # Keep climbing the tree until we reach the top.
                    progline_count += 1 
                    current_halo = Halos_ProgLine['FirstProgenitor'][current_halo] # Move to the next halo.
                count_thistree.append(progline_count) 
            
            Mvir_ProgLine = []
            Spin_ProgLine = []
            z_ProgLine = []
            current_halo = w[np.argmax(count_thistree)] # This will be the index of w that had the longest progenitor line.
            while (Halos_ProgLine['FirstProgenitor'][current_halo] != -1): # Keep climbing the tree until we reach the top.
                Mvir_ProgLine.append(np.log10(Halos_ProgLine['Mvir'][current_halo] * 1.0e10 / AllVars.Hubble_h))
                Spin_Mag = pow(pow(Halos_ProgLine['Spin'][current_halo][0], 2) + pow(Halos_ProgLine['Spin'][current_halo][1], 2) + pow(Halos_ProgLine['Spin'][current_halo][2], 2), 0.5)
                Spin_ProgLine.append(Spin_Mag)
                z_ProgLine.append(AllVars.SnapZ[Halos_ProgLine['SnapNum'][current_halo]])
                
                current_halo = Halos_ProgLine['FirstProgenitor'][current_halo] # Move to the next halo.
            label = "Halo {0}, File {1}".format(offset_progline + w[progline_idx], file_idx)

            ax1.plot(z_ProgLine, Mvir_ProgLine, label = label, color = 'k', lw = 0.8) 
            ax2.plot(z_ProgLine, Spin_ProgLine, label = label, color = 'k', lw = 0.8) 

            ax1.set_xlabel("z")
            ax2.set_xlabel("z")

            ax1.set_ylabel(r'$\log_{10}\ M_{H} \:[M_{\odot}]$', fontsize = PlotScripts.global_fontsize)
            ax2.set_ylabel(r'$|\Lambda|$')

            ax1.set_xlim([min(z_ProgLine) - 0.1, max(z_ProgLine) + 0.1])
            ax2.set_xlim([min(z_ProgLine) - 0.1, max(z_ProgLine) + 0.1])
           
            ax1.set_ylim([8.0, 11.0])
            ax2.set_ylim([0.0, 0.03]) 

     
            '''
            leg = ax1.legend(loc='upper right', numpoints=1, labelspacing=0.1)
            leg.draw_frame(False)  # Don't want a box frame
            for t in leg.get_texts():  # Reduce the size of the text
                t.set_fontsize(PlotScripts.global_legendsize)
            '''
        plt.tight_layout()

        if (simulation == 0):
            simulation_tag = "Britton_Shift"
        elif (simulation == 1):   
            simulation_tag = "Tiamat"
        else:
            simulation_tag = "Manodeep"

        outputFile = "./MainProgLine_{2}_DivideLittleh_files_{0}_{1}.png".format(file_range_low[meta_file_idx], file_range_high[meta_file_idx], simulation_tag)
        plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print('Saved file to {0}'.format(outputFile))
        plt.close()


def compare_halo_numbers(Halos):

    print("Snapshot \t SUBFIND \t Trees \t Difference \t Ratio")

    subfind_total = 0
    tree_total = 0

    for snapshot_idx in range(0, 92):
        fname_subfind_groups = "{0}catalogs/subfind_{1:03d}.catalog_groups_properties/subfind_{1:03d}.catalog_groups_properties.0".format(subfind_dir, snapshot_idx)

        with open(fname_subfind_groups, 'rb') as file_subfind_groups:
            file_number = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)
            n_files = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)
            N_groups_thisfile = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0]
            N_groups_allfile = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)

        w = np.where((Halos['SnapNum'] == snapshot_idx))[0]

        subfind_total += N_groups_allfile
        tree_total += len(w)

        print("{0:6d} \t {1:6d} \t {2:6d} \t {3:6d} \t {4:6.4f}".format(snapshot_idx, N_groups_allfile[0], len(w), N_groups_allfile[0] - len(w), len(w) / N_groups_allfile[0])) 

    subfind_total = int(subfind_total)
    tree_total = int(tree_total)
    print("======================")
    print("========TOTALS========")
    print("======================")
    print("SUBFIND \t Trees \t Difference \t Ratio")
    print("{0:6d} \t {1:6d} \t {2:6d} \t {3:6.4f}".format(subfind_total, tree_total, subfind_total - tree_total, tree_total / subfind_total)) 


def check_pointers(Halos, simulation, treedir=None):

    if (simulation == 1):
        Halos = read_trees_smallarray(treedir, 14, simulation) 

    w_firstprog = np.where((Halos['FirstProgenitor'] != -1))[0]
    w_nextprog = np.where((Halos['NextProgenitor'] != -1))[0]
    Nhalos = len(Halos)

    print("There are {0} halos with non -1 FirstProgenitor pointers and {1} with non -1 NextProgenitor pointers".format(len(w_firstprog), len(w_nextprog)))
    print("This is out of a total of {0} halos".format(Nhalos))
    print("This means the ratio of halos without -1 FirstProgenitor and NextProgenitor pointers are {0:.4f} and {1:.4f} respectively.".format(len(w_firstprog) / Nhalos, len(w_nextprog) / Nhalos))

    print("==========================")
    print("========SUBSETTING========")
    print("==========================")

    print("Snapshot \t NHalos above Snapshot \t Ratio of TotalHalos above Snapshot \t NHalos With Non -1 NextProg Above Snapshot \t Ratio of Halos with Non -1 NextProg above Snapshot") 
    
    for snapshot_idx in range(0, 92):

        w_firstprog = np.where((Halos['FirstProgenitor'] != -1) & (Halos['SnapNum'] > snapshot_idx))[0]
        w_firstprog_thissnap = np.where((Halos['FirstProgenitor'] != -1) & (Halos['SnapNum'] == snapshot_idx))[0]
        w_snaphigh = np.where((Halos['SnapNum'] > snapshot_idx))[0]  
        w_thissnap = np.where((Halos['SnapNum'] == snapshot_idx))[0]  
        

        print("{0:7d} \t {1:7d} \t {2:7.4f} \t {3:7d} \t {4:7.4f}".format(snapshot_idx, len(w_snaphigh), len(w_snaphigh) / Nhalos, len(w_firstprog), len(w_firstprog) / len(w_snaphigh)))
  
    exit() 
    print("Snapshot \t Halos at this Snapshot \t NHalos With Non -1 FirstProg at this Snapshot \t Ratio of Halos with non -1 FirstProg at this Snapshot")
    for snapshot_idx in range(20, 91):

        w_firstprog = np.where((Halos['FirstProgenitor'] != -1) & (Halos['SnapNum'] > snapshot_idx))[0]
        w_firstprog_thissnap = np.where((Halos['FirstProgenitor'] != -1) & (Halos['SnapNum'] == snapshot_idx))[0]
        w_snaphigh = np.where((Halos['SnapNum'] > snapshot_idx))[0]  
        w_thissnap = np.where((Halos['SnapNum'] == snapshot_idx))[0]  

        if len(w_thissnap != 0):
            print("{0:7d} \t {1:7d} \t {2:7d} \t {3:7.4f}".format(snapshot_idx, len(w_thissnap), len(w_firstprog_thissnap), len(w_firstprog_thissnap) / len(w_thissnap))) 
if __name__ == '__main__':

    if (len(sys.argv) != 4):
        print("Usage: python3 analysis_trees.py <snaplow> <snaphigh> <0 for Britton Sim, 1 for Tiamat Extended>")
        exit()
     
    snaplow = int(sys.argv[1])
    snaphigh = int(sys.argv[2])
    simulation = int(sys.argv[3])

    print("Running with snaplow = {0}, snaphigh = {1}, Simulation = {2}".format(snaplow, snaphigh, simulation))

    if (simulation == 0):
        treedir = "/lustre/projects/p134_swin/jseiler/subfind_britton/trees/britton_shifted/vertical"
        cosmol = AllVars.Set_Params_Britton()
    elif (simulation == 1):
        treedir = "/lustre/projects/p070_astro/gpoole/Simulations/Tiamat/trees/version_nominal_res/vertical/"
        cosmol = AllVars.Set_Params_Tiamat_extended()
    elif (simulation == 2):
        treedir = "/lustre/projects/p004_swin/jseiler/Rockstar_output/1024_Halos_noLL/Ltrees/"
        cosmol = AllVars.Set_Params_Mysim() 
    else:
        print("Only simulation 0 (Britton Sim), 1 (Tiamat extended) and 2 (Mandeep's Simulation) are supported.")
        exit()

    PlotScripts.Set_Params_Plot()

#    if (simulation == 0):
#        Halos = read_trees_onearray(treedir)

    #compare_halo_numbers(Halos)

   
    ''' 
    if (simulation == 0):
        check_pointers(Halos, simulation) 
    else:
        check_pointers(0, simulation, treedir)
    
    if (simulation == 0):
        plot_progline(0, simulation, treedir)
    else:
        plot_progline(0, simulation, treedir)
    '''
    
    plot_hmf_comparison(snaplow, snaphigh, ["/lustre/projects/p134_swin/jseiler/subfind_britton/trees/britton_shifted/vertical", "/lustre/projects/p070_astro/gpoole/Simulations/Tiamat/trees/version_nominal_res/vertical/"]) 

#    plot_tree_hmf(snaplow, snaphigh, cosmol, treedir, simulation)    

