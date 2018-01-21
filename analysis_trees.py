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

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

group_dir = '/lustre/projects/p004_swin/bsmith/1.6Gpc/means/halo_1721673/dm_gadget/data/'
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

Rockstar_Halo_Desc_full= [
('ID',              np.int64),
('pos',             (np.float32, 6)),
('coreval',         (np.float32, 3)),
('bulkvel',         (np.float32, 3)),
('mass',               np.float32),
('r',               np.float32),
('child_r',         np.float32),
('vmax_r',          np.float32),
('mgrav',           np.float32),
('vmax',            np.float32),
('rvmax',           np.float32),
('rs',              np.float32),
('klypin_rs',       np.float32),
('vrms',            np.float32),
('J',               (np.float32, 3)),
('energy',          np.float32),
('spin',            np.float32),
('alt_m',           (np.float32, 4)),
('Xoff',            np.float32),
('Voff',            np.float32),
('b_to_a',          np.float32),
('c_to_a',          np.float32),
('A',               (np.float32, 3)),
('b_to_a2',         np.float32),
('c_to_a2',         np.float32),
('A2',              (np.float32, 3)),
('bullock_spin',    np.float32),
('kin_to_pot',      np.float32),
('m_pe_b',          np.float32),
('NumPart',         np.int64),
('NumChildPart',    np.int64),
('p_start',         np.int64),
('desc',            np.int64),
('flags',           np.int64),
('n_core',          np.int64),
('min_pos_err',     np.float32),
('min_vel_err',     np.float32),
('min_bulkvel_err', np.float32)
                 ] 

names = [Rockstar_Halo_Desc_full[i][0] for i in range(len(Rockstar_Halo_Desc_full))]
formats = [Rockstar_Halo_Desc_full[i][1] for i in range(len(Rockstar_Halo_Desc_full))]
Rockstar_Halo_Desc = np.dtype({'names':names, 'formats':formats}, align=True)

Rockstar_header_full= [
('magic',           np.int64),
('SnapNum',         np.int64),
('chunk',           np.int64),
('scale',           np.float32),
('Omega_m',         np.float32),
('Omega_l',         np.float32),
('Hubble_h',         np.float32),
('bounds',          (np.float32, 6)),
('NumHalos',        np.int64),
('NumPart',         np.int64), 
('BoxSize',         np.float32),
('PartMass',        np.float32),
('PartType',        np.int64),
('format_revision', np.int32),
('rockstar_version', (np.int8, 12)),
('unused',          (np.int8, 144))
                 ] 
names = [Rockstar_header_full[i][0] for i in range(len(Rockstar_header_full))]
formats = [Rockstar_header_full[i][1] for i in range(len(Rockstar_header_full))]
Rockstar_header = np.dtype({'names':names, 'formats':formats}, align=True)

def get_tree_offsets(treedir, file_idx, simulation):
    """
    This function determines the number of halos per tree. Since the file is read in full this will represent the halo offset for each tree.
    Assumes the tree are named as '<tree_dir>/subgroup_trees_<file_idx>.dat' where file_idx is padded out to 3 digits or '<tree_dir>/lhalotree.bin.<file_idx>' depending on the simulation. 

    Parameters
    ==========

    treedir : string
        Base directory path for the trees.
    file_idx : int
        File number we are reeding in.   
    simulation : int
        Simulation we are reading for.  Determines the naming convention for the files.
        0 : Pip (Britton's Simulation)
        1 : Tiamat
        2 : Manodeep's 1024 Simulation
 
    Returns
    =======

    HalosPerTree : array of ints
        Number of halos within each tree of the file.
    """

    if (simulation == 0 or simulation == 1):
        fname = "{0}/subgroup_trees_{1:03d}.dat".format(treedir, file_idx)
    elif (simulation == 2):
        fname = "{0}/lhalotree.bin.{1}".format(treedir, file_idx)
    else:
        raise ValueError("Invalid simulation option chosen.")
    
    print("Reading for file {0}".format(fname)) 
    fin = open(fname, 'rb')  # Open the file

    trees_thisfile = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file        
    halos_thisfile = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of gals in file.

    HalosPerTree = np.fromfile(fin, np.dtype((np.int32, trees_thisfile)),1)[0] # Read the number of Halos in each tree

    return HalosPerTree 

def read_trees_onearray(treedir, simulation):
    """
    Reads all halos of a simulation into an array. Only recommended for small simulations where memory is not an issue. 
    Assumes the tree are named as '<tree_dir>/subgroup_trees_<file_idx>.dat' where file_idx is padded out to 3 digits or '<tree_dir>/lhalotree.bin.<file_idx>' depending on the simulation. 

    Parameters
    ==========

    treedir : string
        Base directory path for the trees.
    simulation : int
        Simulation we are reading for.  Determines the naming convention for the files.
        0 : Pip (Britton's Simulation)
        1 : Tiamat
        2 : Manodeep's 1024 Simulation

    Returns
    =======

    Halos : array of halos with data-type specified by 'Halo_Desc_full'
        The read in halos for this simulation. 
    
    Notes
    =====

    This function should only be used for small simulations in which memory is not an issue. 
    """ 


    total_halos = 0
    total_trees = 0

    for file_idx in range(num_files):
        if (simulation == 0 or simulation == 1):
            fname = "{0}/subgroup_trees_{1:03d}.dat".format(treedir, file_idx)
        elif (simulation == 2):
            fname = "{0}/lhalotree.bin.{1}".format(treedir, file_idx)
        else:
            raise ValueError("Invalid simulation option chosen.")

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
   

def read_fof_tab_mass(group_dir, snapshot_idx, use_PartMass=0):
    """
    Reads the friends-of-friends HDF5 files of Pip. 
    Assumes the input mass is in units of 1e10 Msun/h.

    Parameters
    ==========

    group_dir : string
        Directory name for the files. Path is assumed to be '<group_dir>/groups_<snapshot_idx>_xx.hdf5' where the snapshot_idx is padded to 3 digits and 'xx' is the subfile number (0 to num_cores). 
    snapshot_idx : int
        The snapshot number we are reading.
    use_partmass : int 
        Toggles whether we want to use the virial mass of the group (use_PartMass = 0) or the number of particles in the group times the particle mass (use_PartMass = 1).
        Default 0.

    Returns
    =======

    mass_fof_tab : array of floats 
        Array containing the mass of each group for this snapshot. 
        Units is the same as the input (usually 1e10 Msun/h). 
    """ 

    mass_fof_tab = []
    for core_idx in range(num_cores):
        fname_fof_tab_ids = "{1}groups_{0:03d}/fof_tab_{0:03d}.{2}.hdf5".format(snapshot_idx, group_dir, core_idx)

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

    return mass_fof_tab

def hmf_fof_tab(snapshot_idx, group_dir, m_low, m_high, use_partMass=0):
    """
    Creates the halo mass function for the friends-of-friends HDF5 files of Pip.
    The width of the bins is given by 'bin_width'.

    Parameters
    ==========

    snapshot_idx : int
        The snapshot number we are reading.
    group_dir : string
        Directory name for the files. Path is assumed to be '<group_dir>/groups_<snapshot_idx>_xx.hdf5' where the snapshot_idx is padded to 3 digits and 'xx' is the subfile number (0 to num_cores). 
    m_low, m_high : floats
        Lower and upper bounds for the halo mass binning.
    use_partmass : int 
        Toggles whether we want to use the virial mass of the group (use_PartMass = 0) or the number of particles in the group times the particle mass (use_PartMass = 1).
        Default 0.

    Returns
    =======

    counts_fof_tab : array of ints 
        Number count of halos within each mass bin. Bins defined by m_low, m_high with bin width given by 'bin_width'. 
    """ 

   
    mass_fof_tab = read_fof_tab_mass(group_dir, snapshot_idx, use_partMass) 
    (counts_fof_tab, bin_edges, bin_middle) = AllVars.Calculate_Histogram(mass_fof_tab, bin_width, 0, m_low, m_high)

    return counts_fof_tab        

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
            Halos_ThisFile, HalosPerTree = ReadScripts.read_trees_smallarray(treedir, file_idx, simulation) 
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


def compare_halo_numbers(snaplow, snaphigh, rockstar_dir, tiamat_dir):
    """
    TODO: Update this function.         
    """ 


    subfind_total = 0
    Rockstar_total = 0
    Tiamat_total = 0

    print("Note that Tiamat and Britton's simulation do not align perfectly with redshift values and we are looking at the SNAPSHOT numbers.")
    for snapshot_idx in range(snaplow, snaphigh + 1):
        fname_subfind_groups = "{0}catalogs/subfind_{1:03d}.catalog_groups_properties/subfind_{1:03d}.catalog_groups_properties.0".format(subfind_dir, snapshot_idx)

        ## SUBFIND ##
        with open(fname_subfind_groups, 'rb') as file_subfind_groups:
            file_number = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)
            n_files = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)
            N_groups_thisfile = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0]
            N_groups_allfile = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)

        subfind_total += N_groups_allfile

        ## Rockstar ##

        num_cores = 64
        for core_idx in range(num_cores):
            fname = "{0}halos_{2}.{1}.ascii".format(rockstar_dir, core_idx, snapshot_idx)
            print("Reading from file {0}".format(fname))

            ID, num_p, mvir, mbound_vir, rvir, vmax, rvmax, vrms, x, y, z, vx, vy, vz, Jx, Jy, Jz, E, Spin, PosUncertainty, VelUncertainty, bulk_vx, bulk_vy, bulk_vz, BulkVelUnc, n_core, m200b, m200c, m500c, m2500c, Xoff, Voff, spin_bullock, b_to_a, c_to_a, A_x, A_y, A_z, b_to_a_500c, c_to_a_500c, A_x_500c, A_y_500c, A_z_500c, Rs, Rs_Klypin, T_on_U, M_pe_Behroozi, M_pe_Diemer, idx, i_so, i_ph, num_cp, mmetric = np.loadtxt(fname, skiprows=20, unpack = True)

            Rockstar_total += len(mvir) 
       
        ## Tiamat ##

        for file_idx in range(27): 
            Halos_Tiamat, dummy = ReadScripts.read_trees_smallarray(tiamat_dir, file_idx, 1)
            w = np.where((Halos_Tiamat['SnapNum'] == snapshot_idx))[0]  
            Tiamat_total += len(w)
            print(len(w))
        
        print("Snapshot \t SUBFIND \t Rockstar \t Tiamat") 
        print("{0:6d} \t {1:6d} \t {2:6d} \t {3:6d}".format(snapshot_idx, int(subfind_total), int(Rockstar_total), int(Tiamat_total))) 

    '''
    subfind_total = int(subfind_total)
    tree_total = int(tree_total)
    print("======================")
    print("========TOTALS========")
    print("======================")
    print("SUBFIND \t Trees \t Difference \t Ratio")
    print("{0:6d} \t {1:6d} \t {2:6d} \t {3:6.4f}".format(subfind_total, tree_total, subfind_total - tree_total, tree_total / subfind_total)) 
    '''
    print("Done")

def check_pointers(Halos, simulation, treedir=None):
    """
    Checks the existence of -1 pointers (for NextProgenitor and FirstProgenitor) for the specified simulation.   

    Parameters
    ==========

    Halos : array with data-type given by 'Halo_Desc_full'
        Array containing the halos we want to check.       

    TODO : Properly update this function. 
    """



    if (simulation == 1):
        Halos = ReadScripts.read_trees_smallarray(treedir, 14, simulation) 

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
  

    print("Snapshot \t Halos at this Snapshot \t NHalos With Non -1 FirstProg at this Snapshot \t Ratio of Halos with non -1 FirstProg at this Snapshot")
    for snapshot_idx in range(20, 91):

        w_firstprog = np.where((Halos['FirstProgenitor'] != -1) & (Halos['SnapNum'] > snapshot_idx))[0]
        w_firstprog_thissnap = np.where((Halos['FirstProgenitor'] != -1) & (Halos['SnapNum'] == snapshot_idx))[0]
        w_snaphigh = np.where((Halos['SnapNum'] > snapshot_idx))[0]  
        w_thissnap = np.where((Halos['SnapNum'] == snapshot_idx))[0]  

        if len(w_thissnap != 0):
            print("{0:7d} \t {1:7d} \t {2:7d} \t {3:7.4f}".format(snapshot_idx, len(w_thissnap), len(w_firstprog_thissnap), len(w_firstprog_thissnap) / len(w_thissnap))) 

def load_rockstar_file(rockstar_dir, snapshot_idx):

    fname = "{0}out_{1}.list".format(rockstar_dir, snapshot_idx)
    print("Reading from file {0}".format(fname))

    ID, Mvir, pos_x, pos_y, pos_z = np.loadtxt(fname, usecols = (0, 2, 8, 9, 10), comments = "#", skiprows = 16, unpack = True)

    return ID, Mvir, pos_x, pos_y, pos_z

def create_converted_subfind_hmf(subfind_dir, snapshot_idx, m_low, m_high):

    fname_subfind_groups = "{0}_{1:03d}.catalog_groups".format(subfind_dir, snapshot_idx)

    with open(fname_subfind_groups, 'rb') as file_subfind_groups:
        N_groups = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0]       
        n_byte_group_offsets = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0]
        grouplen = np.fromfile(file_subfind_groups, np.dtype(np.int32), N_groups)
        group_offsets = np.fromfile(file_subfind_groups, np.dtype(np.int32), N_groups)
        group_sub = np.fromfile(file_subfind_groups, np.dtype(np.int32), N_groups)

        print("N_groups = {0}\t grouplen = {1}".format(N_groups, grouplen))
    
    halo_mass_subfind = np.log10(grouplen * AllVars.PartMass / AllVars.Hubble_h)
        
    (halo_counts_subfind, bin_edges, bin_middle) = AllVars.Calculate_Histogram(halo_mass_subfind, bin_width, 0, m_low, m_high)

    print(halo_counts_subfind)
    return halo_counts_subfind, bin_middle

def create_rockstar_hmf(rockstar_dir, snapshot_idx, m_low, m_high):
    """
    Generates the halo mass function for ROCKSTAR files.

    Parameters
    ==========

    rockstar_dir : int 
        Base path for the location of the ROCKSTAR files. 
    snapshot_idx : int  
        Snapshot number we are creating the function for.
    m_low, m_high : floats
        Lower and upper bounds for the halo mass binning.

    Returns
    =======

    halo_counts : array of ints 
        Number count of halos within each mass bin. Bins defined by m_low, m_high with bin width given by 'bin_width'. 
    bin_middle : array of floats
        Location of the middle of each mass bin.

    Notes
    =====

    Units within the ROCKSTAR files are assumed to be Msun/h. 
    """ 

    ID, Mvir, pos_x, pos_y, pos_z = load_rockstar_file(rockstar_dir, snapshot_idx)

    mvir = np.log10(Mvir)
    (halo_counts, bin_edges, bin_middle) = AllVars.Calculate_Histogram(mvir, bin_width, 0, m_low, m_high)
    
    return halo_counts, bin_middle

def plot_rockstar_only(rockstar_dir, snapshot_idx, label, m_low, m_high, ax):
    """
    Plots the halo mass function for halos found using ROCKSTAR. 

    Parameters
    ==========

    rockstar_dir : int 
        Base path for the location of the ROCKSTAR files. 
    snapshot_idx : int  
        Snapshot number we are creating the function for.
    label : string
        Label to be plotted for this line.
    m_low, m_high : floats
        Lower and upper bounds for the halo mass binning.
    ax : axis
        Axis object that we are plotting on.

    Returns
    =======

    None, the graph is plotted on the passed axis.

    Notes
    =====

    Units plotting on the x/y axis are log10(Msun/h) and Mpc^-3 dex^-1.
    """ 

    halo_counts, bin_middle = create_rockstar_hmf(rockstar_dir, snapshot_idx, m_low, m_high)
    ax.plot(bin_middle, np.multiply(halo_counts / pow(AllVars.BoxSize,3) * pow(AllVars.Hubble_h, 3) / bin_width, bin_middle), label = label)

def plot_converted_subfind_only(subfind_dir, snapshot_idx, label, m_low, m_high, ax):

    halo_counts, bin_middle = create_converted_subfind_hmf(subfind_dir, snapshot_idx, m_low, m_high)
    ax.plot(bin_middle, np.multiply(halo_counts / pow(AllVars.BoxSize,3) * pow(AllVars.Hubble_h, 3) / bin_width, bin_middle), label = label)

def create_pip_tree_hmf(pip_dir, snaplow, snaphigh, m_low, m_high, num_files, use_PartMass=0):

    NB = int((m_high - m_low)/bin_width) 
    halo_counts = np.zeros((snaphigh - snaplow + 1, NB), dtype = np.int64)    

    for file_idx in range(num_files):
       
        Halos_Pip, dummy = ReadScripts.read_trees_smallarray(tiamat_dir, file_idx, 1)

        print("File {0}".format(file_idx))
        count = 0
        for snapshot_idx in range(snaplow, snaphigh + 1):
      
            print("Snapshot {0}".format(snapshot_idx)) 

            cosmol = AllVars.Set_Params_Tiamat_extended()
                                   
            w = np.where((Halos_Tiamat['SnapNum'] == map_pip_to_tiamat[snapshot_idx]))[0] # Find those halos whose redshift is closest to Pip's redshift at this snapshot. 
            if (use_PartMass == 0):
                halo_mass = np.log10(Halos_Tiamat['Mvir'][w] * 1.0e10 / AllVars.Hubble_h)
            else:
                halo_mass = np.log10(Halos_Tiamat['Len'][w] * AllVars.PartMass / AllVars.Hubble_h)
            
            (halo_counts_trees, bin_edges, bin_middle) = AllVars.Calculate_Histogram(halo_mass, bin_width, 0, m_low, m_high)
            halo_counts[count] += halo_counts_trees
            count += 1

    return halo_counts


def create_horizontal_trees_hmf(tiamat_dir, snaplow, snaphigh, map_pip_to_tiamat, m_low, m_high, num_files, simulation, use_PartMass=0):
    """
    Generates the halo mass function from the Tiamat trees. 
    Since we need to read in the entire set of Tiamat files regardless of the snapshot number we call this function once for the entire snapshot range. 

    Parameters
    ==========

    tiamat_dir : int 
        Base path for the location of the Tiamat trees. 
    snaplow, snaphigh : int 
        Bounds for the snapshots we are creating the HMF for.
    map_pip_to_tiamat : array of ints
        Since we often want to plot a comparison between Pip and Tiamat, this array provides the snapshots of Tiamat that most closely match the redshifts of Pip.
        If we want to plot only Tiamat this array should be a range(number Tiamat Snapshots). 
    m_low, m_high : floats
        Lower and upper bounds for the halo mass binning.
    num_files : int
        The number of files in this simulation.
    simulation : int
        Simulation we are reading for.  Determines the naming convention for the files.
        0 : Pip (Britton's Simulation) built using Greg's code
        1 : Tiamat
        2 : Manodeep's 1024 Simulation
        3 : Pip built using Rockstar.
        4 : Kali built using Greg's code.
    use_partmass : int 
        Toggles whether we want to use the virial mass of the group (use_PartMass = 0) or the number of particles in the group times the particle mass (use_PartMass = 1).
        Default 0.

    Returns
    =======

    halo_counts : array of ints 
        Number count of halos within each mass bin. Bins defined by m_low, m_high with bin width given by 'bin_width'. 

    Notes
    =====

    Units within the Tiamat trees are assumed to be 1e10 Msun/h. 
    """ 

    NB = int((m_high - m_low)/bin_width) 
    halo_counts = np.zeros((snaphigh - snaplow + 1, NB), dtype = np.int64)    

    for file_idx in range(num_files):
       
        Halos_Tiamat, dummy = ReadScripts.read_trees_smallarray(tiamat_dir, file_idx, simulation)

        print("File {0}".format(file_idx))
        count = 0
        for snapshot_idx in range(snaplow, snaphigh + 1):
      
            print("Snapshot {0}".format(snapshot_idx)) 

            if (simulation == 1):
                cosmol = AllVars.Set_Params_Tiamat_extended()
            elif (simulation == 3):
                cosmol = AllVars.Set_Params_Britton()
            elif (simulation == 4):
                cosmol = AllVars.Set_Params_Kali()
            else:
                print("Set simulation variable correctly for create horizontal_hmf")
                exit()

            w = np.where((Halos_Tiamat['SnapNum'] == map_pip_to_tiamat[snapshot_idx]))[0] # Find those halos whose redshift is closest to Pip's redshift at this snapshot. 
            if (use_PartMass == 0):
                halo_mass = np.log10(Halos_Tiamat['Mvir'][w] * 1.0e10 / AllVars.Hubble_h)
            else:
                halo_mass = np.log10(Halos_Tiamat['Len'][w] * AllVars.PartMass / AllVars.Hubble_h)
            
            (halo_counts_trees, bin_edges, bin_middle) = AllVars.Calculate_Histogram(halo_mass, bin_width, 0, m_low, m_high)
            halo_counts[count] += halo_counts_trees
            count += 1

    return halo_counts

def plot_horizontal_tree_hmf_only(halo_counts, snapshot_number, label, m_low, m_high, ax): 
    """
    Plots the halo mass function for Tiamat halos. 

    Parameters
    ==========

    halo_counts : Nested 2D array of ints with length equal to the number of snapshots (only one snapshot plotted in this function).
       Number of halos within each mass bin. 
    snapshot_number : int 
       The snapshot number we are plotting for. This specifies which array index of halo_counts we are plotting.          
    label : string
        Label to be plotted for this line.
    m_low, m_high : floats
        Lower and upper bounds for the halo mass binning.
    ax : axis
        Axis object that we are plotting on.

    Returns
    =======

    None, the graph is plotted on the passed axis.

    Notes
    =====

    Units plotting on the x/y axis are log10(Msun/h) and Mpc^-3 dex^-1.
    """ 
 
    bin_edges = np.arange(m_low, m_high + bin_width, bin_width)  
    bin_middle = bin_edges[:-1] + 0.5 * bin_width
        
    ax.plot(bin_middle, np.multiply(halo_counts[snapshot_number] / pow(AllVars.BoxSize,3) * pow(AllVars.Hubble_h, 3) / bin_width, bin_middle), label = label)

def plot_fof_tab_only(group_dir, snapshot_idx, label, m_low, m_high, ax):
    """
    Plots the halo mass function for Pip friends-of-friends HDF5 files. 

    Parameters
    ==========

    group_dir : string
        Directory name for the files. Path is assumed to be '<group_dir>/groups_<snapshot_idx>_xx.hdf5' where the snapshot_idx is padded to 3 digits and 'xx' is the subfile number (0 to num_cores). 
    snapshot_idx : int
        The snapshot number we are reading. 
    label : string
        Label to be plotted for this line.
    m_low, m_high : floats
        Lower and upper bounds for the halo mass binning.
    ax : axis
        Axis object that we are plotting on.

    Returns
    =======

    None, the graph is plotted on the passed axis.

    Notes
    =====

    Units plotting on the x/y axis are log10(Msun/h) and Mpc^-3 dex^-1.
    """ 

    fof_tab_count = hmf_fof_tab(snapshot_idx, group_dir, m_low, m_high) 
    bin_edges = np.arange(m_low, m_high + bin_width, bin_width)  
    bin_middle = bin_edges[:-1] + 0.5 * bin_width

    ax.plot(bin_middle, np.multiply(fof_tab_count / pow(AllVars.BoxSize,3) * pow(AllVars.Hubble_h, 3) / bin_width, bin_middle), label = label, color = 'r')
  
def create_hmf_nbodykit_fof(fof_dir, snapshot_idx, m_low, m_high):

    fname = "{0}groups_{1:03d}/nbodykit_fof_tab_{1:03d}.hdf5".format(fof_dir, snapshot_idx)

    with h5py.File(fname, "r") as f:
        mass = f['Mass'][:]

    Mvir = np.log10(mass * 1.0e10 / AllVars.Hubble_h)
    (halo_counts, bin_edges, bin_middle) = AllVars.Calculate_Histogram(Mvir, bin_width, 0, m_low, m_high)

    return halo_counts, bin_middle
 
def plot_nbodykit_fof_only(fof_dir, snapshot_idx, label, m_low, m_high, ax):

    fof_nbodykit_count, bin_middle = create_hmf_nbodykit_fof(fof_dir, snapshot_idx, m_low, m_high)

    ax.plot(bin_middle, np.multiply(fof_nbodykit_count / pow(AllVars.BoxSize,3) * pow(AllVars.Hubble_h, 3) / bin_width, bin_middle), label = label, color = 'r')

def plot_thibault_hmf(ax):

    bin_middle, phi = np.loadtxt("/home/jseiler/hmf_z9.5.dat", unpack = True, skiprows = 2)  
    ax.plot(bin_middle, phi, label = "Thibault z = 9.5")   
   

def create_subfind_hmf(subfind_dir, snapshot_idx, m_low, m_high):

    fname_subfind_groups = "{0}/catalogs/subfind_{1:03d}.catalog_groups_properties/subfind_{1:03d}.catalog_groups_properties.0".format(subfind_dir, snapshot_idx)

    with open(fname_subfind_groups, 'rb') as file_subfind_groups:
        i_file = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0]       
        n_file= np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0]
        N_groups_thisfile= np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0]
        N_groups_total = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0] 

    print("Subfind properties are split up over {0} files".format(n_file))
    Mvir = []
    for file_idx in range(n_file):
        
        fname = "{0}/catalogs/subfind_{1:03d}.catalog_groups_properties/subfind_{1:03d}.catalog_groups_properties.{2}".format(subfind_dir, snapshot_idx, file_idx)
        with open(fname, 'rb') as file_subfind_groups:
            i_file = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0]       
            n_file= np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0]
            N_groups_thisfile= np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0]
            N_groups_total = np.fromfile(file_subfind_groups, np.dtype(np.int32), 1)[0] 

            print("File {0} has {1} groups".format(file_idx, N_groups_thisfile))
            Halo_thisfile = np.fromfile(file_subfind_groups, SUBFIND_Halo_Desc, N_groups_thisfile) 
            
            for i in range(N_groups_thisfile):
                Mvir.append(np.log10(Halo_thisfile[i]['M_vir'] / AllVars.Hubble_h))     
        
    (halo_counts_subfind, bin_edges, bin_middle) = AllVars.Calculate_Histogram(Mvir, bin_width, 0, m_low, m_high)

    print(halo_counts_subfind)
    return halo_counts_subfind, bin_middle
    
 
def plot_kali_subfind_only(kali_subfind_dir, snapshot_idx, label, m_low, m_high, ax):
     
    halo_counts, bin_middle = create_subfind_hmf(kali_subfind_dir, snapshot_idx, m_low, m_high)
    ax.plot(bin_middle, np.multiply(halo_counts / pow(AllVars.BoxSize,3) * pow(AllVars.Hubble_h, 3) / bin_width, bin_middle), label = label)
      
def plot_rockstar_tiamat_comparison(rockstar_dir, tiamat_dir, group_dir, subfind_dir, rockstar_tree_dir, nbodykit_fof_dir, kali_subfind_dir, kali_dir, snaplow, snaphigh, SnapZ_Pip, SnapZ_Tiamat, SnapZ_Kali):
    """
    Plots the halo mass function over a snapshot range (1 for each snapshot)  for the halos of Tiamat and Pip found by a) ROCKSTAR and b) the friends-of-friends HDF5 files on a single graph. 
    Since the redshifts of Pip and Tiamat are not identical, we plot the snapshots of the closest redshifts.
    The snapshot range is the range of Pip snapshots and the saved graph will have redshift name using the Pip redshift. 

    Parameters
    ==========

    rockstar_dir, tiamat_dir, group_dir : strings
        Directory names for the files. See the documentation for the individual 'create_xxx_hmf' functions for the assumed format of the file names.  
    snaplow, snaphigh : int 
        Bounds for the snapshots we are creating the HMF for.       
    SnapZ_Pip, SnapZ_Tiamat : array of floats
        Redshift for each snapshot of Pip and Tiamat.

    Returns
    =======

    None, the graph is saved with name specified by 'outputFile' 

    Notes
    =====

    Units plotting on the x/y axis are log10(Msun/h) and Mpc^-3 dex^-1.
    """ 
    m_low = 6
    m_high = 12

    map_pip_to_tiamat = AllVars.find_nearest_redshifts(SnapZ_Pip, SnapZ_Tiamat)       
    map_kali_to_tiamat = AllVars.find_nearest_redshifts(SnapZ_Kali, SnapZ_Tiamat)   
 
    tiamat_halo_counts = create_horizontal_trees_hmf(tiamat_dir, snaplow, snaphigh, map_kali_to_tiamat, m_low, m_high, 27, 1)
    #pip_halo_counts = create_horizontal_trees_hmf(pip_dir, snaplow, snaphigh, np.arange(0, 92), m_low, m_high, 64, 3)
    kali_halo_counts = create_horizontal_trees_hmf(kali_dir, snaplow, snaphigh, np.arange(0, 106), m_low, m_high, 64, 4)        

    count = 0
    for snapshot_idx in range(snaplow, snaphigh + 1):
        
        ax1 = plt.subplot(111)
        
        cosmo = AllVars.Set_Params_Kali()
        redshift = AllVars.SnapZ[snapshot_idx]
        label = "SUBFIND Kali Tab z = {0:.2f}".format(redshift)
        plot_kali_subfind_only(kali_subfind_dir, snapshot_idx, label, m_low, m_high, ax1)


        label = "Kali Tree Tab z = {0:.2f}".format(redshift)
        plot_horizontal_tree_hmf_only(kali_halo_counts, count, label, m_low, m_high, ax1)
   
        #cosmol = AllVars.Set_Params_Britton()
        #redshift = AllVars.SnapZ[snapshot_idx]
        #label = "FoF Fraction = 0.7 z = {0:.3f}".format(AllVars.SnapZ[-7])
        #plot_rockstar_only(rockstar_dir, len(SnapZ_Pip) - 7, label, m_low, m_high, ax1)
        #ax1.axvline(9, color = 'k', ls = '--')

        #label = "FoF Tab z = {0:.3f}".format(AllVars.SnapZ[snapshot_idx])
        #plot_fof_tab_only(group_dir, snapshot_idx, label, m_low, m_high, ax1)
    
        #plot_thibault_hmf(ax1)
         
        #label = "Nbodykit FoF z = {0:.3f}".format(AllVars.SnapZ[snapshot_idx])
        #plot_nbodykit_fof_only(nbodykit_fof_dir, snapshot_idx, label, m_low, m_high, ax1)
 
        #cosmol = AllVars.Set_Params_Britton()
        #redshift = AllVars.SnapZ[snapshot_idx]
        #label = "SUBFIND Britton z = {0:.3f}".format(AllVars.SnapZ[snapshot_idx])
        #plot_converted_subfind_only(subfind_dir, snapshot_idx, label, m_low, m_high, ax1)
        #
        cosmol = AllVars.Set_Params_Tiamat_extended()
        label = "Tiamat z = {0:.3f}".format(AllVars.SnapZ[map_kali_to_tiamat[snapshot_idx]])        
        plot_horizontal_tree_hmf_only(tiamat_halo_counts, count, label, m_low, m_high, ax1)

        #cosmol = AllVars.Set_Params_Britton()
        #label = "Pip Trees z = {0:.3f}".format(AllVars.SnapZ[map_pip_to_tiamat[snapshot_idx]])        
        #plot_horizontal_tree_hmf_only(pip_halo_counts, count, label, m_low, m_high, ax1)

        ax1.set_xlabel(r'$\log_{10}\ M_{H} \:[M_{\odot}]$', fontsize = PlotScripts.global_fontsize)
        ax1.set_ylabel(r'$\left(\frac{dn}{d\log{M}}\right) \ [\mathrm{Mpc}^{-3}]$', fontsize = PlotScripts.global_fontsize)
        ax1.set_ylim([1e-4, 1e3])
        ax1.set_yscale('log', nonposy='clip')

        leg = ax1.legend(loc='upper right', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(PlotScripts.global_legendsize)

        outputFile = "./HMF_Kali_z{0:.3f}.png".format(redshift)
        #outputFile = "./HMF_Kali_z5.78.png"
        plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print('Saved file to {0}'.format(outputFile))
        plt.close()

        count += 1

def read_rockstar_halos_binary_small(rockstar_dir, snapshot_idx, core_idx):

    fname = "{0}halos_{1}.{2}.bin".format(rockstar_dir, snapshot_idx, core_idx)
    print("Reading file {0}".format(fname))
    fin = open(fname, 'rb')
    header = np.fromfile(fname, Rockstar_header, 1)
    print(header)
    
    print("There is {0} halos for this file".format(header['NumHalos'][0]))
   
    Halos = np.fromfile(fname, Rockstar_Halo_Desc, header['NumHalos'][0])
    
    return Halos 

   
def load_rockstar_particles(rockstar_dir, snapshot_idx, core_idx):

    fname = "{0}halos_{1}.{2}.particles".format(rockstar_dir, snapshot_idx, core_idx)
    print("reading from file {0}".format(fname))
    x, y, z, particle_id, halo_id = np.loadtxt(fname, comments = "#", usecols = (0, 1, 2, 6, 7), unpack = True)

    return x, y, z, particle_id, halo_id
 
def check_rockstar_halos(rockstar_dir, snaplow, snaphigh):
       
    for snapshot_idx in range(snaplow, snaphigh + 1): 

        ID, Mvir, pos_x, pos_y, pos_z = load_rockstar_file(rockstar_dir, snapshot_idx) 

        print(ID)
        exit()
        unique_halos, unique_halos_idx = np.unique(pos_x, return_index = True)
        idx = np.arange(0, len(ID))
        duplicate_halos_idx = [int(x) for x in idx if x not in unique_halos_idx]
        
        print(duplicate_halos_idx) 
        for dup_idx in range(1, len(duplicate_halos_idx)):
            if (math.isclose(pos_x[int(duplicate_halos_idx[0])], pos_x[int(duplicate_halos_idx[dup_idx])], abs_tol = 0.01) == True):
                print(pos_x[duplicate_halos_idx[dup_idx]])
                print(pos_x[duplicate_halos_idx[0]])

        print("For snapshot {0} there are {1} (out of {2}) duplicate halos".format(snapshot_idx, len(duplicate_halos_idx), len(pos_x)))
  
        ax1 = plt.subplot(111)
        w = np.where((np.log10(Mvir / AllVars.Hubble_h) > 8.2))[0]
        print("There is {0} halos above Mass of 8.2^10 Msun".format(len(w)))

        ax1.plot(pos_x[w], pos_y[w], alpha = 0.5) 

        outputFile = "./halos_check_{0:.3f}.png".format(AllVars.SnapZ[snapshot_idx])
        plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print('Saved file to {0}'.format(outputFile))
        plt.close()

def check_rockstar_particles(rockstar_dir, snaplow, snaphigh):

    for snapshot_idx in range(snaplow, snaphigh + 1): 

        x, y, z, particle_id, halo_id = load_rockstar_particles(rockstar_dir, snapshot_idx, 0)
        print("Read in {0} particles".format(len(x)))
      
        unique_part = np.unique(particle_id)
 
        print("For snapshot{0} there are {1} repeat particles (out of {2})".format(snapshot_idx, len(particle_id) - len(unique_part), len(particle_id)))

if __name__ == '__main__':

    if (len(sys.argv) != 4):
        print("Usage: python3 analysis_trees.py <snaplow> <snaphigh> <0 for Britton's Sim, 1 for Tiamat Extended, 2 for Manodeep's Sim>")
        exit()
     
    snaplow = int(sys.argv[1])
    snaphigh = int(sys.argv[2])
    simulation = int(sys.argv[3])

    rockstar_dir = "/lustre/projects/p004_swin/jseiler/Rockstar_output/Britton_Halos_final_noLL/"
    tiamat_dir = "/lustre/projects/p070_astro/gpoole/Simulations/Tiamat/trees/version_nominal_res/vertical/"
    kali_subfind_dir = "/lustre/projects/p134_swin/jseiler/kali/subfind_results"
    kali_dir = "/lustre/projects/p134_swin/jseiler/kali/subfind_results/trees/trees/vertical/"
    subfind_dir = "/lustre/projects/p004_swin/jseiler/rockstar_to_subfind/rockstar"
    pip_dir = "/lustre/projects/p004_swin/jseiler/Rockstar_output/Britton_Halos_final_noLL/Lhalos/"
    nbodykit_fof_dir = "/lustre/projects/p134_swin/jseiler/nbodykit_fofs/"

    cosmol = AllVars.Set_Params_Britton()
    SnapZ_Pip = AllVars.SnapZ
    cosmol = AllVars.Set_Params_Tiamat_extended()
    SnapZ_Tiamat = AllVars.SnapZ
    cosmol = AllVars.Set_Params_Kali()
    SnapZ_Kali = AllVars.SnapZ
  
   
    PlotScripts.Set_Params_Plot()

    #compare_halo_numbers(snaplow, snaphigh, rockstar_dir, tiamat_dir)

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

    plot_rockstar_tiamat_comparison(rockstar_dir, tiamat_dir, group_dir, subfind_dir, pip_dir, nbodykit_fof_dir, kali_subfind_dir, kali_dir, snaplow, snaphigh, SnapZ_Pip, SnapZ_Tiamat, SnapZ_Kali)
    #check_rockstar_halos(rockstar_dir, snaplow, snaphigh)
    #check_rockstar_particles("/lustre/projects/p004_swin/jseiler/Rockstar_output/Britton_checking_dup/", snaplow, snaphigh)
    
