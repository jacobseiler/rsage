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
import matplotlib.cm as cm

from Corrfunc.theory.xi import xi

from mpl_toolkits.axes_grid1 import make_axes_locatable

from scipy.integrate import quad

output_format = ".png"

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

def plot_fesc_z(halo_bins, fesc_fej, fesc_quasar, output_tag, plot_time=1):

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax3 = fig.add_subplot(212)

    fig2 = plt.figure()
    ax5 = fig2.add_subplot(211)
    ax7 = fig2.add_subplot(212)

    cosmo = cosmology.FlatLambdaCDM(H0 = AllVars.Hubble_h*100, Om0 = AllVars.Omega_m) 
    t_BigBang = cosmo.lookback_time(100000).value # Lookback time to the Big Bang in Gyr.

    max_fej = []
    max_bins_fej = []   

    max_quasar = []
    max_bins_quasar = []   
    
    z = []
    lookback = []
 
    for i in range(20, 101): 

        max_fej_snap = max(fesc_fej[i][5::])        
        max_bins_fej_snap = halo_bins[i][[k for k,j in enumerate(fesc_fej[i]) if j == max_fej_snap][0]] 

        max_fej.append(max_fej_snap)
        max_bins_fej.append(max_bins_fej_snap)       
      
        max_quasar_snap = max(fesc_quasar[i][5::])           
        max_bins_quasar_snap = halo_bins[i][[k for k,j in enumerate(fesc_quasar[i]) if j == max_quasar_snap][0]] 

        if(max_bins_quasar_snap > 11.0): # This is our outlier range. Beyond this there is only a few galaxies in the bin.
            max_quasar_snap = heapq.nlargest(2, (fesc_quasar[i][5::]))[1] 
            max_bins_quasar_snap = halo_bins[i][[k for k,j in enumerate(fesc_quasar[i]) if j == max_quasar_snap][0]] 
 
        max_quasar.append(max_quasar_snap)
        max_bins_quasar.append(max_bins_quasar_snap)       

        z.append(AllVars.SnapZ[i])
        lookback.append((t_BigBang - AllVars.Lookback_Time[i]) * 1.0e3)

        #print("Snap {0} Max_Bin_Quasar {1}".format(i, max_bins_quasar_snap))

    if (plot_time == 1):
        im = ax1.scatter(lookback, max_fej, c = max_bins_fej, cmap = 'viridis') 
        im = ax3.scatter(lookback, max_quasar, c = max_bins_quasar, cmap = 'viridis') 

        ax1.xaxis.set_minor_locator(mtick.MultipleLocator(PlotScripts.time_tickinterval))
        ax1.set_xlim(PlotScripts.time_xlim)
        ax3.xaxis.set_minor_locator(mtick.MultipleLocator(PlotScripts.time_tickinterval))
        ax3.set_xlim(PlotScripts.time_xlim)

        ax1.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_fontsize)   
        ax3.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_fontsize)

        ax2 = ax1.twiny()
        ax4 = ax3.twiny()

        t_plot = (t_BigBang - cosmo.lookback_time(PlotScripts.z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
        z_labels = ["$%d$" % x for x in PlotScripts.z_plot] # Properly Latex-ize the labels.

        ax2.set_xlabel(r"$z$", size = PlotScripts.global_labelsize)
        ax2.set_xlim(PlotScripts.time_xlim)
        ax2.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
        ax2.set_xticklabels(z_labels) # But label them as redshifts.

        ax4.set_xlabel(r"$z$", size = PlotScripts.global_labelsize)
        ax4.set_xlim(PlotScripts.time_xlim)
        ax4.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
        ax4.set_xticklabels(z_labels) # But label them as redshifts.  

        im2 = ax5.scatter(lookback, max_bins_fej, c = max_fej, cmap = 'viridis') 
        im2 = ax7.scatter(lookback, max_bins_quasar, c = max_quasar, cmap = 'viridis') 

        ax5.xaxis.set_minor_locator(mtick.MultipleLocator(PlotScripts.time_tickinterval))
        ax5.set_xlim(PlotScripts.time_xlim)
        ax7.xaxis.set_minor_locator(mtick.MultipleLocator(PlotScripts.time_tickinterval))
        ax7.set_xlim(PlotScripts.time_xlim)

        ax5.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_fontsize)   
        ax7.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_fontsize)

        ax6 = ax5.twiny()
        ax8 = ax7.twiny()

        t_plot = (t_BigBang - cosmo.lookback_time(PlotScripts.z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
        z_labels = ["$%d$" % x for x in PlotScripts.z_plot] # Properly Latex-ize the labels.

        ax6.set_xlabel(r"$z$", size = PlotScripts.global_labelsize)
        ax6.set_xlim(PlotScripts.time_xlim)
        ax6.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
        ax6.set_xticklabels(z_labels) # But label them as redshifts.

        ax8.set_xlabel(r"$z$", size = PlotScripts.global_labelsize)
        ax8.set_xlim(PlotScripts.time_xlim)
        ax8.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
        ax8.set_xticklabels(z_labels) # But label them as redshifts.  

    else:

        im = ax1.scatter(z, max_fej, c = max_bins_fej, cmap = 'viridis') 
        im = ax3.scatter(z, max_quasar, c = max_bins_quasar, cmap = 'viridis') 

        ax1.set_xlabel(r"$z$", size = PlotScripts.global_fontsize)   
        ax3.set_xlabel(r"$z$", size = PlotScripts.global_fontsize) 
 
        ax1.set_xlim([PlotScripts.z_plot[0], PlotScripts.z_plot[-1]])
        ax3.set_xlim([PlotScripts.z_plot[0], PlotScripts.z_plot[-1]])

        ax1.xaxis.set_minor_locator(mtick.MultipleLocator(0.25))
        ax3.xaxis.set_minor_locator(mtick.MultipleLocator(0.25))

        ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.025))
        ax3.yaxis.set_minor_locator(mtick.MultipleLocator(0.025))

        im2 = ax5.scatter(z, max_bins_fej, c = max_fej, cmap = 'viridis') 
        im2 = ax7.scatter(z, max_bins_quasar, c = max_quasar, cmap = 'viridis') 

        ax5.set_xlabel(r"$z$", size = PlotScripts.global_fontsize)   
        ax7.set_xlabel(r"$z$", size = PlotScripts.global_fontsize) 
 
        ax5.set_xlim([PlotScripts.z_plot[0], PlotScripts.z_plot[-1]])
        ax7.set_xlim([PlotScripts.z_plot[0], PlotScripts.z_plot[-1]])

        ax5.xaxis.set_minor_locator(mtick.MultipleLocator(0.25))
        ax7.xaxis.set_minor_locator(mtick.MultipleLocator(0.25))

        ax5.yaxis.set_minor_locator(mtick.MultipleLocator(0.25))
        ax7.yaxis.set_minor_locator(mtick.MultipleLocator(0.25))

    ax1.set_ylabel(r'$\mathrm{Peak} \: f_\mathrm{esc}$', fontsize = PlotScripts.global_fontsize) 
    ax3.set_ylabel(r'$\mathrm{Peak} \: f_\mathrm{esc}$', fontsize = PlotScripts.global_fontsize) 

    ax5.set_ylabel(r'$\log_{10} \mathrm{Peak} \: M_\mathrm{vir} \: [M_\odot]$', fontsize = PlotScripts.global_fontsize) 
    ax7.set_ylabel(r'$\log_{10} \mathrm{Peak} \: M_\mathrm{vir} \: [M_\odot]$', fontsize = PlotScripts.global_fontsize) 

    cbaxes = fig.add_axes([0.78, 0.15, 0.03, 0.8])  
    cb = plt.colorbar(im, cax = cbaxes) 
    cb.set_label(r"$\log_{10} \mathrm{Peak} \: M_\mathrm{vir} \: [M_\odot]$")

    cbaxes = fig2.add_axes([0.78, 0.15, 0.03, 0.8])  
    cb = plt.colorbar(im2, cax = cbaxes) 
    cb.set_label(r"$\mathrm{Peak} \: f_\mathrm{esc}$")

    fig.tight_layout(rect = [0, 0, 0.8, 1])
    fig2.tight_layout(rect = [0, 0, 0.8, 1])
    if (plot_time == 1):
        outputFile1 = './fesc_{0}_Time{1}'.format(output_tag, output_format)
        outputFile2 = './mvir_{0}_Time{1}'.format(output_tag, output_format)
    else:
        outputFile1 = './peak_fesc_{0}_redshift{1}'.format(output_tag, output_format)
        outputFile2 = './peak_mvir_{0}_redshift{1}'.format(output_tag, output_format)

    fig.savefig(outputFile1)  # Save the figure
    fig2.savefig(outputFile2)  # Save the figure

    print("Saved to {0}".format(outputFile1))
    print("Saved to {0}".format(outputFile2))

    plt.close(fig)
    plt.close(fig2)

def determine_autocorr(horizontal_outdir, SnapNum):    

    nbins = 10
    bins = np.linspace(0.1, 10.0, nbins + 1)
     
    Halos = read_horizontal_tree(horizontal_outdir, SnapNum)

    print("Snapshot {0} there is {1} halos".format(SnapNum, len(Halos)))
   
    x = Halos['Pos'][:,0]
    y = Halos['Pos'][:,1]
    z = Halos['Pos'][:,2]
    
    xi_counts = xi(AllVars.BoxSize, 2, bins, Halos['Pos'][:,0], Halos['Pos'][:,1], Halos['Pos'][:,2]) 

    return xi_counts

def plot_fesc_snap(halo_bins, fesc_fej, fesc_quasar, output_tag, snaplist):

    ax1 = plt.subplot(111)

    for snap_idx in range(len(snaplist)):
        snapshot = snaplist[snap_idx]
        redshift = "z = {0:.2f}".format(AllVars.SnapZ[snapshot])
        ax1.plot(halo_bins[snapshot], fesc_fej[snapshot], label = redshift, color = PlotScripts.colors[snap_idx], ls = '-') 
        ax1.plot(halo_bins[snapshot], fesc_quasar[snapshot], color = PlotScripts.colors[snap_idx], ls = '--') 

    ax1.plot(np.nan, np.nan, color = 'k', ls = '-', label = "Ejected Fraction")
    ax1.plot(np.nan, np.nan, color = 'k', ls = '--', label = "Quasar Fraction")

    leg = ax1.legend(loc=1, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')

    ## Output ##

    outputFile = './{0}{1}'.format(output_tag, output_format)
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))

def write_horizontal_tree(tiamat_dir, horizontal_dir, SnapNum):

    snap_halo = []

    for file_idx in range(27):
        Tiamat_Halos, _ = ReadScripts.read_trees_smallarray(tiamat_dir, file_idx, 1) 

        w = np.where(Tiamat_Halos['SnapNum'] == SnapNum)[0]

        print("There are {0} Halos in this file.  Of these, there are {1} Halos at snapshot {2}".format(len(Tiamat_Halos), len(w), SnapNum))
        for i in range(len(w)):
            snap_halo.append(Tiamat_Halos[w[i]])

    fname = "{1}/horizontal_snap{0}".format(SnapNum, horizontal_dir)
    output_file = open(fname, "wb")
   
    output_file.write(np.int32(len(snap_halo)))
    for halo_idx in range(len(snap_halo)):
        for field in names:
            output_file.write(snap_halo[halo_idx][field])

    output_file.close()
 
def read_horizontal_tree(horizontal_dir, SnapNum):

    fname = "{1}/horizontal_snap{0}".format(SnapNum, horizontal_dir)
   
    fin = open(fname, "rb") 
 
    halos_thisfile = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of halos in file.   
    print("Reading {0} Halos".format(halos_thisfile)) 
    Halos = np.fromfile(fin, Halo_Desc, halos_thisfile)  # Read in the halos.

    fin.close()

    return Halos

def calculate_density_field(tiamat_dir, horizontal_outdir, SnapList):

    def get_xi_value(xi, r):

        for r_val in range(len(xi)):
            if (r > xi['rmin'][r_val] and r < xi['rmax'][r_val]):  
                xi_value = xi['xi'][r_val]
                return xi_value
       
        print("Requested an r value of {0}.  This is outside of the binned values for xi.".format(r))
        raise ValueError()
 
    def integrand(xi, k, r):
        # Here xi will be the auto-correlation structure array from autocorr.  
        # Depending upon the value of r, we will call a helper function to determine the correct value of xi(r).

        xi_value = get_xi_value(xi, r)

    for snap in SnapList:
        
        xi = determine_autocorr(horizontal_outdir, snap)
        print(xi['ravg'])
        print(xi['xi'])
        xi_value = get_xi_value(xi, 1.324) 
        exit() 
if __name__ == '__main__':

    tiamat_dir = "/lustre/projects/p070_astro/gpoole/Simulations/Tiamat/trees/version_nominal_res/vertical/"
 
    halo_bins = np.loadtxt("/lustre/projects/p004_swin/jseiler/tiamat_analysis/halo_bins_$f_\mathrm{esc} \: \propto \: f_\mathrm{ej}$.txt")

    fesc_fej = np.loadtxt("/lustre/projects/p004_swin/jseiler/tiamat_analysis/mean_fesc_halo_$f_\mathrm{esc} \: \propto \: f_\mathrm{ej}$.txt")
    fesc_quasar = np.loadtxt("/lustre/projects/p004_swin/jseiler/tiamat_analysis/mean_fesc_halo_$f_\mathrm{esc} \: \propto \: \mathrm{quasar} (\mathrm{boost} = 1.0, N_\mathrm{dynamical} = 1)$.txt")

    horizontal_outdir = "/lustre/projects/p004_swin/jseiler"

    cosmo = AllVars.Set_Params_Tiamat_extended()
    PlotScripts.Set_Params_Plot()
    # [79, 64, 52] = z = [6, 7, 8]

    #write_horizontal_tree(tiamat_dir, horizontal_outdir, 79)  
    #exit() 

    #SnapList = [40, 52, 64, 79, 89, 95]
    SnapList = [79] 
 
    #plot_fesc_z(halo_bins, fesc_fej, fesc_quasar, "fesc_analysis", 0)     
    #plot_fesc_snap(halo_bins, fesc_fej, fesc_quasar, "fesc_snapshot", SnapList)

    calculate_density_field(tiamat_dir, horizontal_outdir, SnapList) 
 
