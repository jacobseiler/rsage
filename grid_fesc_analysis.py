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

m_low = 7.0 # We only sum the photons coming from halos within the mass range m_low < Halo Mass < m_high
m_high = 15.0

output_format = ".pdf"

def plot_fesc_z(filename_base, PlotSnapshot, model_tags, output_tag):

    ax1 = plt.subplot(111)

    for model_number in range(len(filename_base)):
        count = 0
        for snapshot_idx in PlotSnapshot[model_number]:
            fname = "{0}_{1}.txt".format(filename_base[model_number], snapshot_idx)
                                                                              
            snap, fesc, Mvir, Ngamma, Ngammafesc = np.loadtxt(fname, unpack=True)

            fname_out = "{0}_{1}".format(filename_base[model_number], snapshot_idx)
            exit()

            Mvir = np.log10(Mvir * 1.0e10 / AllVars.Hubble_h)                   
            mean, std, N, bins_mid = AllVars.Calculate_2D_Mean(Mvir, fesc, 0.2, m_low, m_high) 

            w = np.where((N > 0.0))
            if model_number == 0:
                label = "z = {0:.2f}".format(AllVars.SnapZ[snapshot_idx])
            else:
                label = ""
            ax1.plot(bins_mid[w], mean[w], color = PlotScripts.colors[count], ls = PlotScripts.linestyles[model_number], label = label, lw = PlotScripts.global_linewidth)
            #ax1.fill_between(bins_mid[w], np.subtract(mean[w],std[w]), np.add(mean[w],std[w]), color = PlotScripts.colors[count], alpha = 0.25)
            count += 1
            print("Snapshot {0}, model {1} plotted".format(snapshot_idx, model_number)) 

        #ax1.plot(np.nan, np.nan, color = 'k', ls = PlotScripts.linestyles[model_number], lw = PlotScripts.global_linewidth, label = model_tags[model_number])


    ax1.axhline(0.30, 0, 100, color ='k', linewidth = PlotScripts.global_linewidth, linestyle = '-.')
    ax1.axvline(np.log10(32.0*AllVars.PartMass / AllVars.Hubble_h), color = 'k', linewidth = PlotScripts.global_linewidth, linestyle = '-.')   
    ax1.text(10.7, 0.37, r"$f_\mathrm{esc} = 0.35$", color = 'k', size = PlotScripts.global_fontsize)

    ax1.set_xlabel(r'$\log_{10}\ M_\mathrm{vir}\ [M_{\odot}]$', size = PlotScripts.global_fontsize) 
    ax1.set_ylabel(r'$f_\mathrm{esc}$', size = PlotScripts.global_fontsize)
    ax1.set_xlim([8.6, 11.75])
    ax1.set_ylim([-0.05, 1.0])   

    ax1.set_xticks(np.arange(9.0, 12.0))  
    ax1.xaxis.set_minor_locator(mtick.MultipleLocator(0.25))
    ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.05))
   
    leg = ax1.legend(loc='upper left', bbox_to_anchor=(0.3, 1.02), numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize)

 
    outputFile = './{0}{1}'.format(output_tag, output_format)
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()

def plot_total_nion(SnapList, base_grid_name, GridSize, simulation_norm, plot_time, labels, OutputDir, output_tag):

    print("")
    print("Plotting the total number of ionizing photons") 

    ax1 = plt.subplot(111)

    nion = np.empty((len(total_nion), len(total_nion[0])))

    for model_number in range(0, len(total_nion)):
        if simulation_norm[model_number] == 0:
            cosmo = AllVars.Set_Params_Mysim()
        elif simulation_norm[model_number] == 2:
            cosmo = AllVars.Set_Params_Tiamat()
        elif simulation_norm[model_number] == 6:
            cosmo = AllVars.Set_Params_Kali()
        else:
            raise ValueError("Set simulation norm option in plot_total_nion")
        for snapshot_idx in SnapList: 
            nion_path = "{0}_{1:03d}".format(base_grid_name[model_number], snapshot_idx) 
            nion_snap = ReadScripts.read_binary_grid(nion_path, GridSize[model_number], 1) * 1.0e50
            
            nion[model_number][snapshot_idx] = np.log10(np.sum(nion_snap) / ((AllVars.BoxSize / AllVars.Hubble_h)**3))	

        print("The sum of ionizing photons for model {0} is {1}".format(model_number, sum(nion[model_number]))) 
        if plot_time == 1:
            ax1.plot((AllVars.t_BigBang- cosmo.lookback_time(ZZ).value) * 1.0e3, nion[model_number], color = PlotScripts.colors[model_number], label = labels[model_number], ls = PlotScripts.linestyles[model_number], lw = 3)
        else:	
            ax1.plot(ZZ, nion[model_number], color = PlotScripts.colors[model_number], label = labels[model_number], ls = PlotScripts.linestyles[model_number], lw = 3)

    bouwens_z = np.arange(6,16) # Redshift range for the observations.
    bouwens_t = (AllVars.t_BigBang - cosmo.lookback_time(bouwens_z).value) * 1.0e3 # Corresponding values for what we will plot on the x-axis.

    bouwens_1sigma_lower = [50.81, 50.73, 50.60, 50.41, 50.21, 50.00, 49.80, 49.60, 49.39, 49.18] # 68% Confidence Intervals for the ionizing emissitivity from Bouwens 2015.
    bouwens_1sigma_upper = [51.04, 50.85, 50.71, 50.62, 50.56, 50.49, 50.43, 50.36, 50.29, 50.23]

    bouwens_2sigma_lower = [50.72, 50.69, 50.52, 50.27, 50.01, 49.75, 49.51, 49.24, 48.99, 48.74] # 95% CI. 
    bouwens_2sigma_upper = [51.11, 50.90, 50.74, 50.69, 50.66, 50.64, 50.61, 50.59, 50.57, 50.55]

    if plot_time == 1:
        ax1.fill_between(bouwens_t, bouwens_1sigma_lower, bouwens_1sigma_upper, color = 'k', alpha = 0.2)
        ax1.fill_between(bouwens_t, bouwens_2sigma_lower, bouwens_2sigma_upper, color = 'k', alpha = 0.4, label = r"$\mathrm{Bouwens \: et \: al. \: (2015)}$")
    else:
        ax1.fill_between(bouwens_z, bouwens_1sigma_lower, bouwens_1sigma_upper, color = 'k', alpha = 0.2)
        ax1.fill_between(bouwens_z, bouwens_2sigma_lower, bouwens_2sigma_upper, color = 'k', alpha = 0.4, label = r"$\mathrm{Bouwens \: et \: al. \: (2015)}$")

    if plot_time == 1:
        ax1.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = label_size)
        ax1.xaxis.set_minor_locator(mtick.MultipleLocator(PlotScripts.time_tickinterval))
        ax1.set_xlim(PlotScripts.time_xlim)

        ax1.set_ylabel(r'$\sum f_\mathrm{esc}\dot{N}_\gamma \: [\mathrm{s}^{-1}\mathrm{Mpc}^{-3}]$', fontsize = PlotScripts.global_fontsize)                 
        #ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        ax1.set_ylim([48.5, 51.5])

        ax2 = ax1.twiny()

        t_plot = (AllVars.t_BigBang - cosmo.lookback_time(PlotScripts.z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
        z_labels = ["$%d$" % x for x in PlotScripts.z_plot] # Properly Latex-ize the labels.

        print(z_labels)
        ax2.set_xlabel(r"$z$", fontsize = label_size)
        ax2.set_xlim(PlotScripts.time_xlim)
        ax2.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
        ax2.set_xticklabels(z_labels) # But label them as redshifts.
        
    else:
        ax1.set_xlabel(r"$z$", size = label_size)
        ax1.set_ylabel(r"$\mathrm{log}_{10} \: \dot{N}_{\mathrm{HI}} \: [\mathrm{s}^{-1}\mathrm{Mpc}^{-3}]$", fontsize = label_size)
        ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        ax1.set_ylim([48.5, 51.5])

        '''
        ax3.set_xlim(time_xlim)
        ax3.set_ylim([-0.3, 0.30])
        ax3.set_yticks(np.arange(-0.20, 0.30, 0.20))
        ax3.yaxis.set_minor_locator(mtick.MultipleLocator(0.05))


        for tick in ax3.yaxis.get_major_ticks():
        tick.label.set_fontsize(label_size - 5)
        '''
    ax1.text(350, 50.0, r"$68\%$", horizontalalignment='center', verticalalignment = 'center', fontsize = PlotScripts.global_labelsize)
    ax1.text(350, 50.8, r"$95\%$", horizontalalignment='center', verticalalignment = 'center', fontsize = PlotScripts.global_labelsize)
    leg = ax1.legend(loc='lower right', numpoints=1,
                     labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(legend_size + 3)

    plt.tight_layout()

    outputFile = OutputDir + output_tag + output_format 			
    plt.savefig(outputFile)  # Save the figure	
    print('Saved file to {0}'.format(outputFile))		
        
    plt.close()	

def grid_output_to_binary(SnapList, filename_base, filename_out):

    for snapshot in SnapList:

        fname = "{0}_{1}.txt".format(filename_base[0], snapshot)        
        print("Opening {0}".format(fname))
        snap, fesc, Mvir, Ngamma, Ngamma_fesc = np.loadtxt(fname, unpack=True)
               
        fname_out = "{0}snapshot_{1}".format(filename_out, snapshot)
        np.savez(fname_out, Mvir=Mvir, Ngamma=Ngamma)
        
        name = "{0}.npz".format(fname_out)
        test = np.load(name)
        print(test['Mvir'][0])

        print("Wrote to {0}.npz".format(fname_out))    

if __name__ == '__main__':

    AllVars.Set_Params_Kali()
    PlotScripts.Set_Params_Plot()
    PlotSnapshot = [[93, 76, 64, 55, 42]]
    #PlotSnapshot = [[50, 64, 76, 93], [50, 64, 76, 93]]
    #model_tags = [r"$f_\mathrm{esc} \: \propto \: \mathrm{SN}$", r"$f_\mathrm{esc} \: \propto \: \mathrm{Quasar}$"]
    model_tags = [r"$f_\mathrm{esc} \: \propto \: \mathrm{Quasar}$", r"$\mathrm{Quasar \: Hot \: Cold}$"]

    #fname = ["/lustre/projects/p004_swin/jseiler/kali/grids/kali_delayedSN_grid256_densfieldgridding_Ejected_alpha0.700beta0.000_fescproperties", "/lustre/projects/p004_swin/jseiler/kali/grids/kali_delayedSN_grid256_densfieldgridding_quasar_0.15_1.00_1.00_fescproperties"]
    #fname = ["/lustre/projects/p004_swin/jseiler/kali/grids/kali_delayedSN_grid256_densfieldgridding_quasar_0.15_1.00_1.00_fescproperties", "/lustre/projects/p004_swin/jseiler/kali/grids/kali_coldhot_quasar_0.15_1.00_1.00_fescproperties"]
    #fname = ["/lustre/projects/p004_swin/jseiler/kali/grids/kali_delayedSN_grid256_densfieldgridding_quasar_0.15_1.00_1.00_fescproperties"]
    fname = ["/lustre/projects/p004_swin/jseiler/kali/grids/kali_fiducial_grid256_fesc0.30_HaloPartCut32_fescproperties"]

    grid_output_to_binary(np.arange(20, 98), fname, "/lustre/projects/p004_swin/jseiler/kali/anne/")
    #plot_fesc_z(fname, PlotSnapshot, model_tags, "grid_fesc_z_test")     


    
   
