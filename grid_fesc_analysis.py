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

m_gal_low = 5.0 
m_gal_high = 15.0

output_format = ".pdf"

def my_load_data(fname):

    data = ReadScripts.load_data(fname)
    
    snap = data[:,0]
    fesc = data[:,1]
    Mvir = data[:,2] # Units of 1.0e10 Msun/h.
    Mstar = data[:,3] # Units of 1.0e10 Msun/h.
    Ngamma = data[:,4] # Units of 1.0e50 Photons/s.
    Ngammafesc = data[:,5] # Units of 1.0e50 Photons/s.
    

    Mvir = np.log10(Mvir * 1.0e10 / AllVars.Hubble_h)
    Mstar = np.log10(Mstar * 1.0e10 / AllVars.Hubble_h)
    Ngamma = Ngamma * 1.0e50
    Ngammafesc = Ngammafesc * 1.0e50
    
    return snap, fesc, Mvir, Mstar, Ngamma, Ngammafesc

def plot_fescNgamma_z(filename_base, PlotSnapshot, model_tags, output_tag):

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)

    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111) 
 
    for model_number in range(len(filename_base)):        
        total_photons = []
        for count, snapshot_idx in enumerate(PlotSnapshot[model_number]):
            cumulative_nion_ratio_above = []
            cumulative_nion_ratio_below = []
            
            fname = "{0}_{1}".format(filename_base[model_number], snapshot_idx)                            
            snap, fesc, Mvir, Mstar, Ngamma, Ngammafesc = my_load_data(fname)
                                                  
            mean, std, N, sum_, bins_mid = AllVars.Calculate_2D_Mean(Mvir, Ngammafesc, 0.05, m_low, m_high) 

            w = np.where((N > 0.0))
            if model_number == 0:
                label = "z = {0:.2f}".format(AllVars.SnapZ[snapshot_idx])
            else:
                label = ""

            bins_cut = bins_mid[w]
            photons_cut = sum_[w] 
    
            total_photons.append(np.sum(photons_cut))
            
            for bin_idx in range(len(photons_cut)):                
                sum_photons_above = np.sum(photons_cut[bin_idx:-1])
                cumulative_nion_ratio_above.append(sum_photons_above / total_photons[count])

                sum_photons_below = np.sum(photons_cut[0:bin_idx])
                cumulative_nion_ratio_below.append(sum_photons_below / total_photons[count])
        
                if (cumulative_nion_ratio_above[bin_idx] < 0.5):
                    print(bins_cut[bin_idx -1])
                 
            ax1.plot(bins_cut, np.log10(photons_cut), color = PlotScripts.colors[count], ls = PlotScripts.linestyles[model_number], label = label, lw = PlotScripts.global_linewidth)
            ax2.plot(bins_cut, cumulative_nion_ratio_above, color = PlotScripts.colors[count], ls = PlotScripts.linestyles[model_number], label = label, lw = PlotScripts.global_linewidth)
            ax3.plot(bins_cut, photons_cut / total_photons[count], color = PlotScripts.colors[count], ls = PlotScripts.linestyles[model_number], label = label, lw = PlotScripts.global_linewidth)
  
            #print("Total number of ionizing photons for Snapshot {0}, model {1} is {2}".format(snapshot_idx, model_number, total_photons)) 
   
        ax1.plot(np.nan, np.nan, color = 'k', ls = PlotScripts.linestyles[model_number], lw = PlotScripts.global_linewidth, label = model_tags[model_number])
        ax2.plot(np.nan, np.nan, color = 'k', ls = PlotScripts.linestyles[model_number], lw = PlotScripts.global_linewidth, label = model_tags[model_number])

    resolution = np.log10(32.0*AllVars.PartMass / AllVars.Hubble_h)
    #ax1.axvline(resolution, color = 'k', linewidth = PlotScripts.global_linewidth, linestyle = '-.')   
    ax2.axhline(0.5, color = 'k', linewidth = PlotScripts.global_linewidth, linestyle = '-.')   

    ax1.set_xlabel(r'$\log_{10}\ M_\mathrm{vir}\ [M_{\odot}]$', size = PlotScripts.global_fontsize) 
    ax1.set_ylabel(r'$\log_{10} \Sigma \left(\dot{N}_\gamma f_\mathrm{esc}\right) \: [\mathrm{s}^{-1}]$', size = PlotScripts.global_fontsize)
    ax1.set_xlim([resolution, 11.75])
    ax1.set_ylim([53.0, 57.0])   

    ax2.set_xlabel(r'$\log_{10}\ M_\mathrm{vir}\ [M_{\odot}]$', size = PlotScripts.global_fontsize) 
    ax2.set_ylabel(r'$\mathrm{Fraction \: Of} \: \Sigma \dot{N}_{\gamma, \mathrm{esc}} \: \left[\mathrm{M>M}_\mathrm{vir}\right]$', size = PlotScripts.global_fontsize)
    ax2.set_xlim([resolution, 11.75])
    ax2.set_ylim([0.0, 1.0])   

    ax3.set_xlabel(r'$\log_{10}\ M_\mathrm{vir}\ [M_{\odot}]$', size = PlotScripts.global_fontsize) 
    ax3.set_ylabel(r'$\mathrm{Fraction \: Of} \: \Sigma \left(\dot{N}_\gamma f_\mathrm{esc}\right) \: \mathrm{M=M}_\mathrm{vir}$', size = PlotScripts.global_fontsize)
    ax3.set_xlim([resolution, 11.75])
    ax3.set_ylim([0.0, 1.0])   

    ax1.set_xticks(np.arange(9.0, 12.0))  
    ax1.xaxis.set_minor_locator(mtick.MultipleLocator(0.25))
    ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.25))

    ax2.set_xticks(np.arange(9.0, 12.0))  
    ax2.xaxis.set_minor_locator(mtick.MultipleLocator(0.25))
    ax2.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
   
    leg = ax1.legend(loc='lower left', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize)

    leg = ax2.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize)

    leg = ax3.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize)
 
    outputFile1 = './{0}{1}'.format(output_tag, output_format)
    outputFile2 = './{0}_cumratio{1}'.format(output_tag, output_format)
    outputFile3 = './{0}_ratio{1}'.format(output_tag, output_format)
    
    fig.savefig(outputFile1, bbox_inches='tight')
    fig2.savefig(outputFile2, bbox_inches='tight')
    fig3.savefig(outputFile3, bbox_inches='tight')

    print('Saved file to {0}'.format(outputFile1))
    print('Saved file to {0}'.format(outputFile2))
    print('Saved file to {0}'.format(outputFile3))

    plt.close()

   
def plot_fesc_z(filename_base, PlotSnapshot, model_tags, output_tag):

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)

    for model_number in range(len(filename_base)):
        count = 0
        for count, snapshot_idx in enumerate(PlotSnapshot[model_number]):
            fname = "{0}_{1}".format(filename_base[model_number], snapshot_idx)
                            
            snap, fesc, Mvir, Mstar, Ngamma, Ngammafesc = my_load_data(fname)

            print("There are {0} galaxies at Snapshot {1}".format(len(Mvir), snapshot_idx))
                                          
            mean_Mvir, std_Mvir, N_Mvir, sum_fesc_Mvir, bins_mid_Mvir = AllVars.Calculate_2D_Mean(Mvir, fesc, 0.2, m_low, m_high)             
            mean_Mstar, std_Mstar, N_Mstar, sum_fesc_Mstar, bins_mid_Mstar = AllVars.Calculate_2D_Mean(Mstar, fesc, 0.2, m_gal_low, m_gal_high)

            w_halo = np.where((N_Mvir > 0.0))
            w_gal = np.where((N_Mstar > 0.0))
            if model_number == 0:
                label = "z = {0:d}".format(int(round(AllVars.SnapZ[snapshot_idx])))
            else:
                label = ""
            ax1.plot(bins_mid_Mvir[w_halo], mean_Mvir[w_halo], color = PlotScripts.colors[count], ls = PlotScripts.linestyles[model_number], label = label, lw = PlotScripts.global_linewidth)                        
            ax2.plot(bins_mid_Mstar[w_gal], mean_Mstar[w_gal], color = PlotScripts.colors[count], ls = PlotScripts.linestyles[model_number], label = label, lw = PlotScripts.global_linewidth)
            #ax1.fill_between(bins_mid[w], np.subtract(mean[w],std[w]), np.add(mean[w],std[w]), color = PlotScripts.colors[count], alpha = 0.25)

        #ax1.plot(np.nan, np.nan, color = 'k', ls = PlotScripts.linestyles[model_number], lw = PlotScripts.global_linewidth, label = model_tags[model_number])

    resolution = np.log10(32.0*AllVars.PartMass / AllVars.Hubble_h)
    
    ax1.axhline(0.20, 0, 100, color ='k', linewidth = PlotScripts.global_linewidth, linestyle = '--')
    ax1.axhline(0.10, 0, 100, color ='k', linewidth = PlotScripts.global_linewidth, linestyle = '-.')

    ax2.axhline(0.20, 0, 100, color ='k', linewidth = PlotScripts.global_linewidth, linestyle = '--')
    ax2.axhline(0.10, 0, 100, color ='k', linewidth = PlotScripts.global_linewidth, linestyle = '-.')

    #ax1.axvline(resolution, color = 'k', linewidth = PlotScripts.global_linewidth, linestyle = '-.')   
    ax1.text(11.0, 0.22, r"$f_\mathrm{esc} = 0.20$", color = 'k', size = PlotScripts.global_fontsize - 2)
    ax1.text(9.25, 0.11, r"$f_\mathrm{esc, \: base}$", color = 'k', size = PlotScripts.global_fontsize - 2)

    ax2.text(7.6, 0.21, r"$f_\mathrm{esc} = 0.20$", color = 'k', size = PlotScripts.global_fontsize - 2)
    ax2.text(7.6, 0.11, r"$f_\mathrm{esc, \: base}$", color = 'k', size = PlotScripts.global_fontsize - 2)

    ax1.set_xlabel(r'$\log_{10}\ \mathrm{M}_\mathrm{vir}\ [\mathrm{M}_{\odot}]$', size = PlotScripts.global_fontsize) 
    ax1.set_ylabel(r'$\mathrm{Mean} \: f_\mathrm{esc}$', size = PlotScripts.global_fontsize)
    ax1.set_xlim([resolution, 11.75])
    ax1.set_ylim([-0.01, 0.50])

    ax2.set_xlabel(r'$\log_{10}\ \mathrm{M}_\mathrm{*}\ [\mathrm{M}_{\odot}]$', size = PlotScripts.global_fontsize) 
    ax2.set_ylabel(r'$\mathrm{Mean} \: f_\mathrm{esc}$', size = PlotScripts.global_fontsize)
    #ax2.set_xlim([resolution, 11.75])
    ax2.set_ylim([0.05, 0.46])   

    ax1.set_xticks(np.arange(9.0, 12.0))  
    ax1.xaxis.set_minor_locator(mtick.MultipleLocator(0.25))
    ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.05))

    ax2.set_xlim([6.8, 9.8])
    ax2.set_xticks(np.arange(7.0, 11.0))  
    ax2.xaxis.set_minor_locator(mtick.MultipleLocator(0.25))
    ax2.yaxis.set_minor_locator(mtick.MultipleLocator(0.05))
   
    leg = ax1.legend(loc='upper left', bbox_to_anchor=(0.3, 1.02), numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize)

    leg = ax2.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize)
 
    outputFile1 = './{0}_halo{1}'.format(output_tag, output_format)
    outputFile2 = './{0}_galaxy{1}'.format(output_tag, output_format)

    fig.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    fig2.savefig(outputFile2, bbox_inches='tight')  # Save the figure

    print('Saved file to {0}'.format(outputFile1))
    print('Saved file to {0}'.format(outputFile2))
    plt.close()

def plot_nion(base_grid_name, SnapList, GridSize, simulation_norm, plot_time, labels, output_tag):

    print("")
    print("Plotting the total number of ionizing photons") 

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    fig2 = plt.figure()
    ax3 = fig2.add_subplot(111)

    nion = np.empty((len(base_grid_name), len(SnapList[0])))
    dt_nion = np.empty((len(base_grid_name), len(SnapList[0])-1))
    t_dt = np.empty((len(base_grid_name), len(SnapList[0])-1))

    for model_number in range(0, len(base_grid_name)):
        if simulation_norm[model_number] == 0:
            cosmo = AllVars.Set_Params_Mysim()
        elif simulation_norm[model_number] == 2:
            cosmo = AllVars.Set_Params_Tiamat()
        elif simulation_norm[model_number] == 6:
            cosmo = AllVars.Set_Params_Kali()
        else:
            raise ValueError("Set simulation norm option in plot_total_nion")
        for count, snapshot_idx in enumerate(SnapList[model_number]):  
            nion_path = "{0}_{1:03d}".format(base_grid_name[model_number], snapshot_idx) 
            nion_snap = ReadScripts.read_binary_grid(nion_path, GridSize[model_number], 1) * 1.0e50
    
            nion[model_number][count] = np.log10(np.sum(nion_snap) / ((AllVars.BoxSize / AllVars.Hubble_h)**3))	

            if (count != 0):
                dt_nion[model_number][count-1] = (nion[model_number][count] - nion[model_number][count-1]) / (AllVars.SnapZ[SnapList[model_number][count-1]] - AllVars.SnapZ[SnapList[model_number][count]]) 
                t_dt[model_number][count-1] = (AllVars.t_BigBang - AllVars.Lookback_Time[SnapList[model_number][count]])*1.0e3
         
        if plot_time == 1:            
            ax1.plot((AllVars.t_BigBang - AllVars.Lookback_Time[SnapList[model_number]])*1.0e3, nion[model_number], color = PlotScripts.colors[model_number], label = labels[model_number], ls = PlotScripts.linestyles[model_number], lw = 3)
            ax3.plot(t_dt[model_number], dt_nion[model_number], color = PlotScripts.colors[model_number], label = labels[model_number], ls = PlotScripts.linestyles[model_number], lw = 3)
                    
        else:	
            ax1.plot(SnapList[model_number], nion[model_number], color = PlotScripts.colors[model_number], label = labels[model_number], ls = PlotScripts.linestyles[model_number], lw = 3)

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
        ax1.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_labelsize)
        ax1.set_xticks(np.arange(300.0, 1000.0, 100.0))  
        ax1.xaxis.set_minor_locator(mtick.MultipleLocator(PlotScripts.time_tickinterval))
        ax1.set_xlim(PlotScripts.time_xlim)

        ax1.set_ylabel(r'$\sum f_\mathrm{esc}\dot{N}_\gamma \: [\mathrm{s}^{-1}\mathrm{Mpc}^{-3}]$', fontsize = PlotScripts.global_fontsize)                 
        ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.25))
        ax1.set_ylim([48.5, 51.5])

        ax3.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_labelsize)
        ax3.xaxis.set_minor_locator(mtick.MultipleLocator(PlotScripts.time_tickinterval))
        ax3.set_xlim(PlotScripts.time_xlim)

        ax3.set_ylabel(r'$d\sum f_\mathrm{esc}\dot{N}_\gamma / dz \: [\mathrm{s}^{-1}\mathrm{Mpc}^{-3}]$', fontsize = PlotScripts.global_fontsize)                 
        #ax3.yaxis.set_minor_locator(mtick.MultipleLocator(0.25))
        ax3.set_ylim([0.0, 1.0])

        ax2 = ax1.twiny()

        t_plot = (AllVars.t_BigBang - cosmo.lookback_time(PlotScripts.z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
        z_labels = ["$%d$" % x for x in PlotScripts.z_plot] # Properly Latex-ize the labels.

        ax2.set_xlabel(r"$z$", fontsize = PlotScripts.global_labelsize)
        ax2.set_xlim(PlotScripts.time_xlim)
        ax2.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
        ax2.set_xticklabels(z_labels) # But label them as redshifts.

        ax4 = ax3.twiny()

        t_plot = (AllVars.t_BigBang - cosmo.lookback_time(PlotScripts.z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
        z_labels = ["$%d$" % x for x in PlotScripts.z_plot] # Properly Latex-ize the labels.

        ax4.set_xlabel(r"$z$", fontsize = PlotScripts.global_labelsize)
        ax4.set_xlim(PlotScripts.time_xlim)
        ax4.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
        ax4.set_xticklabels(z_labels) # But label them as redshifts.
        
    else:
        ax1.set_xlabel(r"$z$", size = PlotScripts.global_labelsize)
        ax1.set_ylabel(r"$\mathrm{log}_{10} \: \dot{N}_{\mathrm{HI}} \: [\mathrm{s}^{-1}\mathrm{Mpc}^{-3}]$", fontsize = PlotScripts.global_labelsize)
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
        t.set_fontsize(PlotScripts.global_legendsize)

    leg = ax3.legend(loc='lower right', numpoints=1,
                     labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize)

    fig.tight_layout()
    fig2.tight_layout()

    outputFile1 = "./{0}{1}".format(output_tag, output_format) 
    outputFile2 = "./{0}_dtdz{1}".format(output_tag, output_format) 

    fig.savefig(outputFile1)  # Save the figure	
    fig2.savefig(outputFile2)  # Save the figure
	
    print('Saved file to {0}'.format(outputFile1))
    print('Saved file to {0}'.format(outputFile2))
        
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

    #PlotSnapshot = [[35, 40, 45, 50], [35, 40, 45, 50]]
    #PlotSnapshot = [[35, 40, 45, 50]]

    #PlotSnapshot = [[64, 76, 80, 93], [64, 76, 80, 93]]
    #PlotSnapshot = [[64, 76, 80, 93]]
       
    PlotSnapshot = [[42, 64, 76, 93], [42, 64, 76, 93]]
    #PlotSnapshot = [[27, 40, 49, 76, 93]]

    #model_tags = [r"$f_\mathrm{esc} \: \propto \: \mathrm{Quasar, 0.05,1.00,2.00}$", r"$f_\mathrm{esc} = 0.10$"]
    model_tags = [r"$f_\mathrm{esc} \: \propto \: \mathrm{Quasar \: Activity}$", r"$f_\mathrm{esc} = 0.20$"]
    #model_tags = [r"$f_\mathrm{esc} = 0.20$"]        
    #model_tags = [r"$Useless$", r"$f_\mathrm{esc} = 0.15, \: \mathrm{SF 0.03}$"]
    #model_tags = [model_tags[0]]
 
    fname = ["/lustre/projects/p004_swin/jseiler/kali/grids/kali_starburst_quasarwind_SF0.01_NoFractional_QuasarEff0.01_CorrectDiskInstability_quasar_0.10_1.00_2.00_HaloPartCut32_fescproperties", "/lustre/projects/p004_swin/jseiler/kali/grids/kali_starburst_quasarwind_SF0.01_NoFractional_QuasarEff0.01_CorrectDiskInstability_fesc0.20_HaloPartCut32_fescproperties"]     
    fname = [fname[0]]

    #nion_fname = ["/lustre/projects/p004_swin/jseiler/kali/grids/kali_starburst_quasarwind_SF0.01_NoFractional_QuasarEff0.01_CorrectDiskInstability_fesc0.20_HaloPartCut32_nionHI", "/lustre/projects/p004_swin/jseiler/kali/grids/kali_starburst_quasarwind_SF0.01_NoFractional_fesc0.10_HaloPartCut32_nionHI"]
    nion_fname = ["/lustre/projects/p004_swin/jseiler/kali/grids/kali_starburst_quasarwind_SF0.01_NoFractional_QuasarEff0.01_CorrectDiskInstability_quasar_0.10_1.00_2.00_HaloPartCut32_nionHI", "/lustre/projects/p004_swin/jseiler/kali/grids/kali_starburst_quasarwind_SF0.01_NoFractional_QuasarEff0.01_CorrectDiskInstability_fesc0.20_HaloPartCut32_nionHI"] 
    #nion_fname = [nion_fname[0]]   
 
    simulation_norm = [6, 6]

    #plot_nion(nion_fname, [np.arange(27, 99), np.arange(27,99)], [256, 256], simulation_norm, 1, model_tags, "nion_QuasarEff0.01_CorrectDiskInstability_quasar_0.10_1.00_2.00_fesc0.20")

    plot_fesc_z(fname, PlotSnapshot, model_tags, "grid_fesc_QuasarEff0.02_quasar_0.10_1.00_2.00")

    #plot_fescNgamma_z(fname, PlotSnapshot, model_tags, "grid_fesc_z_paper_Ngammafesc_allz_starburst_medium_quasarwind_fesc0.10_fboosted1.00_HaloPartCut32_1DynTime")
