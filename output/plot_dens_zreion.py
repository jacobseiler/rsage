#!/usr/bin/env python
from __future__ import print_function

import os
import matplotlib
import matplotlib.ticker as mtick
matplotlib.use('Agg')

import pylab as plt
import numpy as np
from numpy.fft import fftn, ifftn
from matplotlib.colors import LogNorm

import scipy.integrate as integrate
from scipy import stats

import PlotScripts
import ReadScripts
import AllVars
import misc_func as misc 

PlotScripts.Set_Params_Plot()
output_format = ".png"


def plot_zreion_density_hexbin(z_reion, density,
                               output_tag="zreion_density_hexbin"):

    fig1 = plt.figure(figsize=(8,8))

    ax1 = fig1.add_subplot(111)

    w = np.where(density < 5)

    hb = ax1.hexbin(z_reion[w], density[w], cmap='inferno', bins='log')

    cb = fig1.colorbar(hb, ax=ax1)
    cb.set_label('counts')

    outputFile = "./{0}{1}".format(output_tag, output_format) 
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()


def plot_zreion_density_mean(z_reion, density,
                             output_tag="", 
                             density_zreion=True):

    fig1 = plt.figure(figsize=(8,8))

    ax1 = fig1.add_subplot(111)

    if density_zreion:
        density = np.log10(density)
        bins = np.arange(0, 10, 0.1)
        bin_inds = np.digitize(density, bins)

        mean_zreion = np.zeros(len(bins))
        std_zreion = np.zeros(len(bins))
        
    else:
        bins = AllVars.SnapZ[20:97]+0.001
        bin_inds = np.digitize(z_reion, bins)

        mean_density = np.zeros(len(bins))
        std_density = np.zeros(len(bins))

    for bin_idx in range(len(bins)):        
        w = np.where(bin_inds == bin_idx)[0]
        if density_zreion:
            mean_zreion[bin_idx] = np.mean(z_reion[w]) 
            std_zreion[bin_idx] = np.std(z_reion[w]) 
        else:
            mean_density[bin_idx] = np.mean(density[w]) 
            std_density[bin_idx] = np.std(density[w]) 


    if density_zreion:
        ax1.errorbar(bins, mean_zreion, yerr=std_zreion)
        ax1.set_xlabel(r'$\mathbf{log_{10} \delta}$',
                       size = PlotScripts.global_fontsize)

        ax1.set_ylabel(r'$\mathbf{z_{reion}}$',
                       size = PlotScripts.global_fontsize)
        outputFile = "./{0}_density_zreion_mean{1}".format(output_tag, output_format) 
    else:
        ax1.errorbar(AllVars.SnapZ[20:97], mean_density, yerr=std_density)
        ax1.set_ylabel(r'$\mathbf{\delta}$',
                       size = PlotScripts.global_fontsize)
        ax1.set_xlabel(r'$\mathbf{z_{reion}}$',
                       size = PlotScripts.global_fontsize)

        outputFile = "./{0}_zreion_density_mean{1}".format(output_tag, output_format) 
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()



def plot_zreion_density_crosscorr(kbins, crosscorr, model_tags, output_tag):
    
    fig1 = plt.figure(figsize=(8,8))
    ax1 = fig1.add_subplot(111)

    ax1.set_xscale('log')

    for model_number in range(len(kbins)):
        ax1.plot(kbins[model_number][1:-1], crosscorr[model_number][1:-1],
                 label = model_tags[model_number], 
                 color = PlotScripts.colors[model_number])

    ax1.set_xlabel(r"$k [Mpc^{-1}h]$", fontsize = PlotScripts.global_fontsize)
    ax1.set_ylabel(r"$r_{\delta, z_{reion}}$",
                   fontsize = PlotScripts.global_fontsize)

    ax1.set_ylim([0.0, 1.0])

    leg = ax1.legend(loc='upper right', numpoints=1,
                         labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize-2)

    plt.tight_layout()

    outputFile = "./{0}{1}".format(output_tag, output_format) 
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()


def plot_zreion(z_reion_grid, GridSize, output_tag, cut_slice=44,
                cut_thickness=4):

    z_reion_grid = np.reshape(z_reion_grid, (GridSize, GridSize, GridSize))

    fig1 = plt.figure(figsize=(8,8))
    ax1 = fig1.add_subplot(111)

    min_z = min(AllVars.SnapZ[27::])
    max_z = max(AllVars.SnapZ[27::])

    im = ax1.imshow(z_reion_grid[:,:,cut_slice:cut_slice+cut_thickness].mean(axis = -1), 
                    interpolation='nearest', origin='low',
                    vmin=min_z, vmax=max_z, 
                    extent=[0,AllVars.BoxSize,0,AllVars.BoxSize], 
                    cmap='afmhot_r')

    cax = fig1.add_axes()
    cbar = fig1.colorbar(im, cax=cax)
    cbar.ax.set_ylabel(r'$\mathbf{z_{reion}}$', rotation = 90, size = PlotScripts.global_labelsize)
					    
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


def plot_density(density_grid, GridSize, output_tag, cut_slice=127,
                 cut_thickness=3):

    density_grid = np.reshape(density_grid, (GridSize, GridSize, GridSize))

    fig1 = plt.figure(figsize=(8,8))
    ax1 = fig1.add_subplot(111)

    min_dens = np.amin(density_grid)
    max_dens = np.amax(density_grid)
    
    #im = ax1.imshow(np.mean(density_grid[:,cut_slice:cut_slice+cut_thickness,:],axis=1), # Cut xz
    #im = ax1.imshow(np.mean(density_grid[:,:,cut_slice:cut_slice+cut_thickness],axis=2), #Cut xy
    im = ax1.imshow(np.mean(density_grid[cut_slice:cut_slice+cut_thickness,:,:],axis=0), # Cut yz
                    interpolation='nearest', origin='low',
                    norm=LogNorm(vmin=min_dens, vmax=max_dens), 
                    extent=[0,AllVars.BoxSize,0,AllVars.BoxSize], 
                    cmap='afmhot_r')

    cax = fig1.add_axes()
    cbar = fig1.colorbar(im, cax=cax)
    cbar.ax.set_ylabel(r'$\mathbf{\delta}$', rotation = 90, size = PlotScripts.global_labelsize)

    ax1.set_xlabel(r'$\mathrm{y}  (h^{-1}Mpc)$',
                   fontsize=PlotScripts.global_labelsize)  
    ax1.set_ylabel(r'$\mathrm{z}  (h^{-1}Mpc)$',
                   fontsize=PlotScripts.global_labelsize)  

    ax1.set_xlim([0.0, AllVars.BoxSize]) 
    ax1.set_ylim([0.0, AllVars.BoxSize])

    plt.tight_layout()
    #outputFile = OutputDir + output_tag + '_z' + str(z) + output_format 
    outputFile = "./{0}{1}".format(output_tag, output_format) 
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))
            
    plt.close()


def calculate_zreion_dens(zreion_path, dens_path, GridSize, precision):

    kbins = []
    crosscorr = []

    for model_number in range(len(zreion_path)):
        print("Calculating Cross-Corr for model {0}".format(model_number))
        z_reion_grid = ReadScripts.read_binary_grid(zreion_path[model_number],
                                                    GridSize[model_number],
                                                    precision[model_number], 
                                                    True)
 
        density_grid= ReadScripts.read_binary_grid(dens_path[model_number],
                                                   GridSize[model_number],
                                                   precision[model_number], 
                                                   True) 

        density_grid -= 1.0
        z_reion_grid = z_reion_grid/np.mean(z_reion_grid) - 1

        kmid_bins, cross_pspec, pspec1, pspec2 = AllVars.calc_cross_corr(density_grid,
                                                                         z_reion_grid,
                                                                         AllVars.BoxSize)

        crosscoeff = cross_pspec / (pspec1* pspec2)**0.5

        kbins.append(kmid_bins)
        crosscorr.append(crosscoeff)

    return kbins, crosscorr


if __name__ == '__main__':

    AllVars.Set_Params_Kali()

    z_reion_model1="/fred/oz004/jseiler/kali/self_consistent_output/shifted_constant/grids/cifog/new_constant_fesc0.2_reionization_redshift"
    z_reion_model2="/fred/oz004/jseiler/kali/self_consistent_output/shifted_fej/grids/cifog/shifted_fej_alpha0.4_beta0.0_reionization_redshift"
    z_reion_model3="/fred/oz004/jseiler/kali/self_consistent_output/shifted_SFR/grids/cifog/shifted_SFR_alpha0.4_beta0.0_reionization_redshift"
    z_reion_model4="/fred/oz004/jseiler/kali/self_consistent_output/shifted_MHneg/grids/cifog/shifted_MHneg_1e8_1e12_0.99_0.01_reionization_redshift"
    z_reion_model5="/fred/oz004/jseiler/kali/self_consistent_output/shifted_MHpos/grids/cifog/shifted_MHpos_1e8_1e12_0.01_0.50_reionization_redshift"

    zreion_path = [z_reion_model1,
                   z_reion_model2,
                   z_reion_model3,
                   z_reion_model4,
                   z_reion_model5]

    density_path="/fred/oz004/jseiler/kali/density_fields/1024_subsampled_256/snap097.dens.dat"
    dens_path = [density_path,
                 density_path,
                 density_path,
                 density_path,
                 density_path]   
 
    GridSize = [256, 256, 256, 256, 256]
    precision = [2, 2, 2, 2, 2]

    kbins, crosscorr = calculate_zreion_dens(zreion_path, dens_path, 
                                            GridSize, precision)

    model_tags = ["Constant", "fej", "SFR", "MHneg", "MHpos"]
    plot_zreion_density_crosscorr(kbins, crosscorr, model_tags,
                                  "full_zreion_dens_crosscorr")
                                  
