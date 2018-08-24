#!/usr/bin/env python
from __future__ import print_function

import sys
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import pylab as plt
import time
import random
import itertools
import matplotlib.patheffects as PathEffects
from matplotlib import patches

import PlotScripts
import ReadScripts
import AllVars
import misc_func as misc

output_format = ".png"
matplotlib.rcdefaults()

plt.rc('text', usetex=True)


def plot_nine_panel(fname_ionized, SnapList, GridSize, HI_fraction_target,
                    OutputDir, output_tag):

    fig, ax = plt.subplots(nrows=len(HI_fraction_target), ncols=len(SnapList), 
                           sharey=False, sharex=False, figsize=(12, 12))
    fig.subplots_adjust(wspace = 0.02, hspace = 0.02)

    for count, snapnum in enumerate(range(len(SnapList[0]))):
        for model_number in range(len(fname_ionized)):

            this_snapnum = SnapList[model_number][snapnum]
            this_redshift = Redshift[model_number][snapnum]
            this_ax = ax[snapnum][model_number]
 
            XHII_fname = "{0}_{1:03d}".format(fname_ionized[model_number], 
                                              this_snapnum) 

            XHII = ReadScripts.read_binary_grid(XHII_fname,
                                                GridSize[model_number],
                                                precision[model_number])

            ionized_cells = np.log10(1 - XHII)

            index_cut = int(cut_slice * (AllVars.BoxSize/100.0)*(GridSize[model_number]/GridSize[0])) # Wish to cut the box at the same spatial point for all models.  So normalize the index that this corresponds to to model1.
            thickness_cut = int(np.ceil(1 * (AllVars.BoxSize/100.0)*(GridSize[model_number]/GridSize[0]))) # We will take a thickness of 1 cell for model1 and then normalize the number of cells we need to get the same thickness for the other models.

            thickness_cut = 1

            print("For model {0} we begin our cut at index {1} (corresponding to physical position of {2:.4f}) and cut with a thickness of {3} cells.".format(model_number, index_cut, index_cut * AllVars.BoxSize/float(GridSize[model_number]), thickness_cut))
           
            im = this_ax.imshow(ionized_cells[:,:,index_cut:index_cut+thickness_cut].mean(axis = -1), interpolation='none', origin='low', extent =[0,AllVars.BoxSize,0,AllVars.BoxSize], vmin = -8, vmax = 0, cmap = 'afmhot_r')

            this_ax.axis('off')
            
            this_ax.set_xlim([0.0, AllVars.BoxSize]) 
            this_ax.set_ylim([0.0, AllVars.BoxSize])

            #rect = patches.Rectangle((80, 30), 20, 20, lw=3, 
            #                         edgecolor='r', facecolor='none')
            #this_ax.add_patch(rect)

            if (count == 0): 
                title = r"%s" %(model_tags[model_number])
                this_ax.set_title(title, 
                                  size = PlotScripts.global_labelsize - 4) 

            if (model_number == 0): 
                label = (r"$\mathbf{\langle \chi_{HI}\rangle = %.2f}$") %(HI_fraction_target[snapnum])
                this_ax.text(-0.2,0.8, label, transform = this_ax.transAxes, 
                             size = PlotScripts.global_labelsize - 10,
                             rotation=90) 

            tmp = r"$z = %.2f$" %(this_redshift)
            txt = this_ax.text(0.55,0.9, tmp, transform = this_ax.transAxes, 
                          size = PlotScripts.global_labelsize - 16, color = 'k')
            txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
            plt.draw()

    cax = fig.add_axes([0.91, 0.11, 0.03, 0.77])
    ticks = np.arange(-8.0, 1.0, 1.0)
    cbar = fig.colorbar(im, cax=cax, ticks=ticks) 
    cbar.ax.set_yticklabels([r"$\mathbf{%d}$" % x for x in ticks], 
                            fontsize = PlotScripts.global_legendsize+10)
    cbar.ax.set_ylabel(r'$\mathbf{log_{10}\left(\chi_{HI}\right)}$', rotation = 90, size = PlotScripts.global_labelsize)
    #cbar.ax.tick_params(labelsize = PlotScripts.global_legendsize + 10)

    outputFile = "{0}/{1}{2}".format(OutputDir, output_tag, output_format)
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))	
            
    plt.close()


def plot_diff(fname_ionized, SnapList, GridSize, HI_fraction_target,
              OutputDir, output_tag):

    fig, ax = plt.subplots(nrows=len(HI_fraction_target), ncols=len(SnapList), 
                           sharey=False, sharex=False, figsize=(12, 12))
    fig.subplots_adjust(wspace = 0.02, hspace = 0.02)

    for count, snapnum in enumerate(range(len(SnapList[0]))):

        XHII_fname = "{0}_{1:03d}".format(fname_ionized[0], 
                                          SnapList[0][snapnum]) 

        XHII_fid = ReadScripts.read_binary_grid(XHII_fname,
                                                GridSize[0],
                                                precision[0])
       
        XHII_fid[XHII_fid > 0.5] = 1.0

        for model_number in range(len(fname_ionized)):

            this_snapnum = SnapList[model_number][snapnum]
            this_redshift = Redshift[model_number][snapnum]
            this_ax = ax[snapnum][model_number]
 
            XHII_fname = "{0}_{1:03d}".format(fname_ionized[model_number], 
                                              this_snapnum) 

            XHII = ReadScripts.read_binary_grid(XHII_fname,
                                                GridSize[model_number],
                                                precision[model_number])
           
            XHII[XHII > 0.5] = 1.0

            XHII_diff = XHII - XHII_fid
            index_cut = int(cut_slice * (AllVars.BoxSize/100.0)*(GridSize[model_number]/GridSize[0])) # Wish to cut the box at the same spatial point for all models.  So normalize the index that this corresponds to to model1.
            thickness_cut = int(np.ceil(1 * (AllVars.BoxSize/100.0)*(GridSize[model_number]/GridSize[0]))) # We will take a thickness of 1 cell for model1 and then normalize the number of cells we need to get the same thickness for the other models.

            thickness_cut = 1

            print("For model {0} we begin our cut at index {1} (corresponding to physical position of {2:.4f}) and cut with a thickness of {3} cells.".format(model_number, index_cut, index_cut * AllVars.BoxSize/float(GridSize[model_number]), thickness_cut))
           
            im = this_ax.imshow(XHII_diff[:,:,index_cut:index_cut+thickness_cut].mean(axis = -1), 
                                interpolation='none', origin='low', 
                                extent =[0,AllVars.BoxSize,0,AllVars.BoxSize], 
                                vmin = -1, vmax = 1, cmap = 'afmhot_r')

            #this_ax.axis('off')
            
            this_ax.set_xlim([0.0, AllVars.BoxSize]) 
            this_ax.set_ylim([0.0, AllVars.BoxSize])

            if (count == 0): 
                title = r"%s" %(model_tags[model_number])
                this_ax.set_title(title, 
                                  size = PlotScripts.global_labelsize - 4) 

            if (model_number == 0): 
                label = (r"$\mathbf{\langle \chi_{HI}\rangle = %.2f}$") %(HI_fraction_target[snapnum])
                this_ax.text(-0.2,0.8, label, transform = this_ax.transAxes, 
                             size = PlotScripts.global_labelsize - 10,
                             rotation=90) 

            tmp = r"$z = %.2f$" %(this_redshift)
            txt = this_ax.text(0.55,0.9, tmp, transform = this_ax.transAxes, 
                          size = PlotScripts.global_labelsize - 16, color = 'k')
            txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
            plt.draw()

            frame1 = plt.gca()
            frame1.axes.get_xaxis().set_visible(False)
            frame1.axes.get_yaxis().set_visible(False)

            this_ax.tick_params(labelbottom=False)
            this_ax.tick_params(labelleft=False)

    cax = fig.add_axes([0.91, 0.11, 0.03, 0.77])
    ticks = [-1, 0, 1]
    cbar = fig.colorbar(im, cax=cax, ticks=ticks) 
    cbar.ax.set_yticklabels(["Fiducial Ionized",
                             "Both Ionized",
                             "Model Ionized"])
                            
    cbar.ax.set_ylabel(r'', rotation = 90, size = PlotScripts.global_labelsize)    

    outputFile = "{0}/{1}{2}".format(OutputDir, output_tag, output_format)
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))	
            
    plt.close()

if __name__ == '__main__':

    ###########################   

    PlotScripts.Set_Params_Plot()      
 
    model_tags = [r"$\mathbf{f_{esc} = 0.30}$",
                  r"$\mathbf{f_{esc} \: \propto \: M_H^{-1}}$",
                  r"$\mathbf{f_{esc} \: \propto \: M_H}$",
                  r"$\mathbf{f_{esc} \: \propto \: f_{ej}}$",
                  r"$\mathbf{f_{esc} \: \propto \: SFR}$"]

    output_tags = ["Constant",
                   "MH_Neg",
                   "MH_Pos",
                   "Ejected",
                   "SFR"]

    number_models = 5 

    model = 'rsage'

    OutputDir = "./ionization_plots/" + model + '/'
    if not os.path.exists(OutputDir):
        os.makedirs(OutputDir)

    GridSize_model1 = 256
    precision_model1 = 2

    filepath_model1="/fred/oz004/jseiler/kali/self_consistent_output/rsage_constant/grids/cifog/const_0.3_XHII"
    filepath_model2="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHneg/grids/cifog/MHneg_1e8_1e12_0.99_0.05_XHII"
    filepath_model3="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHpos/grids/cifog/MHpos_1e8_1e12_0.01_0.50_XHII"
    filepath_model4="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/grids/cifog/fej_alpha0.50_beta0.0_XHII"
    filepath_model5="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/grids/cifog/SFR_XHII"

    fname_density=["/fred/oz004/jseiler/kali/density_fields/1024_subsampled_256/snap",
                   "/fred/oz004/jseiler/kali/density_fields/1024_subsampled_256/snap",
                   "/fred/oz004/jseiler/kali/density_fields/1024_subsampled_256/snap",
                   "/fred/oz004/jseiler/kali/density_fields/1024_subsampled_256/snap",
                   "/fred/oz004/jseiler/kali/density_fields/1024_subsampled_256/snap"]

    precision = [precision_model1, 
                 precision_model1, 
                 precision_model1, 
                 precision_model1, 
                 precision_model1]

    GridSize = [GridSize_model1, 
                GridSize_model1, 
                GridSize_model1, 
                GridSize_model1, 
                GridSize_model1]

    fname_ionized = [filepath_model1,
                     filepath_model2,
                     filepath_model3,
                     filepath_model4,
                     filepath_model5]

    cosmo = AllVars.Set_Params_Kali() #  Let's just assume we're using Kali.

    SnapList = [np.arange(28, 98),
                np.arange(28, 98),
                np.arange(28, 98),
                np.arange(28, 98),
                np.arange(28, 98)]

    HI_fraction_target = [0.90, 0.75, 0.50, 0.25, 0.10] 

    SnapList = [[18, 25, 33, 41, 47],
                [16, 23, 31, 40, 47],
                [21, 28, 36, 43, 49],
                [17, 23, 31, 40, 46],
                [20, 26, 34, 42, 47]]
    Redshift = np.zeros_like(SnapList, dtype=np.float32)

    cut_slice = 256.0/2.0

    for snap in range(len(SnapList)):
        for inner_snap in range(len(SnapList[snap])):
            SnapList[snap][inner_snap] += 28    
            Redshift[snap][inner_snap] = AllVars.SnapZ[SnapList[snap][inner_snap]]

    plot_nine_panel(fname_ionized, SnapList, GridSize, HI_fraction_target,
                    OutputDir, "nine_panel")

    #plot_diff(fname_ionized, SnapList, GridSize, HI_fraction_target,
    #          OutputDir, "diff")
