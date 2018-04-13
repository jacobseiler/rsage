#!/usr/bin/env python
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

import ReadScripts
import PlotScripts

import matplotlib.colors as colors
import matplotlib.cm as cm

PlotScripts.Set_Params_Plot()

output_format = ".png"

def my_load_data(filebase, snapshot, numfiles):

    for file_num in range(numfiles):
        fname = "{0}snap{1:03d}_{2}".format(filebase, snapshot, file_num)        
        data = ReadScripts.load_data(fname)
        if file_num == 0:
            array = data 
        else:
            array = np.append(array, data)
    return array 

def plot_hex(x, y, xlabel, ylabel, output_tag):

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    hb = ax1.hexbin(x, y, bins = 'log', cmap = 'inferno')

    ax1.set_xlabel(xlabel, fontsize = PlotScripts.global_labelsize)
    ax1.set_ylabel(ylabel, fontsize = PlotScripts.global_labelsize)

    cb = fig1.colorbar(hb, ax=ax1)
    cb.set_label('log10(N)')

    outputFile1 = "./{0}{1}".format(output_tag, output_format)
    fig1.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig1)


if __name__ == '__main__':

    nfile = 64
    directory = "kali/dust"
    snapshot = 98
    output_suff = "centrals_snap{0:03d}".format(snapshot)

    filebase="/lustre/projects/p004_swin/jseiler/{0}/npz_files/stellarmass_dust_".format(directory)
    stellarmass = my_load_data(filebase, snapshot, nfile)

    print("Successfully loaded array with length {0}".format(len(stellarmass)))

    filebase="/lustre/projects/p004_swin/jseiler/{0}/npz_files/halomass_dust_".format(directory)
    halomass = my_load_data(filebase, snapshot, nfile)

 
    filebase="/lustre/projects/p004_swin/jseiler/{0}/npz_files/dustmass_dust_".format(directory)
    dustmass = my_load_data(filebase, snapshot, nfile)

    xlabel = r'$\mathbf{log_{10} \: M_{vir} \:[M_{\odot}]}$'
    ylabel = r'$\mathbf{log_{10} \: M_{Dust}}$'
    plot_hex(halomass, dustmass, xlabel, ylabel, "Halo_Dust_"+output_suff) 

    xlabel = r'$\mathbf{log_{10} \: M_{vir} \:[M_{\odot}]}$'
    ylabel = r'$\mathbf{log_{10} \: M_{*}}$'
    plot_hex(halomass, stellarmass, xlabel, ylabel,
             "Halo_StellarMass_"+output_suff) 

