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

    nfile = 63
    filebase ="/lustre/projects/p004_swin/jseiler/kali/npz_files/stellarmass_dust_"
    stellarmass = my_load_data(filebase, 93, nfile)

    print("Successfully loaded array with length {0}".format(len(stellarmass)))

    filebase ="/lustre/projects/p004_swin/jseiler/kali/npz_files/halomass_dust_"
    halomass = my_load_data(filebase, 93, nfile)

 
    filebase ="/lustre/projects/p004_swin/jseiler/kali/npz_files/dustmass_dust_"
    dustmass = my_load_data(filebase, 93, nfile)

    xlabel = r'$\mathbf{log_{10} \: M_{vir} \:[M_{\odot}]}$'
    ylabel = r'$\mathbf{log_{10} \: M_{Dust}}$'
    plot_hex(halomass, dustmass, xlabel, ylabel, "Halo_Dust") 

    xlabel = r'$\mathbf{log_{10} \: M_{vir} \:[M_{\odot}]}$'
    ylabel = r'$\mathbf{log_{10} \: M_{*}}$'
    plot_hex(halomass, stellarmass, xlabel, ylabel, "Halo_StellarMass") 

