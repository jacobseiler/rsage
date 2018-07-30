import os
import matplotlib
matplotlib.use('Agg')

import pylab as plt

import numpy as np
from numpy.fft import fftn, ifftn

import ReadScripts
import AllVars

output_format = ".png"


def determine_cross_corr(path1, path2, gridsize1, gridsize2, 
                         precision1, precision2):

    grid1 = ReadScripts.read_binary_grid(path1, gridsize1, precision1)
    grid2 = ReadScripts.read_binary_grid(path2, gridsize2, precision2)
    
    grid1 = np.swapaxes(grid1, 1, 2)
    modes1 = ifftn(grid1)
    modes2 = ifftn(grid2)

    kmid_bins, pspec1, p_err = AllVars.modes_to_pspec(modes1,
                                                      AllVars.BoxSize)

    kmid_bins, pspec2, p_err = AllVars.modes_to_pspec(modes2,
                                                      AllVars.BoxSize)

    kmid_bins, cross_pspec, p_err = AllVars.two_modes_to_pspec(modes1, 
                                                             modes2, 
                                                             AllVars.BoxSize)

    return kmid_bins, cross_pspec, pspec1, pspec2


if __name__ == '__main__':

    AllVars.Set_Params_Kali()

    density_path="/fred/oz004/jseiler/kali/density_fields/1024_subsampled_256/snap097.dens.dat"
    nion_path="/fred/oz004/jseiler/kali/self_consistent_output/new_fej/grids/nion/fej_alpha0.4_beta0.0_ejected_0.400_0.000_HaloPartCut32_nionHI_096"
    XHII_path="/fred/oz004/jseiler/kali/self_consistent_output/new_fej/grids/cifog/fej_alpha0.4_beta0.0_XHII_098"
    output_tag="nion_dens_crosscorr_1_2_swap_crosspspec"

    kmid_bins, cross_pspec, pspec1, pspec2= determine_cross_corr(nion_path, density_path,
                                                                256, 256, 1, 2)
    
    crosscorreff = cross_pspec / (pspec1* pspec2)**0.5

    conversion = (2*np.pi)**3/(2*np.pi**2)

    fig1 = plt.figure(figsize=(8,8))
    ax1 = fig1.add_subplot(111)

    pspec = cross_pspec*kmid_bins**3*conversion 

    #ax1.set_ylim([min(pspec[1:-1]), max(pspec[1:-1])])
    #ax1.set_ylim([min(pspec[1:-1]), max(pspec[1:-1])])
    ax1.set_xscale("log")
    #ax1.set_yscale("symlog")

    ax1.plot(kmid_bins[1:-1], pspec[1:-1]) 
    #ax1.plot(kmid_bins[1:-1], crosscorreff[1:-1]) 

    ax1.set_xlabel(r"$h [Mpc^{-1}h]$")
    #ax1.set_ylabel(r"$r_{\delta, Nion}$")
    ax1.set_ylabel(r"$Cross Power$")

    outputFile = "./{0}{1}".format(output_tag,
                                    output_format)                                    

    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))

    plt.close()

