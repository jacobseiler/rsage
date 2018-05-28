#!/usr/bin/env python
from __future__ import print_function

import os
import matplotlib
matplotlib.use('Agg')

import pylab as plt
import numpy as np
from numpy.fft import fftn, ifftn

import PlotScripts
import ReadScripts
import AllVars

output_format = ".png"

def T_naught(z, h, OM, OB):
    """
    Calculates the 21cm brightness temperature for specified cosmology +
    redshift. 

    Parameters
    ---------
    z: Float. Required. 
        Redshift we're calculating at.

    h: Float. Required.
        The Hubble h parameter defined as H = 100*h.

    OM: Float. Required.
        Critical matter density. 

    OB: Float. Required. 
        Critical Baryon density.
 
    Returns
    -------
    T0: Float. 
        21cm brightness temperature in units of mK.

    Units
    -----
    All input parameters are unitless.
    Output brightness temperature has units of mK.
    """

    T0 = 28.5 * ((1.0+z)/10.0)**(0.5) * OB/0.042 * h/0.73 * (OM/0.24)**(-0.5)
    return T0


def calculate_HI_frac(XHII, density):
    """
    Calculates the mass-weighted fraction of ionized hydrogen for a given 
    ionization grid.

    Parameters
    ---------
    XHII: 3-Dimensional Array of floats. Required. 
        Grid that contains the fraction of ionized hydrogen (XHII) in each 
        cell. 

    density: 3-Dimensional Array of floats. Required. 
        Grid that contains the overdensity (rho/<rho>) of dark matter in each
        cell. 
 
    Returns
    -------
    HI: Float.
        Fraction of ionized hydrogen. 

    Units
    -----
    XHII and HI are unitless.
    Density is unitless (overdensity, rho/<rho>).    
    """

    HI = 1.0 - np.sum(XHII * density / np.sum(density))
       
    print("")
    print("Mass averaged HI fraction is {0:.4f}".format(HI))
   
    return HI 


def calculate_power_spectrum(XHII, density):
    """
    Calculates the 21cm power spectrum. 

    Parameters
    ---------
    XHII: 3-Dimensional Array of floats. Required. 
        Grid that contains the fraction of ionized hydrogen (XHII) in each 
        cell. 

    density: 3-Dimensional Array of floats. Required. 
        Grid that contains the overdensity (rho/<rho>) of dark matter in each
        cell. 

    Returns
    -------
    kmid_bins: 1-Dimensional array of floats. 
        The middle of each wavenumber bin. 

    powerspec: 1-Dimensional array of floats.
        The 21cm power spectrum in each k bin.

    p_err: 1-Dimensiona array of floats.
        The error on the 21cm power spectrum in each k bin. 

    Units
    -----
    XHII is unitless.
    density is the overdensity, rho/<rho>.
    The wavenumber bins (kmid_bins) is in units of h/Mpc.
    The 21cm power spectrum (and associated error) is in units of Mpc^3/h^3/(2*pi). 
    """

    GridSize = np.shape(XHII)[0]

    meanion = np.mean(XHII)
    Tb = (1.0 - XHII)*density
    modes = ifftn(Tb)

    kmid_bins, powerspec, p_err = AllVars.modes_to_pspec(modes,
                                                         boxsize=AllVars.BoxSize/AllVars.Hubble_h)
    return (kmid_bins, powerspec, p_err)


def plot_power(HI_fraction, k, P21, P21_Error, model_tags, 
               OutputDir, snap):
    """
    Plots the 21cm power spectrum. 

    Parameters
    ---------
    HI_fraction: Numpy array of floats. Required.
        The HI_fraction for each model the power spectrum is calculated for. 

    k: Nested array of floats, k[model_number0] = array of floats. Required.
        The wavenumber (units of h/Mpc) for the power spectrum for each model.

    P21, P21_Error: Nested array of floats, P21[model_number0] = array of floats. Required.
        The dimensionless 21cm power spectrum (units of mK^2) and associated
        error for each model.

    model_tags: Array of strings with length equal to the number of models. Required.
        Tag that will be added to the legend for each model.

    OutputDir: String. Required.
        The directory the plots will be output to.

    snap: Integer. Required. 
        The snapshot number of kali the power spectrum corresponds to.  Used
        for naming purposes.

    Returns
    -------
    No returns.
    Generates and saves the plot (named via output_tag).  

    Units
    -----
    The wavenumber k is in units of h/Mpc.
    The 21cm power spectrum is in units of mK^2. 
    """

    def adjust_21cm_plot(ax):
        """
        Adds text and adjusts the ranges for the 21cm power spectrum. 

        Parameters
        ---------
        ax: Matplotlib axis handle. Required.
           Axis we're adjusting. 
    
        Returns
        -------
        No returns.
        The axis handle is opened by the outer function and passed in here.
        """

        ax.text(0.1, 1.7, r"$\mathrm{Large \: Scales}$", color = 'k', size = PlotScripts.global_fontsize - 2)
        ax.text(1.4, 1.7, r"$\mathrm{Small \: Scales}$", color = 'k', size = PlotScripts.global_fontsize - 2)
        ax.arrow(0.2, 1.65, -0.05, 0.0, head_width = 0.01, head_length = 0.005, fc = 'k', ec = 'k')
        ax.arrow(1.8, 1.65, 0.7, 0.0, head_width = 0.01, head_length = 0.1, fc = 'k', ec = 'k')

        ax.set_xlabel(r'$k \: \left[\mathrm{Mpc}^{-1}h\right]$', 
                       size = PlotScripts.global_labelsize) 
        ax.set_ylabel(r'$\log_{10} \Delta_{21}^2 \left[\mathrm{mK}^2\right]$', 
                       size = PlotScripts.global_labelsize) 

        ax.set_xscale('log')

        leg = ax.legend(loc='lower right', numpoints=1,
             labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(PlotScripts.global_legendsize)

        ax.set_xlim([7e-2, 5.5])
        ax.set_ylim([0.0, 1.8])

    ax1 = plt.subplot(111)

    for model_number in range(len(k)):
        w = np.where(k[model_number] > 9e-2)[0]
       
        label = r"$\langle \chi_\mathrm{HI}\rangle = %.3f$" \
                %(HI_fraction[model_number])

        ax1.plot(k[model_number][w], np.log10(P21[model_number][w]), 
                 ls = PlotScripts.linestyles[model_number], label = label, 
                 lw = 2, rasterized=True) 

    for p in range(0, len(model_tags)):
        ax1.plot(-1, -5, ls = PlotScripts.linestyles[p], label = model_tags[p], color = 'k', lw = 2)

    adjust_21cm_plot(ax1)

    plt.tight_layout()

    outputFile = "{0}/PowerSpec_Snap{1:03d}{2}".format(OutputDir,
                                                       snap,
                                                       output_format)
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))

    plt.close()

if __name__ == '__main__':

    PlotScripts.Set_Params_Plot()

    fname_ionized=["/fred/oz004/jseiler/kali/self_consistent_output/constant/grids/cifog/newphoton_SF0.03_XHII",
                   "/fred/oz004/jseiler/kali/self_consistent_output/quasar/grids/cifog/newphoton_SF0.03_0.25_1.00_2.50_XHII"]

    fname_density=["/fred/oz004/jseiler/kali/density_fields/1024_subsampled_256/snap",
                   "/fred/oz004/jseiler/kali/density_fields/1024_subsampled_256/snap"]
    precision = [2,2]

    GridSize = [256,256]
    model_tags = ["Constant", "Quasar"]
    
    snaplist = np.arange(28, 98)

    cosmo = AllVars.Set_Params_Kali() #  Let's just assume we're always using
                                      #  Kali.

    OutputDir = "./21cm_plots/constant"
    if not os.path.exists(OutputDir):
        os.makedirs(OutputDir)

    for snap in snaplist:
        print("============================")
        print("SNAPSHOT {0}".format(snap))
        print("============================")

        k = []
        P21 = []
        P21_Error = []
        HI_fraction = []

        for model_number in range(len(fname_ionized)):
            XHII_fname = "{0}_{1:03d}".format(fname_ionized[model_number], snap)

            XHII = ReadScripts.read_binary_grid(XHII_fname,
                                                GridSize[model_number],
                                                precision[model_number])

            density_fname = "{0}{1:03d}.dens.dat".format(fname_density[model_number], 
                                                         snap)
            density = ReadScripts.read_binary_grid(density_fname,
                                                   GridSize[model_number],
                                                   precision[model_number])

            HI_frac = calculate_HI_frac(XHII, density)

            T0 = T_naught(AllVars.SnapZ[snap], AllVars.Hubble_h,
                          AllVars.Omega_m, AllVars.Omega_b)
                  
            tmp_k, tmp_PowSpec, tmp_Error = \
            calculate_power_spectrum(XHII, density, GridSize[model_number])

            k.append(tmp_k)
            P21.append(T0*T0 * tmp_PowSpec * tmp_k**3 * 2.0 * np.pi)
            P21_Error.append(tmp_Error)

            HI_fraction.append(HI_frac) 

        plot_power(HI_fraction, k, P21, P21_Error, 
                   model_tags, OutputDir, snap) 

