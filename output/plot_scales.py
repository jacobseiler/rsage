#!/usr/bin/env python
from __future__ import print_function

import os
import matplotlib
import matplotlib.ticker as mtick
matplotlib.use('Agg')

import pylab as plt
import numpy as np
from numpy.fft import fftn, ifftn

import scipy.integrate as integrate
from scipy import stats

import PlotScripts
import ReadScripts
import AllVars
import misc_func as misc 

label_size = 20
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


    The same again except for just the HII power spectrum.

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

    kmid_bins_XHII, pspec_XHII, p_err_XHII = AllVars.modes_to_pspec(ifftn(XHII),
                                                                    boxsize=AllVars.BoxSize/AllVars.Hubble_h)

    return (kmid_bins, powerspec, p_err,
            kmid_bins_XHII, pspec_XHII, p_err_XHII)


def plot_power(XHII_fraction, P_smallscale, P_largescale, model_tags, 
               OutputDir, mode, plot_mode, target_XHII_fraction=None): 
    """
    Plots the 21cm power spectrum. 

    Parameters
    ---------
    HI_fraction: Numpy array of floats. Required.
        The HI_fraction for each model the power spectrum is calculated for. 

    P_smallscale: Numpy array of floats. Required.
        The value of the powerspectrum at wavenumber k = 0.2h Mpc^-1.

    model_tags: Array of strings with length equal to the number of models. Required.
        Tag that will be added to the legend for each model.

    OutputDir: String. Required.
        The directory the plots will be output to.

    mode: Integer, 0 or 1. Required.
        Denotes whether we're plotting the 21cm Power Spectrum, the HII power
        spectrum.

    plot_mode: Integer, 0 or 1. Required.
        Denotes whether we're plotting Power vs Neutral hydrogen fraction, or
        Large scale power vs Small Scale Power.

    Returns
    -------
    No returns.
    Generates and saves the plot. 

    Units
    -----
    The wavenumber k is in units of h/Mpc.
    The 21cm power spectrum is in units of mK^2. 
    """

    def adjust_plot(ax1, mode, plot_mode, ax2=None):
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

        if plot_mode == 0:
            ax1.set_xlabel(r'$\mathbf{\langle\chi_{HI}\rangle}$', size = label_size)
            axes = [ax1, ax2]    
        else:
            axes = [ax1]

        for ax in axes:

            ax.tick_params(which = 'both', direction='in',
                              width = PlotScripts.global_tickwidth)

            ax.tick_params(which = 'major',
                              length = PlotScripts.global_ticklength)

            ax.tick_params(which = 'minor',
                              length = PlotScripts.global_ticklength - 2.5)

        for axis in ['top','bottom','left','right']:
            ax1.spines[axis].set_linewidth(PlotScripts.global_axiswidth)

        #tick_locs = np.arange(-3, 2.5, 1.0)
        #ax[i].set_xticklabels([r"$\mathbf{10^{%d}}$" % x for x in tick_locs],
                            #fontsize = PlotScripts.global_fontsize)

        if mode == 0:            
            ax1.set_ylabel(r'$\mathbf{\Delta_{21}^2 \left[mK^2\right] (k = 0.2Mpc^{-1}h)}$', size = label_size-2)

            if plot_mode == 0:
                ax2.set_ylabel(r'$\mathbf{\Delta_{21}^2 \left[mK^2\right] (k = 1.0Mpc^{-1}h)}$', size = label_size-2)
            else:
                ax1.set_xlabel(r'$\mathbf{\Delta_{21}^2 \left[mK^2\right] (k = 1.0Mpc^{-1}h)}$', size = label_size-2)
        else:
            ax1.set_ylabel(r'$\mathbf{\Delta_{XHII}^2 \left[mK^2\right] (k = 0.2Mpc^{-1}h)}$', size = label_size-2)
            if plot_mode == 0: 
                ax2.set_ylabel(r'$\mathbf{\Delta_{XHII}^2 \left[mK^2\right] (k = 1.0Mpc^{-1}h)}$', size = label_size-2) 
            else:
                ax1.set_xlabel(r'$\mathbf{\Delta_{XHII}^2 \left[mK^2\right] (k = 1.0Mpc^{-1}h)}$', size = label_size-2)
            
        if mode == 0:   	
            if plot_mode == 0:
                ax1.set_xlim([0.0, 1.00])
                ax1.set_ylim([0.0, 33.0])

                ax2.set_ylim([0.0, 39.0])

                ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.05))
            else:
                ax1.set_xlim([0.0, 39.0])
                ax1.set_ylim([0.0, 33.0])

                ax1.xaxis.set_minor_locator(mtick.MultipleLocator(1))
                ax1.yaxis.set_minor_locator(mtick.MultipleLocator(1))

                print("Ax1 xaxis")
                labels = ax1.xaxis.get_ticklabels()
                locs = ax1.xaxis.get_ticklocs()
                for label, loc in zip(labels, locs):
                    print("{0} {1}".format(label, loc)) 
                
                print("Ax1 yaxis")
                labels = ax1.yaxis.get_ticklabels()
                locs = ax1.yaxis.get_ticklocs()
                for label, loc in zip(labels, locs):
                    print("{0} {1}".format(label, loc)) 

                tick_locs = np.arange(0.0, 45, 10.0)
                ax1.set_yticklabels([r"$\mathbf{%d}$" % x for x in tick_locs],
                                     fontsize = PlotScripts.global_fontsize)

                ax1.set_xticklabels([r"$\mathbf{%d}$" % x for x in tick_locs],
                                     fontsize = PlotScripts.global_fontsize)

                ax1.text(15, 2, r"$\mathbf{Start}$", 
                         fontsize = PlotScripts.global_fontsize-4)
                ax1.text(1, 5, r"$\mathbf{End}$",
                         fontsize = PlotScripts.global_fontsize-4)
        else:       
            ax1.set_ylim([0.0, 0.045])
            ax1.set_xlim([0.0, 0.03])

            if plot_mode == 0:
                ax1.set_xlim([0.0, 1.00])
                ax2.set_ylim([0.0, 0.045])


       
    fig = plt.figure() 
    ax1 = fig.add_subplot(111) 

    if plot_mode == 0:
        ax2 = ax1.twinx()

    for model_number in range(len(XHII_fraction)):

        label = model_tags[model_number]  

        if plot_mode == 0:
            ax1.plot(XHII_fraction[model_number],
                     P_largescale[model_number], 
                     color = PlotScripts.colors[model_number],
                     ls = '-', 
                     lw = PlotScripts.global_linewidth - 2, rasterized=True, 
                     label = label)

            ax2.plot(XHII_fraction[model_number],
                     P_smallscale[model_number], 
                     color = PlotScripts.colors[model_number],
                     ls = '--', 
                     lw = PlotScripts.global_linewidth-2, rasterized=True, 
                     label = label)

        else:
            ax1.plot(P_smallscale[model_number],
                     P_largescale[model_number], 
                     color = PlotScripts.colors[model_number],
                     ls = '-', 
                     lw = PlotScripts.global_linewidth-2, rasterized=True, 
                     label = label)

            for count, val in enumerate(target_XHII_fraction):

                HI_string = "{0:.2f}".format(val)

                idx = (np.abs(XHII_fraction[model_number] - val)).argmin()
                ax1.scatter(P_smallscale[model_number][idx],
                            P_largescale[model_number][idx], 
                            color = PlotScripts.colors[model_number], 
                            rasterized=True, 
                            marker = PlotScripts.markers[count])
                    
    if plot_mode == 0:
        ax1.plot(np.nan, np.nan, ls = '-', label = "Large Scales", 
                 lw = 1, color = 'k')
        ax1.plot(np.nan, np.nan, ls = '--', label = "Small Scales", 
                 lw = 1, color = 'k')

    else:
        one_to_one = np.arange(0.0, 40.0, 1e-6) 
        ax1.plot(one_to_one, one_to_one, ls = '--', lw = 0.5, color = 'k')

        for count, val in enumerate(target_XHII_fraction):
            idx = (np.abs(np.array(XHII_fraction[model_number]) - val)).argmin()

            HI_string = "{0:.2f}".format(val)            
            label = r"$\mathbf{\langle \chi_{HI}\rangle = " +HI_string + r"}$"

            ax1.scatter(np.nan, np.nan, color = 'k',
                        marker = PlotScripts.markers[count],
                        label = label)

    if plot_mode == 0:
        adjust_plot(ax1, mode, plot_mode, ax2)
    else: 
        adjust_plot(ax1, mode, plot_mode) 

    if plot_mode == 0:
        leg = ax1.legend(loc='lower center', numpoints=1,
                           labelspacing=0.1, markerscale=6)
    else:
        leg = ax1.legend(loc='upper left', numpoints=1,
                           labelspacing=0.1, markerscale=6)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize-10)

    plt.tight_layout()

    if mode == 0:
        output_tag = "21cm"
    else:
        output_tag = "XHII"

    if plot_mode == 0:
        output_tag += "_ScalesVHI"
    else: 
        output_tag += "_ScaleVScale"

    outputFile = "{0}/PowerSpec{2}{1}".format(OutputDir,
                                              output_format,
                                              output_tag)

    plt.savefig(outputFile)  # Save the figure

    print('Saved file to {0}'.format(outputFile))

    plt.close()


def load_and_plot(save_tag, target_XHII_fraction, OutputDir):

    fnames = ["XHII_fraction", "P21_smallscale", "P21_largescale",
              "PHII_smallscale", "PHII_largescale", "model_tags"]

    for count, tag in enumerate(fnames):
        fname = "./{0}_{1}.npz".format(save_tag, tag) 

        if count == 0:
            XHII_fraction = np.load(fname)["arr_0"]

        elif count == 1:
            P21_smallscale = np.load(fname)["arr_0"]

        elif count == 2:
            P21_largescale = np.load(fname)["arr_0"]

        elif count == 3:
            PHII_smallscale = np.load(fname)["arr_0"]

        elif count == 4:
            PHII_largescale = np.load(fname)["arr_0"]

        elif count == 5:
            model_tags = np.load(fname)["arr_0"]

    plot_power(XHII_fraction, P21_smallscale, P21_largescale, model_tags, 
               OutputDir, 0, 0)

    plot_power(XHII_fraction, PHII_smallscale, PHII_largescale, model_tags, 
               OutputDir, 1, 0)

    plot_power(XHII_fraction, P21_smallscale, P21_largescale, model_tags, 
               OutputDir, 0, 1, target_XHII_fraction)

    plot_power(XHII_fraction, PHII_smallscale, PHII_largescale, model_tags, 
               OutputDir, 1, 1, target_XHII_fraction)


if __name__ == '__main__':

    PlotScripts.Set_Params_Plot()
    AllVars.Set_Constants()

    filepath_model1="/fred/oz004/jseiler/kali/self_consistent_output/shifted_constant/grids/cifog/new_constant_fesc0.2_XHII"
    filepath_model2="/fred/oz004/jseiler/kali/self_consistent_output/shifted_fej/grids/cifog/shifted_fej_alpha0.4_beta0.0_XHII"
    filepath_model3="/fred/oz004/jseiler/kali/self_consistent_output/shifted_SFR/grids/cifog/shifted_SFR_alpha0.4_beta0.0_XHII"
    filepath_model4="/fred/oz004/jseiler/kali/self_consistent_output/shifted_MHneg/grids/cifog/shifted_MHneg_1e8_1e12_0.99_0.05_XHII"
    filepath_model5="/fred/oz004/jseiler/kali/self_consistent_output/shifted_MHpos/grids/cifog/shifted_MHpos_1e8_1e12_0.01_0.50_XHII"

    fname_ionized=[filepath_model1,
                   filepath_model2,
                   filepath_model3,
                   filepath_model4,
                   filepath_model5]

    fname_density=["/fred/oz004/jseiler/kali/density_fields/1024_subsampled_256/snap",
                   "/fred/oz004/jseiler/kali/density_fields/1024_subsampled_256/snap",
                   "/fred/oz004/jseiler/kali/density_fields/1024_subsampled_256/snap",
                   "/fred/oz004/jseiler/kali/density_fields/1024_subsampled_256/snap",
                   "/fred/oz004/jseiler/kali/density_fields/1024_subsampled_256/snap"]

    precision = [2, 2, 2, 2, 2]

    GridSize = [256, 256, 256, 256, 256]

    model_tags = [r"$\mathbf{f_\mathrm{esc} = 0.20}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: SFR}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}^{-1}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}}$"]

    endsnap = 98 
    SnapList = [np.arange(28, endsnap), np.arange(28, endsnap), 
                np.arange(28, endsnap), np.arange(28, endsnap), 
                np.arange(28, endsnap)]

    cosmo = AllVars.Set_Params_Kali() #  Let's just assume we're always using
                                      #  Kali.

    OutputDir = "./21cm_plots/shifted"
    if not os.path.exists(OutputDir):
        os.makedirs(OutputDir)

    save_tag = "shifted"
    target_HI_fraction = [0.90, 0.75, 0.50, 0.25, 0.10]

    #misc.determine_close_idx(fname_ionized, fname_density, SnapList, GridSize,
    #                         precision, target_HI_fraction, model_tags) 
    #exit()

    XHII_fraction = [[] for x in range(len(SnapList))]

    ##
    have_data = 1
    P21_smallscale = [[] for x in range(len(SnapList))]
    P21_largescale = [[] for x in range(len(SnapList))]

    PHII_smallscale = [[] for x in range(len(SnapList))]
    PHII_largescale = [[] for x in range(len(SnapList))]
    ##

    if have_data == 1:
        load_and_plot(save_tag, target_HI_fraction, OutputDir)
        exit()

    for snapnum in range(len(SnapList[0])):
        print("Snapshot {0}".format(snapnum)) 
        for model_number in range(len(fname_ionized)):
            XHII_fname = "{0}_{1:03d}".format(fname_ionized[model_number], 
                                              SnapList[model_number][snapnum])

            XHII = ReadScripts.read_binary_grid(XHII_fname,
                                                GridSize[model_number],
                                                precision[model_number])

            density_fname = "{0}{1:03d}.dens.dat".format(fname_density[model_number], 
                                                         SnapList[model_number][snapnum]) 
            density = ReadScripts.read_binary_grid(density_fname,
                                                   GridSize[model_number],
                                                   precision[model_number])

            HI_frac = calculate_HI_frac(XHII, density)
            XHII_fraction[model_number].append(HI_frac)
                                               
            if have_data == 0:

                T0 = T_naught(AllVars.SnapZ[SnapList[model_number][snapnum]], 
                              AllVars.Hubble_h, AllVars.Omega_m, AllVars.Omega_b)
                      
                tmp_k, tmp_PowSpec, tmp_Error, \
                tmp_k_XHII, tmp_Pspec_XHII, tmp_Error_XHII = calculate_power_spectrum(XHII, density) 

                tmp_PowSpec *= tmp_k**3 * T0*T0 * 2.0 *np.pi 
                tmp_Pspec_XHII *= tmp_k_XHII**3 * 2.0 *np.pi 

                P21_largescale[model_number].append(tmp_PowSpec[4])
                P21_smallscale[model_number].append(tmp_PowSpec[25])

                PHII_largescale[model_number].append(tmp_Pspec_XHII[4])
                PHII_smallscale[model_number].append(tmp_Pspec_XHII[25])


    fnames = ["XHII_fraction", "P21_smallscale", "P21_largescale",
              "PHII_smallscale", "PHII_largescale", "model_tags"]
    save_arrs = [XHII_fraction, P21_smallscale, P21_largescale,
                 PHII_smallscale, PHII_largescale, model_tags]

    for tag, arr in zip(fnames, save_arrs):
        fname = "./{0}_{1}".format(save_tag, tag) 
        np.savez(fname, arr) 

