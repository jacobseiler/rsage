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

from scipy import stats

label_size = 20
output_format = ".pdf"

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


def plot_power(HI_fraction_target, k, P21, P21_Error, model_tags, 
               OutputDir, mode, output_tag, plot_paper=0): 
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

    mode: Integer, 0 or 1. Required.
        Denotes whether we're plotting the 21cm Power Spectrum or the HII power
        spectrum.

    Returns
    -------
    No returns.
    Generates and saves the plot. 

    Units
    -----
    The wavenumber k is in units of h/Mpc.
    The 21cm power spectrum is in units of mK^2. 
    """

    def adjust_21cm_plot(ax, HI_labels, model_tags, plot_paper): 
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

        if plot_paper == 1:
            for i in range(5):
                ax[i].set_xlabel(r'$\mathbf{k \: \left[Mpc^{-1}h\right]}$', 
                                 size = PlotScripts.global_labelsize)
                ax[i].set_xscale('log')
                ax[i].set_xlim([7e-2, 5.5])

                HI_string = "{0:.2f}".format(HI_labels[i])

                label = r"$\mathbf{\langle \chi_{HI}\rangle = " +HI_string + r"}$"
                ax[i].text(0.05, 0.9, label, transform = ax[i].transAxes,
                           fontsize = PlotScripts.global_fontsize)

                ax[i].tick_params(which = 'both', direction='in',
                                  width = PlotScripts.global_tickwidth)

                ax[i].tick_params(which = 'major',
                                  length = PlotScripts.global_ticklength)

                ax[i].tick_params(which = 'minor',
                                  length = PlotScripts.global_ticklength - 2.5)

                for axis in ['top','bottom','left','right']:
                    ax[i].spines[axis].set_linewidth(PlotScripts.global_axiswidth)

                tick_locs = np.arange(-3, 2.5, 1.0)
                ax[i].set_xticklabels([r"$\mathbf{10^{%d}}$" % x for x in tick_locs],
                                    fontsize = PlotScripts.global_fontsize)


            ax[0].set_ylabel(r'$\mathbf{\Delta_{21}^2 \left[mK^2\right]}$', 
                             size = PlotScripts.global_labelsize)
            ax[0].set_yscale('log', nonposy='clip')
            #ax[0].set_ylim([0.9, 1.8])

            tick_locs = np.arange(-1.0, 5.5, 1.0)
            ax[0].set_yticklabels([r"$\mathbf{10^{%d}}$" % x for x in tick_locs],
                                 fontsize = PlotScripts.global_fontsize)

        else:
            ax.set_xlabel(r'$\mathbf{k \: \left[Mpc^{-1}h\right]}$', 
                          size = PlotScripts.global_labelsize)
            ax.set_xscale('log')
            ax.set_xlim([7e-2, 5.5])

            ax.set_ylabel(r'$\mathbf{\Delta_{21}^2 \left[mK^2\right]}$', 
                          size = PlotScripts.global_labelsize)

    
            HI_string = "{0:.4f}".format(HI_labels)
            label = r"$\mathbf{\langle \chi_{HI}\rangle = " +HI_string + r"}$"

            ax.text(0.05, 0.9, label, transform = ax.transAxes,
                    fontsize = PlotScripts.global_fontsize)

            ax.tick_params(which = 'both', direction='in',
                           width = PlotScripts.global_tickwidth)

            ax.tick_params(which = 'major',
                           length = PlotScripts.global_ticklength)

            ax.tick_params(which = 'minor',
                           length = PlotScripts.global_ticklength - 2.5)

            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(PlotScripts.global_axiswidth)

            tick_locs = np.arange(-3, 2.5, 1.0)
            ax.set_xticklabels([r"$\mathbf{10^{%d}}$" % x for x in tick_locs],
                               fontsize = PlotScripts.global_fontsize)

            ax.set_yscale('log', nonposy='clip')
            ax.set_ylim([1, 100])

            #tick_locs = np.arange(0.0, 5.5, 1.0)
            #ax.set_yticklabels([r"$\mathbf{10^{%d}}$" % x for x in tick_locs],
            #                   fontsize = PlotScripts.global_fontsize)

    def adjust_HII_plot(ax, HI_labels, model_tags):
        """
        Adds text and adjusts the ranges for the HII power spectrum. 

        Parameters
        ---------
        ax: Matplotlib axis handle. Required.
           Axis we're adjusting. 
    
        Returns
        -------
        No returns.
        The axis handle is opened by the outer function and passed in here.
        """

        for i in range(5):
            ax[i].set_xlabel(r'$\mathbf{k \: \left[Mpc^{-1}h\right]}$', size = label_size)
            ax[i].set_xscale('log')
            ax[i].set_xlim([7e-2, 5.5])

            HI_string = "{0:.2f}".format(HI_labels[i])

            label = r"$\mathbf{\langle \chi_{HI}\rangle = " +HI_string + r"}$"
            ax[i].text(0.05, 0.9, label, transform = ax[i].transAxes,
                       fontsize = PlotScripts.global_fontsize)

            ax[i].tick_params(which = 'both', direction='in',
                              width = PlotScripts.global_tickwidth)

            ax[i].tick_params(which = 'major',
                              length = PlotScripts.global_ticklength)

            ax[i].tick_params(which = 'minor',
                              length = PlotScripts.global_ticklength - 2.5)

            for axis in ['top','bottom','left','right']:
                ax[i].spines[axis].set_linewidth(PlotScripts.global_axiswidth)

            tick_locs = np.arange(-3, 2.5, 1.0)
            ax[i].set_xticklabels([r"$\mathbf{10^{%d}}$" % x for x in tick_locs],
                                fontsize = PlotScripts.global_fontsize)


        ax[0].set_ylabel( r'$\mathbf{log_{10} \Delta_{XHII}^2 \left[mK^2\right]}$', size = label_size)
        ax[0].set_yscale('log', nonposy='clip')
        #ax[0].set_ylim([0.9, 1.8])

        #tick_locs = np.arange(0.0, 5.5, 1.0)
        #ax[0].set_yticklabels([r"$\mathbf{10^{%d}}$" % x for x in tick_locs],
        #                     fontsize = PlotScripts.global_fontsize)

        '''
        labels = ax[0].xaxis.get_ticklabels()
        locs = ax[0].xaxis.get_ticklocs()
        for label, loc in zip(labels, locs):
            print("{0} {1}".format(label, loc)) 
        '''

    if plot_paper == 1:
        fig, ax = plt.subplots(nrows=1, ncols=5, sharey='row', figsize=(16, 6))
    else:
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111) 

    log_cutoff = 15 #  This is the transition where we want to plot the mean
                    #  of the log bins.  Represents where our uncertainty
                    #  because basically 0.

    if plot_paper == 0:
        w = np.where(k > 9e-2)[0]        
        ax.plot(k[w],
                P21[w],
                color = "r", 
                ls = "-", 
                lw = 2, rasterized=True, label = "r")
    else:
        for model_number in range(len(k)):
            for fraction in range(len(k[model_number])):

                if plot_paper == 1:
                    this_ax = ax[fraction]
                    this_k = k[model_number]
                    
                else:
                    this_ax = ax[fraction]

                if mode == 0:
                    w = np.where(k[model_number][fraction] > 9e-2)[0]          
                    bins = np.logspace(np.log10(k[model_number][fraction][w[0]]),
                                       np.log10(k[model_number][fraction][w[-1]]),
                                       num = int(len(k[model_number][fraction][w])/1.5))

                    mean_power, bin_edges, bin_number = stats.binned_statistic(k[fraction][model_number][w],
                                                                           np.log10(P21[fraction][model_number][w]),
                                                                           statistic='mean',
                                                                           bins = bins)


                    this_ax.plot(bin_edges[35:-1], pow(10, mean_power[35:]),
                                 color = PlotScripts.colors[model_number],
                                 ls = PlotScripts.linestyles[model_number], 
                                 lw = 2, rasterized=True)

                    label = model_tags[model_number]  

                    this_ax.plot(k[fraction][model_number][w[0]:w[log_cutoff]],
                                 P21[fraction][model_number][w[0]:w[log_cutoff]],
                                 color = PlotScripts.colors[model_number],
                                 ls = PlotScripts.linestyles[model_number],
                                 lw = 2, rasterized=True, label = label)

                else:
                    label = model_tags[model_number]  
                    this_ax.plot(k[fraction][model_number],
                                 P21[fraction][model_number],
                                 color = PlotScripts.colors[model_number],
                                 ls = PlotScripts.linestyles[model_number],
                                 lw = 2, rasterized=True, label = label)
          
    if mode == 0: 
        adjust_21cm_plot(ax, HI_fraction_target, model_tags, plot_paper)                         
    else:
        adjust_HII_plot(ax, HI_fraction_target, model_tags)

    if plot_paper == 0:
        this_ax = ax
    else:
        this_ax = ax[1]

    leg = this_ax.legend(loc='lower right', numpoints=1,
                         labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize-2)

    plt.tight_layout()
    plt.subplots_adjust(wspace = 0.0, hspace = 0.0)


    if mode == 0:
        final_tag = "21cm_{0}".format(output_tag)
    else:
        final_tag = "XHII_{0}".format(output_tag)

    outputFile = "{0}/PowerSpec{2}{1}".format(OutputDir,
                                              output_format,
                                              final_tag)

    plt.savefig(outputFile)  # Save the figure

    print('Saved file to {0}'.format(outputFile))

    plt.close()

if __name__ == '__main__':

    PlotScripts.Set_Params_Plot()

    filepath_model1="/fred/oz004/jseiler/kali/self_consistent_output/constant/grids/cifog/const0.35_XHII"
    filepath_model2="/fred/oz004/jseiler/kali/self_consistent_output/anne/grids/cifog/1e8_1e12_0.99_0.10_XHII"
    filepath_model3="/fred/oz004/jseiler/kali/self_consistent_output/anne/grids/cifog/1e8_1e12_0.01_0.50_XHII"
    filepath_model4="/fred/oz004/jseiler/kali/self_consistent_output/fej/grids/cifog/newphoton_SF0.03_fej0.7_XHII"
    filepath_model5="/fred/oz004/jseiler/kali/self_consistent_output/SFR/grids/cifog/SFR_0.20_0.30_XHII"

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

    model_tags = [r"$\mathbf{f_\mathrm{esc} = 0.35}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}^{-1}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: SFR}$"]               

    SnapList = [np.arange(28, 98),
                np.arange(28, 98),
                np.arange(28, 98),
                np.arange(28, 98),
                np.arange(28, 98)]


    SnapList = [[19, 26, 34, 42, 48],
                [18, 25, 34, 43, 49],
                [22, 29, 37, 45, 52],
                [17, 24, 33, 41, 48],
                [21, 28, 36, 44, 50]]

    for snap in range(len(SnapList)):
        for inner_snap in range(len(SnapList[snap])):
            SnapList[snap][inner_snap] += 28

    HI_fraction_target = [0.90, 0.75, 0.50, 0.25, 0.10] 

    cosmo = AllVars.Set_Params_Kali() #  Let's just assume we're always using
                                      #  Kali.

    OutputDir = "./21cm_plots/paper_final_2"
    if not os.path.exists(OutputDir):
        os.makedirs(OutputDir)

    k = [[] for x in range(len(SnapList))]
    P21 = [[] for x in range(len(SnapList))]
    P21_Error = [[] for x in range(len(SnapList))]


    k_XHII = [[] for x in range(len(SnapList))]
    P21_XHII = [[] for x in range(len(SnapList))]
    P21_Error_XHII = [[] for x in range(len(SnapList))]

    plot_all = 0

    for snapnum in range(len(SnapList[0])): 
        for model_number in range(len(fname_ionized)):
            snap = SnapList[model_number][snapnum]
            XHII_fname = "{0}_{1:03d}".format(fname_ionized[model_number], 
                                              snap) 

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
                  
            tmp_k, tmp_PowSpec, tmp_Error, \
            tmp_k_XHII, tmp_Pspec_XHII, tmp_Error_XHII = calculate_power_spectrum(XHII, density) 

            if plot_all == 0:
                k[snapnum].append(tmp_k)
                P21[snapnum].append(T0*T0 * tmp_PowSpec * tmp_k**3 * 2.0 * np.pi)
                P21_Error[snapnum].append(tmp_Error)

                k_XHII[snapnum].append(tmp_k_XHII)
                P21_XHII[snapnum].append(tmp_Pspec_XHII)
                P21_Error_XHII[snapnum].append(tmp_Error_XHII)

            else:
                output_tag = "{0}_{1}".format(model_tags[model_number], SnapList[model_number][snapnum])
                plot_power(HI_frac, tmp_k, 
                           T0*T0*tmp_PowSpec*tmp_k**3*2.0*np.pi, 
                           tmp_Error, 
                           model_tags, OutputDir, 0, output_tag, 0)

    plot_power(HI_fraction_target, k, P21, 
               P21_Error, model_tags, OutputDir, 0, "5_panel", 1) 

    plot_power(HI_fraction_target, k_XHII, P21_XHII, 
               P21_Error_XHII, model_tags, OutputDir, 1, "5_panel", 1) 


